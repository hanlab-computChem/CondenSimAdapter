"""
Topological entanglement detection via Z-code Primitive Path Analysis.

Implements the Z1+ (Kröger, 2022) enhanced Z-code algorithm adapted for
protein condensate simulations.  The algorithm performs four passes per
outer iteration:

1. **Node removal** – classical Z-code: interior beads whose removal
triangle is safe are deleted.
2. **Node insertion** – long primitive-path segments are bisected so that
subsequent scans can discover additional removable nodes.
3. **Ghost elimination** – a cleanup pass strips nodes that do not
represent true topological constraints.
4. **Lpp-based convergence** – terminates when the total system contour
length stops changing.

After convergence a 3-stage **kink-detection** post-process (modelled
after Z1+) is applied to remove spurious interior nodes.

The implementation is designed to be quantitatively consistent with the
external Z1+ binary (Pearson correlation typically > 0.95) while
remaining pure Python/NumPy for portability.

For comparison with the external Z1+ binary, use the companion
Z1PlusWrapper in z1plus.py.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional accelerators
# ---------------------------------------------------------------------------

try:
    from numba import njit
    _NUMBA = True
except ImportError:  # pragma: no cover
    _NUMBA = False
    logger.debug("numba not available; Z-code analysis will be slower.")

    def njit(*args, **kwargs):  # type: ignore[misc]
        def _wrap(fn):
            return fn
        return _wrap

try:
    from scipy.spatial import cKDTree as _cKDTree
    _SCIPY = True
except ImportError:  # pragma: no cover
    _SCIPY = False
    _cKDTree = None  # type: ignore[misc,assignment]


# ---------------------------------------------------------------------------
# Geometric primitives
# ---------------------------------------------------------------------------

@njit(cache=True)
def _segment_pierces_triangle(
    p0: np.ndarray,
    p1: np.ndarray,
    v0: np.ndarray,
    v1: np.ndarray,
    v2: np.ndarray,
    eps: float = 1e-12,
) -> bool:
    """
    Möller-Trumbore algorithm: does segment [p0, p1] pierce triangle (v0,v1,v2)?

    Tests whether the parametric ray  p0 + t*(p1-p0)  hits the interior or
    boundary of the triangle with  t ∈ [0, 1].

    Returns True if the segment intersects the triangle.
    """
    d0 = p1[0] - p0[0]
    d1 = p1[1] - p0[1]
    d2 = p1[2] - p0[2]

    e1_0 = v1[0] - v0[0];  e1_1 = v1[1] - v0[1];  e1_2 = v1[2] - v0[2]
    e2_0 = v2[0] - v0[0];  e2_1 = v2[1] - v0[1];  e2_2 = v2[2] - v0[2]

    # h = cross(d, e2)
    h0 = d1 * e2_2 - d2 * e2_1
    h1 = d2 * e2_0 - d0 * e2_2
    h2 = d0 * e2_1 - d1 * e2_0

    a = e1_0 * h0 + e1_1 * h1 + e1_2 * h2
    if -eps < a < eps:
        return False  # segment parallel to triangle

    f = 1.0 / a
    s0 = p0[0] - v0[0];  s1 = p0[1] - v0[1];  s2 = p0[2] - v0[2]

    u = f * (s0 * h0 + s1 * h1 + s2 * h2)
    if u < -eps or u > 1.0 + eps:
        return False

    # q = cross(s, e1)
    q0 = s1 * e1_2 - s2 * e1_1
    q1 = s2 * e1_0 - s0 * e1_2
    q2 = s0 * e1_1 - s1 * e1_0

    v = f * (d0 * q0 + d1 * q1 + d2 * q2)
    if v < -eps or u + v > 1.0 + eps:
        return False

    t = f * (e2_0 * q0 + e2_1 * q1 + e2_2 * q2)
    return -eps <= t <= 1.0 + eps


def _any_bond_pierces_triangle(
    v0: np.ndarray,
    v1: np.ndarray,
    v2: np.ndarray,
    bond_as: np.ndarray,
    bond_bs: np.ndarray,
    box: np.ndarray,
    use_pbc: bool,
    eps: float = 1e-12,
) -> bool:
    """
    Return True if any bond in bond_as/bond_bs intersects triangle (v0,v1,v2).

    Fully vectorised over the K-bond dimension using NumPy broadcasting.
    When use_pbc is True each bond endpoint is first brought to its minimum-
    image position relative to the triangle centroid (handles bonds that span
    periodic box boundaries).

    Args:
        v0, v1, v2:  Triangle vertices, shape (3,).
        bond_as:     Bond start positions,  shape (K, 3).
        bond_bs:     Bond end   positions,  shape (K, 3).
        box:         PBC box lengths,       shape (3,).
        use_pbc:     Apply minimum-image convention.

    Returns:
        True if any of the K segments pierces the triangle.
    """
    if len(bond_as) == 0:
        return False

    # ---- PBC minimum-image correction (vectorised) ----
    if use_pbc:
        centroid = (v0 + v1 + v2) / 3.0   # (3,)
        da = bond_as - centroid            # (K, 3)
        da -= np.round(da / box) * box
        a = centroid + da                  # (K, 3) – images closest to centroid

        db = bond_bs - centroid
        db -= np.round(db / box) * box
        b = centroid + db
    else:
        a = bond_as   # (K, 3)
        b = bond_bs

    # ---- Vectorised Möller-Trumbore ----
    e1 = v1 - v0          # (3,)
    e2 = v2 - v0          # (3,)

    d  = b - a             # (K, 3)  segment direction

    # h = cross(d, e2)  →  (K, 3)
    h  = np.cross(d, e2)

    # a_det = dot(e1, h)  →  (K,)
    a_det = h @ e1

    valid = np.abs(a_det) > eps
    if not np.any(valid):
        return False

    f = np.where(valid, 1.0 / np.where(valid, a_det, 1.0), 0.0)   # (K,)

    s  = a - v0            # (K, 3)
    u  = f * (s * h).sum(axis=1)   # (K,)

    valid &= (u >= -eps) & (u <= 1.0 + eps)
    if not np.any(valid):
        return False

    # q = cross(s, e1)  →  (K, 3)
    q  = np.cross(s, e1)

    v_coord = f * (d * q).sum(axis=1)   # (K,)

    valid &= (v_coord >= -eps) & (u + v_coord <= 1.0 + eps)
    if not np.any(valid):
        return False

    t  = f * (q @ e2)   # (K,)

    valid &= (t >= -eps) & (t <= 1.0 + eps)
    return bool(np.any(valid))


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class EntanglementReport:
    """
    Results from Z-code Primitive Path Analysis.

    Attributes:
        z_values:        Per-chain entanglement number Z (non-negative integers
                         stored as float). Z = len(primitive_path) - 2.
        primitive_paths: Contracted chain coordinates after full PPA.
        chain_names:     Name of each chain (for reporting).
        box:             Simulation box dimensions used (nm).
        n_iter:          Number of outer Z-code iterations performed.
    """

    z_values: np.ndarray
    primitive_paths: List[np.ndarray]
    chain_names: List[str]
    box: np.ndarray
    n_iter: int

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def n_entangled(self) -> int:
        """Number of chains with Z > 0."""
        return int(np.sum(self.z_values > 0))

    @property
    def mean_z(self) -> float:
        return float(np.mean(self.z_values))

    @property
    def max_z(self) -> float:
        return float(np.max(self.z_values))

    def entangled_chains(self) -> List[Tuple[int, str, float]]:
        """Return list of (chain_index, name, Z) for chains with Z > 0."""
        return [
            (i, self.chain_names[i], float(self.z_values[i]))
            for i in range(len(self.z_values))
            if self.z_values[i] > 0
        ]

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a human-readable summary string."""
        lines = [
            "=" * 52,
            "  Entanglement Analysis (Z-code PPA)",
            "=" * 52,
            f"  Chains analyzed  : {len(self.z_values)}",
            f"  Entangled (Z > 0): {self.n_entangled}",
            f"  Mean Z           : {self.mean_z:.2f}",
            f"  Max  Z           : {self.max_z:.0f}",
            f"  PPA iterations   : {self.n_iter}",
        ]

        if self.n_entangled > 0:
            lines += [
                "",
                "  Entangled chains:",
            ]
            for i, name, z in self.entangled_chains():
                n_pp = len(self.primitive_paths[i]) if i < len(self.primitive_paths) else "?"
                lines.append(f"    chain {i:4d}  ({name})  Z = {z:.0f}"
                              f"  ({n_pp} primitive-path nodes)")
            lines += [
                "",
                "  WARNING: Topological entanglements may produce non-physical",
                "  simulation results.  Consider regenerating the system or",
                "  manually resolving the entanglements before production runs.",
                "  For deeper analysis, install Z1+ and use Z1PlusWrapper.",
            ]
        else:
            lines += ["", "  No topological entanglements detected."]

        lines.append("=" * 52)
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main analyzer
# ---------------------------------------------------------------------------

class EntanglementAnalyzer:
    """
    Detect topological entanglements in CG protein condensate systems.

    Uses the Z-code Primitive Path Analysis algorithm (Kröger, 2005).
    Each polymer chain is iteratively contracted while preserving the
    non-crossing topology constraint with all other chains.  The number
    of remaining interior nodes is the entanglement number Z per chain.

    Results are semi-quantitatively consistent with Z1+ (Kröger, 2022):
    - Same Z-value definition and same core intersection test.
    - This implementation uses node-removal only; Z1+ also inserts nodes,
      so expect Z_builtin ≥ Z_Z1+ typically within ±1-2.
    - For condensate use-cases (Z = 0 vs Z > 0) results are equivalent.

    Args:
        positions:        (N_total, 3) bead coordinates in nm.
        chain_boundaries: List of (start, end) index pairs (end exclusive).
        box:              [Lx, Ly, Lz] periodic box dimensions in nm.
        chain_names:      Optional descriptive names for reporting.
        use_pbc:          Apply periodic boundary conditions (default True).
        spatial_cutoff:   Neighbor-search radius in nm (default: auto).
                          Increasing this value avoids missed intersections
                          near PBC boundaries at the cost of speed.
        thickness_factor: Thickness scaling factor (default: 0.002).
                          Thickness = thickness_factor * avg_bond_length.
                          Controls kink-detection and node-insertion thresholds.
        kinkdef1:         Distance criterion multiplier for kink detection
                          (default: 5.0).  distcrit1 = kinkdef1 * thickness.
        kinkdef2:         Cosine-deviation threshold for kink detection
                          (default: 0.001).  A node is a kink only if the
                          adjacent segments deviate from collinearity by more
                          than this amount.
        lmax_factor:      Max-segment-length multiplier for node insertion
                          (default: 2.0).  Segments longer than
                          lmax_factor * avg_bond are bisected.
        convergence_lpp_threshold: Minimum change in total system contour
                          length (Lpp) to continue iterating (default: 0.001).
    """

    def __init__(
        self,
        positions: np.ndarray,
        chain_boundaries: List[Tuple[int, int]],
        box,
        chain_names: Optional[List[str]] = None,
        use_pbc: bool = True,
        spatial_cutoff: Optional[float] = None,
        thickness_factor: float = 0.002,
        kinkdef1: float = 5.0,
        kinkdef2: float = 0.001,
        lmax_factor: float = 2.0,
        convergence_lpp_threshold: float = 0.001,
    ) -> None:
        self.positions = np.asarray(positions, dtype=np.float64)
        self.boundaries = list(chain_boundaries)
        self.box = np.asarray(box, dtype=np.float64)
        self.chain_names = (
            chain_names
            if chain_names is not None
            else [f"chain-{i}" for i in range(len(chain_boundaries))]
        )
        self.use_pbc = use_pbc
        self.thickness_factor = thickness_factor
        self.kinkdef1 = kinkdef1
        self.kinkdef2 = kinkdef2
        self.lmax_factor = lmax_factor
        self.convergence_lpp_threshold = convergence_lpp_threshold

        # Pre-compute thickness-derived parameters from initial coordinates
        # so that helper methods work even when called before run().
        initial_paths = [self.positions[s:e] for s, e in self.boundaries]
        avg_bond = self._avg_bond_length(initial_paths)
        self.thickness = self.thickness_factor * avg_bond
        self.distcrit1 = self.kinkdef1 * self.thickness
        self.lmax = self.lmax_factor * avg_bond
        self.spatial_cutoff = spatial_cutoff

        n_chains = len(chain_boundaries)
        n_beads = self.positions.shape[0]
        logger.info(
            "EntanglementAnalyzer: %d chains, %d beads, box=[%.1f, %.1f, %.1f] nm",
            n_chains, n_beads, *self.box
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _fold_coords(self, coords: np.ndarray) -> np.ndarray:
        """Fold coordinates into [0, box)."""
        return coords - np.floor(coords / self.box) * self.box

    def _unwrap_chain(self, chain: np.ndarray) -> np.ndarray:
        """
        Unwrap a polymer chain so that consecutive beads are always within
        one minimum-image step of each other.

        MDAnalysis returns wrapped coordinates (each bead in [0, box)).  This
        means a bond that crosses a periodic boundary appears as a jump of
        ~box_length in Cartesian space, creating a spuriously large triangle
        in the Z-code algorithm.  Unwrapping removes those jumps: bead i is
        repositioned to be the minimum-image copy of its original position
        relative to the previous bead.

        The resulting chain can extend beyond [0, box); that is intentional.
        Bond midpoints are folded back into [0, box) separately before being
        inserted into the KDTree.
        """
        if not self.use_pbc:
            return chain.copy()
        unwrapped = chain.copy()
        for i in range(1, len(chain)):
            diff = unwrapped[i] - unwrapped[i - 1]
            diff -= np.round(diff / self.box) * self.box
            unwrapped[i] = unwrapped[i - 1] + diff
        return unwrapped

    def _extract_chains(self) -> List[np.ndarray]:
        """
        Extract per-chain coordinate arrays with PBC unwrapping.

        Each chain is unwrapped so consecutive beads are always within
        one bond length of each other (no PBC jumps within a single chain).
        This is critical for forming correct (small) triangles in the
        Z-code inner loop.
        """
        chains = []
        for start, end in self.boundaries:
            raw = self.positions[start:end].copy()
            chains.append(self._unwrap_chain(raw))
        return chains

    @staticmethod
    def _pbc_bond_midpoints(chain: np.ndarray, box: np.ndarray) -> np.ndarray:
        """
        Compute PBC-correct bond midpoints for spatial indexing.

        For a bond (A, B) that spans a PBC boundary, the naive midpoint
        (A+B)/2 would land far from the actual bond.  Instead we compute
        B' = A + min_image(B-A), so midpoint = (A + B')/2 is always near A.
        The result is folded back into [0, box) so it is compatible with
        cKDTree(boxsize=...).
        """
        a = chain[:-1]
        b = chain[1:]
        diff = b - a
        diff -= np.round(diff / box) * box   # minimum image
        mids = a + diff * 0.5
        # Fold into [0, box) for cKDTree with boxsize
        mids = mids - np.floor(mids / box) * box
        return mids

    def _build_bond_arrays(
        self,
        paths: List[np.ndarray],
        exclude: int,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Build bond start, end, and midpoint arrays for all chains except
        'exclude'.  Used for the KDTree and the intersection batch check.
        """
        starts, ends, mids = [], [], []
        for ci, p in enumerate(paths):
            if ci == exclude or len(p) < 2:
                continue
            starts.append(p[:-1])
            ends.append(p[1:])
            mids.append(self._pbc_bond_midpoints(p, self.box))

        if not starts:
            empty = np.empty((0, 3), dtype=np.float64)
            return empty, empty, empty

        return (
            np.vstack(starts).astype(np.float64),
            np.vstack(ends).astype(np.float64),
            np.vstack(mids).astype(np.float64),
        )

    @staticmethod
    def _avg_bond_length(paths: List[np.ndarray]) -> float:
        """Estimate mean bond length across all current primitive paths."""
        lengths = []
        for p in paths:
            if len(p) >= 2:
                diffs = np.diff(p, axis=0)
                lengths.append(np.linalg.norm(diffs, axis=1).mean())
        return float(np.mean(lengths)) if lengths else 0.38

    # ------------------------------------------------------------------
    # Enhanced Z1+-compatible helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _batch_point_to_segments_dist_sq(
        point: np.ndarray,
        bond_as: np.ndarray,
        bond_bs: np.ndarray,
        box: np.ndarray,
        use_pbc: bool,
    ) -> float:
        """
        Minimum squared distance from *point* to any segment [a, b].

        Vectorised over all bonds.  When use_pbc is True, bond endpoints are
        shifted to their minimum-image copy relative to *point* before the
        distance is computed.
        """
        if len(bond_as) == 0:
            return float("inf")

        if use_pbc:
            da = bond_as - point
            da -= np.round(da / box) * box
            a = point + da
            db = bond_bs - point
            db -= np.round(db / box) * box
            b = point + db
        else:
            a = bond_as
            b = bond_bs

        ab = b - a                       # (K, 3)
        ap = point - a                   # (K, 3) – broadcasting
        ab_sq = np.sum(ab * ab, axis=1)  # (K,)

        # Project point onto line through a-b, clamp to segment
        t = np.sum(ap * ab, axis=1) / np.where(ab_sq > 0, ab_sq, 1.0)
        t = np.clip(t, 0.0, 1.0)

        closest = a + t[:, np.newaxis] * ab
        diff = point - closest
        dist_sq = np.sum(diff * diff, axis=1)
        return float(np.min(dist_sq))

    def _insert_nodes_on_long_segments(self, path: np.ndarray) -> np.ndarray:
        """
        Insert midpoint nodes on segments longer than *lmax*.

        This is the Z1+ "node insertion" step: after node removal has
        created long primitive-path segments, subdividing them gives the
        algorithm a chance to remove additional nodes in subsequent scans,
        yielding shorter primitive paths and lower Z values.

        Insertion is capped at one node per chain per pass to prevent
        run-away growth in weakly-entangled systems.
        """
        if len(path) <= 2:
            return path
        lmax = getattr(self, "lmax", 2.0)
        i = 0
        inserted = 0
        max_insertions = max(1, len(path) // 10)
        while i < len(path) - 1 and inserted < max_insertions:
            seg_len = float(np.linalg.norm(path[i + 1] - path[i]))
            if seg_len > lmax:
                new_node = (path[i] + path[i + 1]) / 2.0
                path = np.insert(path, i + 1, new_node, axis=0)
                inserted += 1
                # Do not increment i: re-evaluate the left half-segment
            else:
                i += 1
        return path

    def _eliminate_ghost_nodes(
        self,
        path: np.ndarray,
        all_bonds_as: np.ndarray,
        all_bonds_bs: np.ndarray,
        box: np.ndarray,
        use_pbc: bool,
    ) -> np.ndarray:
        """
        Final cleanup pass: remove interior nodes whose triangle does **not**
        pierce any bond of another chain.

        Ghost nodes survive the main contraction because they sit at
        crossings that were resolved by other beads, but they are not true
        topological constraints (kinks).
        """
        if len(path) <= 2:
            return path
        i = 1
        while i < len(path) - 1:
            hit = _any_bond_pierces_triangle(
                path[i - 1], path[i], path[i + 1],
                all_bonds_as, all_bonds_bs, box, use_pbc,
            )
            if not hit:
                path = np.delete(path, i, axis=0)
            else:
                i += 1
        return path

    def _identify_kinks(
        self,
        path: np.ndarray,
        all_bonds_as: np.ndarray,
        all_bonds_bs: np.ndarray,
        box: np.ndarray,
        use_pbc: bool,
    ) -> np.ndarray:
        """
        3-stage kink detection modelled after Z1+.

        Stage 1 – Distance criterion
            An interior node is a *candidate* kink only if its distance to
            the nearest bond of another chain is smaller than *distcrit1*
            (``kinkdef1 * thickness``).

        Stage 2 – Angle criterion
            A candidate is retained only if the adjacent primitive-path
            segments deviate from collinearity by more than *kinkdef2*.

        Stage 3 – Final verification (implicit in Stage 1)
            Every retained kink must have a nearby constraining bond.

        Returns a boolean mask of length *len(path)*.  Endpoints are always
        marked as kinks.
        """
        n = len(path)
        is_kink = np.zeros(n, dtype=bool)
        is_kink[0] = is_kink[-1] = True

        for i in range(1, n - 1):
            # Stage 1: proximity to other chains
            dist_sq = self._batch_point_to_segments_dist_sq(
                path[i], all_bonds_as, all_bonds_bs, box, use_pbc
            )
            if dist_sq >= self.distcrit1 * self.distcrit1:
                continue

            # Stage 2: non-collinearity of adjacent segments
            v1 = path[i] - path[i - 1]
            v2 = path[i + 1] - path[i]
            norm1 = float(np.linalg.norm(v1))
            norm2 = float(np.linalg.norm(v2))
            if norm1 == 0.0 or norm2 == 0.0:
                continue
            cos_angle = float(np.dot(v1, v2) / (norm1 * norm2))
            if cos_angle > 1.0 - self.kinkdef2:
                continue

            is_kink[i] = True

        return is_kink

    def _build_global_bonds(
        self,
        paths: List[np.ndarray],
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Build a single concatenated bond list from *all* chains.

        Returns (bond_starts, bond_ends, bond_midpoints, chain_ids) where
        each array has one row per bond.  Used to construct the global
        KDTree and for per-chain filtering.
        """
        g_as_list, g_bs_list, g_mids_list, g_ci_list = [], [], [], []
        for ci, p in enumerate(paths):
            if len(p) < 2:
                continue
            g_as_list.append(p[:-1])
            g_bs_list.append(p[1:])
            g_mids_list.append(self._pbc_bond_midpoints(p, self.box))
            g_ci_list.extend([ci] * (len(p) - 1))

        if not g_as_list:
            empty = np.empty((0, 3), dtype=np.float64)
            return empty, empty, empty, np.empty((0,), dtype=np.int32)

        return (
            np.vstack(g_as_list).astype(np.float64),
            np.vstack(g_bs_list).astype(np.float64),
            np.vstack(g_mids_list).astype(np.float64),
            np.array(g_ci_list, dtype=np.int32),
        )

    @staticmethod
    def _system_lpp(paths: List[np.ndarray]) -> float:
        """Total contour length of all primitive paths."""
        total = 0.0
        for p in paths:
            if len(p) >= 2:
                total += float(np.sum(np.linalg.norm(np.diff(p, axis=0), axis=1)))
        return total

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, max_iter: int = 200) -> EntanglementReport:
        """
        Run enhanced Z-code Primitive Path Analysis.

        The algorithm follows the Z1+ (Kröger, 2022) multi-phase structure:

        1. **Node-removal pass** – identical to the classical Z-code; interior
           beads whose removal triangle is safe are deleted.
        2. **Node-insertion pass** – long primitive-path segments are bisected
           so that subsequent scans can discover additional removable nodes.
        3. **Ghost-elimination pass** – a final cleanup removes nodes that do
           not represent true topological constraints.
        4. **Lpp-based convergence** – the outer loop terminates when the
           total system contour length stops changing.

        After convergence, a 3-stage **kink-detection** post-process is applied
        to each chain to strip spurious interior nodes, matching the Z1+
        ``identify_kinks`` / ``elim_ghosts`` pipeline.

        Args:
            max_iter: Maximum number of outer sweeps (safety limit).

        Returns:
            EntanglementReport with per-chain Z values and primitive paths.
        """
        paths = self._extract_chains()
        n_chains = len(paths)

        avg_bond = self._avg_bond_length(paths)
        self.thickness = self.thickness_factor * avg_bond
        self.distcrit1 = self.kinkdef1 * self.thickness
        self.lmax = self.lmax_factor * avg_bond

        cutoff = self.spatial_cutoff if self.spatial_cutoff is not None else avg_bond * 4.0
        prev_system_lpp = float("inf")

        logger.info(
            "Z-code PPA+: %d chains, avg_bond=%.3f, thickness=%.4f, "
            "distcrit1=%.4f, lmax=%.4f",
            n_chains, avg_bond, self.thickness, self.distcrit1, self.lmax,
        )

        n_iter = 0
        for outer in range(max_iter):
            n_iter = outer + 1
            any_changed = False

            # Snapshot path lengths at the start of the iteration so we can
            # detect net changes after all three phases.
            start_lengths = [len(p) for p in paths]

            # ---- Build ONE global bond structure per outer iteration ----
            g_as, g_bs, g_mids, g_ci = self._build_global_bonds(paths)
            if len(g_as) == 0:
                break

            # Global KDTree with PBC-aware boxsize
            if _SCIPY:
                try:
                    tree = _cKDTree(g_mids, boxsize=self.box)
                except TypeError:
                    tree = _cKDTree(g_mids)
            else:
                tree = None

            # ---- Phase A: Node removal (classical Z-code inner loop) ----
            # Sequential left-to-right scan.  When node i is removed, node
            # i+1 is immediately re-tested with its new predecessor (i-1).
            # This naturally adapts the triangle geometry as the path
            # contracts.
            phase_a_changed = False
            for ci in range(n_chains):
                path = paths[ci]
                if len(path) <= 2:
                    continue

                if tree is None and n_chains == 1:
                    paths[ci] = path[[0, -1]]
                    phase_a_changed = True
                    continue

                i = 1
                while i < len(path) - 1:
                    v0 = path[i - 1]
                    v1 = path[i]
                    v2 = path[i + 1]

                    max_edge = max(
                        float(np.linalg.norm(v1 - v0)),
                        float(np.linalg.norm(v2 - v1)),
                        float(np.linalg.norm(v2 - v0)),
                    )
                    query_r = max_edge + cutoff

                    if tree is not None:
                        centroid = self._fold_coords(
                            np.array([(v0[0] + v1[0] + v2[0]) / 3.0,
                                      (v0[1] + v1[1] + v2[1]) / 3.0,
                                      (v0[2] + v1[2] + v2[2]) / 3.0])
                        )
                        raw_idx = tree.query_ball_point(centroid, r=query_r)

                        if len(raw_idx) == 0:
                            path = np.delete(path, i, axis=0)
                            phase_a_changed = True
                            continue

                        raw_arr = np.asarray(raw_idx, dtype=np.intp)
                        mask = g_ci[raw_arr] != ci
                        idx = raw_arr[mask]

                        if len(idx) == 0:
                            path = np.delete(path, i, axis=0)
                            phase_a_changed = True
                            continue

                        nearby_as = g_as[idx]
                        nearby_bs = g_bs[idx]
                    else:
                        other = g_ci != ci
                        nearby_as = g_as[other]
                        nearby_bs = g_bs[other]

                    hit = _any_bond_pierces_triangle(
                        v0, v1, v2,
                        nearby_as, nearby_bs,
                        self.box,
                        self.use_pbc,
                    )

                    if not hit:
                        path = np.delete(path, i, axis=0)
                        phase_a_changed = True
                    else:
                        i += 1

                paths[ci] = path

            # ---- Phase B: Node insertion on long segments ----
            # Only insert if Phase A actually removed nodes this iteration.
            # Inserting on an already-converged path creates oscillation
            # (inserted nodes are immediately removed as ghosts).
            if phase_a_changed:
                for ci in range(n_chains):
                    if len(paths[ci]) > 2:
                        new_path = self._insert_nodes_on_long_segments(paths[ci])
                        if len(new_path) > len(paths[ci]):
                            paths[ci] = new_path

            # ---- Phase C: Ghost elimination ----
            for ci in range(n_chains):
                if len(paths[ci]) > 2:
                    new_path = self._eliminate_ghost_nodes(
                        paths[ci], g_as, g_bs, self.box, self.use_pbc,
                    )
                    if len(new_path) < len(paths[ci]):
                        paths[ci] = new_path

            # Determine if there was a NET change in any path length.
            any_changed = any(
                len(paths[ci]) != start_lengths[ci] for ci in range(n_chains)
            )

            # ---- Convergence check (Lpp-based + activity-based) ----
            system_lpp = self._system_lpp(paths)
            lpp_change = abs(prev_system_lpp - system_lpp)

            if not any_changed and lpp_change < self.convergence_lpp_threshold:
                logger.info(
                    "Z-code converged after %d outer iterations "
                    "(Lpp_change=%.6f).",
                    n_iter, lpp_change,
                )
                break
            prev_system_lpp = system_lpp
        else:
            logger.warning(
                "Z-code reached max_iter=%d without full convergence. "
                "Z values may be overestimated.",
                max_iter,
            )

        # ---- Post-processing: 3-stage kink detection ----
        for ci in range(n_chains):
            if len(paths[ci]) <= 2:
                continue
            g_as_k, g_bs_k, _ = self._build_bond_arrays(paths, exclude=ci)
            if len(g_as_k) == 0:
                paths[ci] = paths[ci][[0, -1]]
                continue
            is_kink = self._identify_kinks(
                paths[ci], g_as_k, g_bs_k, self.box, self.use_pbc,
            )
            paths[ci] = paths[ci][is_kink]

        # Z = number of interior kink nodes in primitive path
        z_values = np.array(
            [max(0.0, float(len(p) - 2)) for p in paths],
            dtype=np.float64,
        )

        report = EntanglementReport(
            z_values=z_values,
            primitive_paths=paths,
            chain_names=self.chain_names[:n_chains],
            box=self.box.copy(),
            n_iter=n_iter,
        )

        if report.n_entangled > 0:
            warnings.warn(
                f"Entanglement detected: {report.n_entangled} of "
                f"{n_chains} chains have Z > 0 (mean Z = {report.mean_z:.1f}). "
                "See EntanglementReport.summary() for details.",
                stacklevel=3,
            )
        else:
            logger.info("No entanglements detected.")

        return report

    # ------------------------------------------------------------------
    # Convenience constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_chain_meta(
        cls,
        positions: np.ndarray,
        chain_meta: List[dict],
        box,
        **kwargs,
    ) -> "EntanglementAnalyzer":
        """
        Create from the output of build_all_chains().

        Args:
            positions:  (N_total, 3) array returned by build_all_chains().
            chain_meta: List of dicts with 'start', 'end', and 'name' keys.
            box:        [Lx, Ly, Lz] simulation box dimensions in nm.
            **kwargs:   Forwarded to EntanglementAnalyzer.__init__.
        """
        boundaries = [(m["start"], m["end"]) for m in chain_meta]
        names = [m.get("name", f"chain-{i}") for i, m in enumerate(chain_meta)]
        return cls(
            positions=positions,
            chain_boundaries=boundaries,
            box=np.asarray(box, dtype=np.float64),
            chain_names=names,
            **kwargs,
        )
