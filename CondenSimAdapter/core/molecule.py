"""
Coordinate builder for CG protein systems.

Handles:
  - Single-chain initial coordinate generation (IDP: spiral/linear; MDP: from PDB)
  - Multi-chain placement (grid, slab, droplet/random)
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import List, Optional, Tuple

from .config import Component, ComponentType, CGConfig, TopologyType

# Minimum safe distance (nm) between any compact-IDR atom and any folded-domain atom
# when building the initial MDP chain coordinates.  Chosen to exceed typical WF sigma
# (~0.47 nm for Gly, ~0.57 nm for Trp) so no WF repulsion spike occurs for any
# intra-chain IDR-folded pair that is not in the 1-2 backbone exclusion list.
_MIN_IDR_CLEARANCE = 0.55


# ---------------------------------------------------------------------------
# Single-chain coordinate generation
# ---------------------------------------------------------------------------

def build_idp_chain(
    n_beads: int,
    method: str = "spiral",
    spacing: float = 0.38,
) -> np.ndarray:
    """
    Generate CA coordinates for a fully disordered chain (IDP).

    Args:
        n_beads:  Number of residues / CA beads.
        method:   'spiral' | 'linear' | 'compact'
        spacing:  Bond length in nm (default 0.38 nm for proteins).

    Returns:
        (N, 3) float64 array in nm, centered at origin.
    """
    if method == "spiral":
        return _build_spiral(n_beads, spacing)
    elif method == "linear":
        return _build_linear(n_beads, spacing)
    elif method == "compact":
        return _build_compact(n_beads, spacing)
    else:
        raise ValueError(f"Unknown IDP build method: {method!r}")


def build_mdp_chain(
    pdb_path: str,
    use_com: bool = False,
) -> Tuple[np.ndarray, str]:
    """
    Extract CG (CA or residue-COM) coordinates from an all-atom or CA PDB.
    Mirrors the behavior of CALVADOS build.geometry_from_pdb.

    Args:
        pdb_path: Path to PDB file.
        use_com:  If True, use residue centre-of-mass instead of CA.

    Returns:
        coords: (N, 3) float64 array in nm, mass-centred at origin.
        sequence: one-letter amino acid string.
    """
    import MDAnalysis as mda
    from MDAnalysis.lib.util import convert_aa_code
    from warnings import catch_warnings, simplefilter

    # Suppress MDAnalysis warnings (matches original CALVADOS behavior)
    with catch_warnings():
        simplefilter("ignore")
        u = mda.Universe(pdb_path)
    
    residues = u.select_atoms("protein").residues

    coords = []
    seq = ""
    for res in residues:
        if use_com:
            xyz = res.atoms.center_of_mass() / 10.0   # Å -> nm
        else:
            ca = res.atoms.select_atoms("name CA")
            if len(ca) == 0:
                ca = res.atoms
            xyz = ca.positions[0] / 10.0
        coords.append(xyz)
        try:
            seq += convert_aa_code(res.resname)
        except Exception:
            seq += "G"

    coords = np.array(coords, dtype=np.float64)
    coords -= coords.mean(axis=0)   # centre at origin
    return coords, seq


def build_mdp_chain_compact_idr(
    pdb_path: str,
    folded_domains: List[Tuple[int, int]],
    n_res: int,
    use_com: bool = False,
    spacing: float = 0.38,
) -> np.ndarray:
    """
    Build MDP initial coordinates with compact IDR regions.

    The folded domain residues are taken directly from the PDB CA positions.
    IDR (disordered) regions are replaced by compact cubic builds anchored near
    the folded domain termini, avoiding the extended IDR conformations that a
    PDB from an MD run may carry (can span > 10 nm) and would cause severe
    inter-chain clashes on the initial slab grid.

    Args:
        pdb_path:       Path to PDB file (all-atom or CA-only).
        folded_domains: List of 1-based inclusive (start, end) domain ranges.
        n_res:          Expected number of residues (from sequence).
        use_com:        If True, use residue centre-of-mass instead of CA.
        spacing:        Bond spacing in nm (default 0.38 nm).

    Returns:
        (n_res, 3) float64 array in nm, centred at origin.
    """
    coords_pdb, _ = build_mdp_chain(pdb_path, use_com=use_com)

    if len(coords_pdb) != n_res or not folded_domains:
        # Residue count mismatch or no domain info: fall back to raw PDB
        return coords_pdb

    coords = coords_pdb.copy()

    # Collect 0-based folded residue indices
    folded_set: set = set()
    for dom_s, dom_e in folded_domains:
        folded_set.update(range(dom_s - 1, dom_e))

    # Folded domain centroid (used to determine outward direction)
    folded_coords = coords_pdb[sorted(folded_set)]
    folded_center = folded_coords.mean(axis=0)

    for seg_start, seg_end in _get_idr_segments(n_res, folded_domains):
        n_seg = seg_end - seg_start
        if n_seg <= 0:
            continue

        # Number of layers in compact cube (same formula as _build_compact)
        N = max(1, int(np.ceil(np.cbrt(n_seg)) - 1))

        # Build compact structure centred at origin
        idr_c = _build_compact(n_seg, spacing)
        idr_c -= idr_c.mean(axis=0)

        # Determine the anchor (folded terminus adjacent to this IDR segment)
        # and the outward direction (away from the folded domain center).
        if seg_start > 0 and (seg_start - 1) in folded_set:
            # C-terminal IDR (or inter-domain IDR after a folded block)
            anchor = coords_pdb[seg_start - 1]
        else:
            # N-terminal IDR (before first folded block)
            anchor = coords_pdb[seg_end]

        outward = anchor - folded_center
        norm = np.linalg.norm(outward)
        outward = outward / norm if norm > 1e-6 else np.array([0.0, 0.0, 1.0])

        # Analytic lower bound for the clearance needed:
        # The compact cube can have an atom at up to (N/2 * sqrt(3) * spacing)
        # in the backward direction from the centroid (cube body diagonal worst case).
        # We start the centroid at anchor + outward * d_start and then shift it
        # outward iteratively until every IDR atom is >= _MIN_IDR_CLEARANCE nm from
        # every folded-domain atom (like CALVADOS check_clash but for intra-chain).
        max_backward = (N / 2.0) * np.sqrt(3.0) * spacing
        d_start = max_backward + spacing   # initial clearance guess

        idr_candidate = idr_c + (anchor + outward * d_start)
        for _ in range(200):
            dists = np.sqrt(
                ((idr_candidate[:, None, :] - folded_coords[None, :, :]) ** 2).sum(axis=-1)
            )
            if dists.min() >= _MIN_IDR_CLEARANCE:
                break
            idr_candidate += outward * (spacing * 0.5)
        coords[seg_start:seg_end] = idr_candidate

    coords -= coords.mean(axis=0)
    return coords


# ---------------------------------------------------------------------------
# Multi-chain placement
# ---------------------------------------------------------------------------

def place_chains_grid(
    chain_coords: List[np.ndarray],
    box: List[float],
) -> np.ndarray:
    """
    Place chains on a 3-D grid inside the given periodic box.

    Grid spacing is chosen so all N chains fit in the box volume,
    with a small uniform offset per site to avoid exact overlap.

    Returns:
        (total_beads, 3) array of absolute positions in nm.
    """
    n_chains = len(chain_coords)
    grid_pts = _build_xyzgrid(n_chains, box)
    return _assemble_chains(chain_coords, grid_pts, box)


def place_chains_slab(
    chain_coords: List[np.ndarray],
    box: List[float],
    slab_width: Optional[float] = None,
    clash_cutoff: float = 0.7,
    max_tries: int = 10_000,
) -> np.ndarray:
    """
    Randomly place chains within a slab centred at z = Lz/2, with PBC-aware
    clash detection.

    NOTE: this function is NOT used by the default slab pipeline.  The default
    path in build_all_chains uses the CALVADOS-style deterministic xyzgrid
    placement (_build_xyzgrid + _assemble_chains), which is orders of magnitude
    faster for large MDP chains.  This function is retained as a utility for
    cases where a fully randomised (no-grid) initial configuration is needed.

    Args:
        chain_coords:  List of (N_i, 3) centered chain coordinate arrays.
        box:           Simulation box dimensions [Lx, Ly, Lz] in nm.
        slab_width:    z-extent of the placement region (nm).
                       Defaults to 0.6 * Lz.
        clash_cutoff:  Minimum allowed inter-chain atom distance (nm).
        max_tries:     Maximum placement attempts per chain before giving up.
    """
    try:
        from MDAnalysis.analysis import distances as mda_dist
        _has_mda = True
    except ImportError:
        _has_mda = False

    if slab_width is None:
        slab_width = box[2] * 0.6

    z_lo = box[2] / 2.0 - slab_width / 2.0
    z_hi = box[2] / 2.0 + slab_width / 2.0

    # MDAnalysis box format: [Lx, Ly, Lz, alpha, beta, gamma] (degrees)
    boxfull = np.array([box[0], box[1], box[2], 90.0, 90.0, 90.0], dtype=np.float32)

    placed: List[np.ndarray] = []
    rng = np.random.default_rng()

    for cidx, chain in enumerate(chain_coords):
        for ntry in range(max_tries):
            # Random translation within slab z-region (CALVADOS: draw_starting_vec)
            trans = np.array([
                rng.uniform(0.0, box[0]),
                rng.uniform(0.0, box[1]),
                rng.uniform(z_lo, z_hi),
            ])
            candidate = chain + trans

            # z-boundary check: all atoms must stay within [0, Lz]
            # (x, y have PBC so no wall check needed there)
            if candidate[:, 2].min() < 0.0 or candidate[:, 2].max() > box[2]:
                continue

            # PBC-aware clash check against already-placed atoms (CALVADOS: check_clash)
            if placed:
                xothers = np.vstack(placed).astype(np.float32)
                cand_f  = candidate.astype(np.float32)
                if _has_mda:
                    d = mda_dist.distance_array(cand_f, xothers, boxfull)
                else:
                    # Fallback without MDAnalysis: naive PBC minimum-image distances
                    diff = cand_f[:, None, :] - xothers[None, :, :]
                    diff -= np.round(diff / boxfull[:3]) * boxfull[:3]
                    d    = np.sqrt((diff ** 2).sum(axis=-1))
                if d.min() < clash_cutoff:
                    continue

            placed.append(candidate)
            break
        else:
            raise ValueError(
                f"Failed to place chain {cidx} ({chain.shape[0]} atoms) "
                f"after {max_tries} attempts. "
                f"Try reducing nmol or increasing box / slab_width."
            )

    return np.vstack(placed)


def place_chains_random(
    chain_coords: List[np.ndarray],
    box: List[float],
    clash_cutoff: float = 0.7,
    max_tries: int = 10_000,
) -> np.ndarray:
    """
    Randomly place chains, rejecting positions that clash with existing beads.

    Used for droplet topology.  Mirrors CALVADOS build.random_placement():
    - check_walls: all atoms must lie within [0, box] (hard walls, no PBC wrapping)
    - clash check: PBC-aware minimum-image distances via MDAnalysis distance_array,
      matching CALVADOS build.check_clash behaviour.
    """
    try:
        from MDAnalysis.analysis import distances as mda_dist
        _has_mda = True
    except ImportError:
        _has_mda = False

    box_arr = np.array(box, dtype=float)
    # MDAnalysis box format: [Lx, Ly, Lz, alpha, beta, gamma]
    boxfull = np.array([box[0], box[1], box[2], 90.0, 90.0, 90.0], dtype=np.float32)

    placed: List[np.ndarray] = []
    rng = np.random.default_rng()

    for cidx, chain in enumerate(chain_coords):
        for ntry in range(max_tries):
            # Random position within full box (CALVADOS: draw_starting_vec)
            trans = rng.uniform(0.0, 1.0, size=3) * box_arr
            candidate = chain + trans

            # Hard-wall check: all atoms must be inside [0, box] (CALVADOS: check_walls)
            if candidate.min() < 0.0 or (box_arr - candidate).min() < 0.0:
                continue

            # PBC-aware clash check against already-placed atoms (CALVADOS: check_clash)
            if placed:
                xothers = np.vstack(placed).astype(np.float32)
                cand_f  = candidate.astype(np.float32)
                if _has_mda:
                    d = mda_dist.distance_array(cand_f, xothers, boxfull)
                else:
                    diff = cand_f[:, None, :] - xothers[None, :, :]
                    diff -= np.round(diff / boxfull[:3]) * boxfull[:3]
                    d    = np.sqrt((diff ** 2).sum(axis=-1))
                if d.min() < clash_cutoff:
                    continue

            placed.append(candidate)
            break
        else:
            raise ValueError(
                f"Failed to place chain {cidx} ({chain.shape[0]} atoms) "
                f"after {max_tries} attempts. "
                f"Try reducing nmol or increasing box."
            )

    return np.vstack(placed)


# ---------------------------------------------------------------------------
# High-level builder
# ---------------------------------------------------------------------------

def build_all_chains(config: CGConfig) -> Tuple[np.ndarray, List[dict]]:
    """
    Build all chains from a CGConfig.

    Returns:
        positions: (total_atoms, 3) float64 array in nm.
        chain_meta: list of dicts with 'name', 'start', 'end', 'sequence',
                    'folded_domains', 'comp_type'.
    """
    all_chains: List[np.ndarray] = []
    chain_meta: List[dict] = []
    offset = 0

    # Determine if we should use COM mapping (CALVADOS3 only)
    use_com = config.resolved_force_field == "calvados3"
    
    for comp in config.components:
        seq = comp.get_sequence()
        for _mol_idx in range(comp.nmol):
            if comp.comp_type == ComponentType.MDP and comp.pdb_path:
                if comp.folded_domains:
                    # Use PDB for the folded domain, compact builds for IDR tails.
                    # Raw PDB IDR conformations from MD can span > 10 nm and cause
                    # inter-chain clashes at the ~4 nm initial slab grid spacing.
                    coords = build_mdp_chain_compact_idr(
                        comp.pdb_path, comp.folded_domains, len(seq), use_com=use_com
                    )
                else:
                    coords, _ = build_mdp_chain(comp.pdb_path, use_com=use_com)
            else:
                coords = build_idp_chain(len(seq), method="compact")
            all_chains.append(coords)
            chain_meta.append({
                "name"          : comp.name,
                "start"         : offset,
                "end"           : offset + len(seq),
                "sequence"      : seq,
                "folded_domains": comp.folded_domains,
                "comp_type"     : comp.comp_type,
                "pdb_path"      : comp.pdb_path if comp.comp_type == ComponentType.MDP else None,
            })
            offset += len(seq)

    # Assemble positions according to topology
    if config.topology == TopologyType.SLAB:
        # CALVADOS-style deterministic placement: build a 3-D grid over the slab
        # volume and translate each chain to its grid point.  This mirrors CALVADOS
        # sim.py (build_xyzgrid + grid_counter) and avoids the O(N²) random-
        # placement + clash-detection loop that becomes prohibitively slow for large
        # MDP chains.  Any residual inter-chain overlaps are resolved by the
        # subsequent energy minimisation step.
        slab_w = config.slab_width if config.slab_width is not None else config.box[2] * 0.6
        grid_pts = _build_xyzgrid(len(all_chains), [config.box[0], config.box[1], slab_w])
        grid_pts += np.array([0.0, 0.0, config.box[2] / 2.0 - slab_w / 2.0])
        positions = _assemble_chains(all_chains, grid_pts, config.box)
    elif config.topology == TopologyType.DROPLET:
        r = config.droplet_radius or (min(config.box) * 0.4)
        droplet_box = [r * 2.0, r * 2.0, r * 2.0]
        positions = place_chains_random(all_chains, droplet_box)
        # place_chains_random returns coordinates in the local [0, 2r] droplet
        # box. Shift that local box so its centre coincides with the real box.
        centre = np.array(config.box, dtype=float) / 2.0
        local_centre = np.array(droplet_box, dtype=float) / 2.0
        positions += centre - local_centre
    else:
        positions = place_chains_grid(all_chains, config.box)

    return positions.astype(np.float64), chain_meta


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _get_idr_segments(
    n_res: int,
    folded_domains: List[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """
    Return IDR segment ranges as (start, end) tuples (0-indexed, end exclusive).

    Example: n_res=193, folded_domains=[(22, 96)]
      → [(0, 21), (96, 193)]
    """
    sorted_doms = sorted(folded_domains)
    segments: List[Tuple[int, int]] = []
    prev_end = 0
    for dom_s, dom_e in sorted_doms:
        if prev_end < dom_s - 1:
            segments.append((prev_end, dom_s - 1))
        prev_end = dom_e
    if prev_end < n_res:
        segments.append((prev_end, n_res))
    return segments


def _build_spiral(n: int, d: float = 0.38) -> np.ndarray:
    """Archimedean spiral in xy, linear in z."""
    coords = np.zeros((n, 3))
    for i in range(n):
        theta = np.sqrt(i / max(n, 1)) * 2.0 * np.pi
        r     = d * np.sqrt(i)
        coords[i, 0] = r * np.cos(theta)
        coords[i, 1] = r * np.sin(theta)
    coords[:, 2] = np.linspace(-n * d / 2.0, n * d / 2.0, n)
    coords -= coords.mean(axis=0)
    return coords


def _build_linear(n: int, d: float = 0.38) -> np.ndarray:
    """Straight chain along z axis."""
    coords = np.zeros((n, 3))
    coords[:, 2] = d * np.arange(n)
    coords -= coords.mean(axis=0)
    return coords


def _build_compact(n: int, d: float = 0.38) -> np.ndarray:
    """
    Simple cubic lattice filling.
    Matches CALVADOS build.build_compact implementation.
    """
    N = int(np.ceil(np.cbrt(n)) - 1)
    xs = []
    i, j, k = 0, 0, 0
    di, dj, dk = 1, 1, 1  # direction
    cti, ctj, ctk = 0, 0, 0

    for idx in range(n):
        xs.append([i, j, k])
        if ctk == N:
            if ctj == N:
                i += di
                cti += 1
                ctj = 0
                dj *= -1
            else:
                j += dj
                ctj += 1
            ctk = 0
            dk *= -1
        else:
            k += dk
            ctk += 1
    xs = (np.array(xs) - 0.5 * N) * d
    return xs


def _build_xyzgrid(n: int, box: List[float]) -> np.ndarray:
    """
    3-D staggered grid of N points, scaled proportionally to the box dimensions.
    Mirrors the algorithm in CALVADOS build.build_xyzgrid with staggered offsets.
    
    The staggered pattern:
    - Adjacent z-planes are offset by (dx/2, dy/2) in xy
    - Adjacent points within xy-plane have alternating z-offset (dz/2)
    """
    n = int(np.ceil(n))
    if n < 1:
        n = 1

    box = np.array(box, dtype=float)
    # Calculate grid dimensions proportionally to box dimensions
    r = box / np.sum(box)
    a = np.cbrt(n / np.prod(r))
    n_float = a * r
    nxyz = np.maximum(np.floor(n_float), 1).astype(int)
    
    # Adjust grid dimensions to fit all N points
    while np.prod(nxyz) < n:
        ndeviation = n_float / nxyz
        devmax = np.argmax(ndeviation)
        nxyz[devmax] += 1
    while np.prod(nxyz) > n:
        nmax = np.argmax(nxyz)
        nxyz[nmax] -= 1
        if np.prod(nxyz) < n:
            nxyz[nmax] += 1
            break

    # Grid spacing
    dx = box[0] / nxyz[0]
    dy = box[1] / nxyz[1]
    dz = box[2] / nxyz[2]

    pts = []
    x, y, z = 0.0, 0.0, 0.0
    ctx, cty, ctz = 0, 0, 0
    zplane = 1  # Controls xy-offset for alternating z-planes
    xyplane = 1  # Controls z-offset within xy-plane

    for _ in np.arange(n):
        # Staggered offset: alternate z-planes offset by (dx/2, dy/2)
        if zplane > 0:
            xshift = 0
            yshift = 0
        else:
            xshift = dx / 2
            yshift = dy / 2

        # Staggered offset: alternating z-offset within xy-plane
        if xyplane < 0:
            zshift = dz / 2
        else:
            zshift = 0

        pts.append([x + xshift, y + yshift, z + zshift])

        # Move to next grid point in x
        ctx += 1
        x += dx

        # Update xyplane based on position parity (alternating pattern)
        if (ctx % 2 == cty % 2):
            xyplane = 1
        else:
            xyplane = -1

        # Move to next row in y
        if ctx == nxyz[0]:
            ctx = 0
            x = 0.0

            cty += 1
            y += dy

            # Move to next layer in z
            if cty == nxyz[1]:
                ctx = 0
                cty = 0
                x = 0.0
                y = 0.0

                ctz += 1
                z += dz

                # Flip zplane for staggered offset in next xy-layer
                zplane = -zplane

    return np.array(pts)


def _assemble_chains(
    chain_coords: List[np.ndarray],
    grid_pts: np.ndarray,
    box: List[float],
) -> np.ndarray:
    """Translate each chain to its grid position and concatenate."""
    all_pos = []
    for i, chain in enumerate(chain_coords):
        pos = chain + grid_pts[i]
        all_pos.append(pos)
    return np.vstack(all_pos)
