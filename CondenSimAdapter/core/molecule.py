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

    Args:
        pdb_path: Path to PDB file.
        use_com:  If True, use residue centre-of-mass instead of CA.

    Returns:
        coords: (N, 3) float64 array in nm, mass-centred at origin.
        sequence: one-letter amino acid string.
    """
    import MDAnalysis as mda
    from MDAnalysis.lib.util import convert_aa_code

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
) -> np.ndarray:
    """
    Place chains on a grid within a slab centred at z = Lz/2.

    The slab width defaults to min(Lx, Ly) / 2 if not provided.
    """
    if slab_width is None:
        slab_width = min(box[0], box[1]) / 2.0
    slab_box = [box[0], box[1], slab_width]
    grid_pts = _build_xyzgrid(len(chain_coords), slab_box)
    # shift to centre of z axis
    grid_pts += np.array([0.0, 0.0, box[2] / 2.0 - slab_width / 2.0])
    return _assemble_chains(chain_coords, grid_pts, box)


def place_chains_random(
    chain_coords: List[np.ndarray],
    box: List[float],
    clash_cutoff: float = 0.8,
    max_tries: int = 10_000,
) -> np.ndarray:
    """
    Randomly place chains, rejecting positions that clash with existing beads.

    Used for droplet topology (no periodic boundary during placement).
    """
    placed_positions: List[np.ndarray] = []
    all_coords: List[np.ndarray] = []

    for chain in chain_coords:
        chain_size = np.ptp(chain, axis=0).max() if len(chain) > 1 else 0.0

        for _ in range(max_tries):
            # random translation within box, keeping chain away from walls
            margin = chain_size / 2.0 + 0.5
            lo = margin
            hi = np.array(box) - margin
            if np.any(hi <= lo):
                hi = np.array(box) * 0.9
                lo = np.array(box) * 0.1
            trans = np.random.uniform(lo, hi)
            candidate = chain + trans

            # clash check against already-placed beads
            if all_coords:
                all_placed = np.vstack(all_coords)
                diffs = candidate[:, None, :] - all_placed[None, :, :]
                dists = np.sqrt((diffs ** 2).sum(axis=-1))
                if dists.min() < clash_cutoff:
                    continue
            all_coords.append(candidate)
            placed_positions.append(candidate)
            break
        else:
            # fallback: place without clash check
            trans = np.random.uniform(0, box, size=3)
            candidate = chain + trans
            all_coords.append(candidate)
            placed_positions.append(candidate)

    return np.vstack(placed_positions)


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

    for comp in config.components:
        seq = comp.get_sequence()
        for _mol_idx in range(comp.nmol):
            if comp.comp_type == ComponentType.MDP and comp.pdb_path:
                coords, _ = build_mdp_chain(comp.pdb_path)
            else:
                coords = build_idp_chain(len(seq), method="spiral")
            all_chains.append(coords)
            chain_meta.append({
                "name"          : comp.name,
                "start"         : offset,
                "end"           : offset + len(seq),
                "sequence"      : seq,
                "folded_domains": comp.folded_domains,
                "comp_type"     : comp.comp_type,
            })
            offset += len(seq)

    # Assemble positions according to topology
    if config.topology == TopologyType.SLAB:
        positions = place_chains_slab(all_chains, config.box)
    elif config.topology == TopologyType.DROPLET:
        r = config.droplet_radius or (min(config.box) * 0.4)
        droplet_box = [r * 2.0, r * 2.0, r * 2.0]
        positions = place_chains_random(all_chains, droplet_box)
        # centre in actual box
        centre = np.array(config.box) / 2.0
        positions += centre
    else:
        positions = place_chains_grid(all_chains, config.box)

    return positions.astype(np.float64), chain_meta


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

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
    """Simple cubic lattice filling."""
    side = int(np.ceil(n ** (1.0 / 3.0)))
    coords = []
    for ix in range(side):
        for iy in range(side):
            for iz in range(side):
                coords.append([ix * d, iy * d, iz * d])
                if len(coords) >= n:
                    break
            if len(coords) >= n:
                break
        if len(coords) >= n:
            break
    coords = np.array(coords[:n], dtype=np.float64)
    coords -= coords.mean(axis=0)
    return coords


def _build_xyzgrid(n: int, box: List[float]) -> np.ndarray:
    """
    3-D grid of N points, scaled proportionally to the box dimensions.
    Mirrors the algorithm in CALVADOS build.build_xyzgrid.
    """
    n = int(np.ceil(n))
    box = np.array(box, dtype=float)
    r = box / box.sum()
    a = (n / r.prod()) ** (1.0 / 3.0)
    nxyz = np.floor(a * r).astype(int)
    # at least 1 in each dimension
    nxyz = np.maximum(nxyz, 1)
    # generate enough points
    pts = []
    xshift = box[0] / (2.0 * nxyz[0])
    yshift = box[1] / (2.0 * nxyz[1])
    zshift = box[2] / (2.0 * nxyz[2])
    for xi in range(nxyz[0]):
        for yi in range(nxyz[1]):
            for zi in range(nxyz[2]):
                x = xi * box[0] / nxyz[0] + xshift
                y = yi * box[1] / nxyz[1] + yshift
                z = zi * box[2] / nxyz[2] + zshift
                pts.append([x, y, z])
    pts = np.array(pts)
    if len(pts) < n:
        # pad with random
        extra = np.random.rand(n - len(pts), 3) * box
        pts = np.vstack([pts, extra])
    return pts[:n]


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
