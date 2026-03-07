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
    Place chains on a staggered 3D grid within a slab centred at z = Lz/2.

    - Uses staggered grid placement (alternate layers offset by half grid spacing)
    - Slab width defaults to 0.75 * Lz if not provided (assuming z is the long axis)
    - Grid is centered in z at box[2] / 2
    """
    if slab_width is None:
        slab_width = box[2] * 0.6  # Default: 60% of z box dimension
    slab_box = [box[0], box[1], slab_width]
    grid_pts = _build_xyzgrid(len(chain_coords), slab_box)
    # shift to centre of z axis
    grid_pts += np.array([0.0, 0.0, box[2] / 2.0 - slab_width / 2.0])
    return _assemble_chains(chain_coords, grid_pts, box)


def place_chains_random(
    chain_coords: List[np.ndarray],
    box: List[float],
    clash_cutoff: float = 0.7,
    max_tries: int = 10_000,
) -> np.ndarray:
    """
    Randomly place chains, rejecting positions that clash with existing beads.

    Used for droplet topology (no periodic boundary during placement).
    Mirrors the algorithm in CALVADOS build.random_placement().
    """
    placed_positions: List[np.ndarray] = []
    all_coords: List[np.ndarray] = []

    for chain in chain_coords:
        for ntry in range(max_tries):
            # random translation within full box (like draw_starting_vec)
            trans = np.random.uniform(0, box, size=3)
            candidate = chain + trans

            # check if outside box (like check_walls)
            if np.min(candidate) < 0:
                continue
            if np.min(np.array(box) - candidate) < 0:
                continue

            # clash check against already-placed beads (like check_clash)
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
            raise ValueError(
                f"Tried {max_tries}x to add molecule. Giving up."
            )

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

    # Determine if we should use COM mapping (CALVADOS3 only)
    use_com = config.resolved_force_field == "calvados3"
    
    for comp in config.components:
        seq = comp.get_sequence()
        for _mol_idx in range(comp.nmol):
            if comp.comp_type == ComponentType.MDP and comp.pdb_path:
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
        positions = place_chains_slab(all_chains, config.box, slab_width=config.slab_width)
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
