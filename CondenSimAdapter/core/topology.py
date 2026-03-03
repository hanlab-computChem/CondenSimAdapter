"""
Build an OpenMM Topology for CG protein systems (one CA bead per residue).
"""

from __future__ import annotations

import numpy as np
from typing import Dict, List, Optional, Tuple

import openmm.app as app
import openmm.unit as unit


# Three-letter -> one-letter conversion table (proteins only)
THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
ONE_TO_THREE: Dict[str, str] = {v: k for k, v in THREE_TO_ONE.items()}

# Residue masses (g/mol)
RESIDUE_MASS: Dict[str, float] = {
    "A": 71.08,  "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14,
    "E": 129.11, "Q": 128.13, "G": 57.05,  "H": 137.14, "I": 113.16,
    "L": 113.16, "K": 128.17, "M": 131.20, "F": 147.18, "P": 97.12,
    "S": 87.08,  "T": 101.11, "W": 186.22, "Y": 163.18, "V": 99.13,
}


def build_topology(
    chain_meta: List[dict],
    positions: np.ndarray,
    box: List[float],
) -> Tuple[app.Topology, np.ndarray]:
    """
    Construct an OpenMM Topology for a CG system (one CA per residue).

    Args:
        chain_meta: list of dicts from molecule.build_all_chains(), each with
                    keys 'name', 'start', 'end', 'sequence', 'folded_domains'.
        positions:  (total_atoms, 3) float64 in nm.
        box:        [Lx, Ly, Lz] in nm.

    Returns:
        topology:   OpenMM Topology object.
        positions:  (N, 3) Quantity in nm (same data, wrapped in unit).
    """
    top = app.Topology()

    # Use the virtual 'element' for CG beads; we pick a lightweight element
    # and override the mass later at the System level.
    ca_elem = app.element.carbon   # placeholder; mass overridden in force field

    for meta in chain_meta:
        chain_obj = top.addChain()
        seq = meta["sequence"]
        for res_idx, aa in enumerate(seq):
            three = ONE_TO_THREE.get(aa, "GLY")
            res = top.addResidue(three, chain_obj)
            top.addAtom("CA", ca_elem, res)

    # Add backbone bonds (CA i -- CA i+1 within each chain)
    atoms = list(top.atoms())
    for meta in chain_meta:
        for i in range(meta["start"], meta["end"] - 1):
            top.addBond(atoms[i], atoms[i + 1])

    # Periodic box
    top.setPeriodicBoxVectors([
        unit.Quantity([box[0], 0, 0], unit.nanometer),
        unit.Quantity([0, box[1], 0], unit.nanometer),
        unit.Quantity([0, 0, box[2]], unit.nanometer),
    ])

    pos_quantity = positions * unit.nanometer
    return top, pos_quantity


def get_masses(chain_meta: List[dict]) -> np.ndarray:
    """Return per-bead masses (amu) in the same order as the topology."""
    masses = []
    for meta in chain_meta:
        for aa in meta["sequence"]:
            masses.append(RESIDUE_MASS.get(aa, 57.05))
    return np.array(masses, dtype=np.float64)


def get_folded_atom_ranges(chain_meta: List[dict]) -> List[Tuple[int, int]]:
    """
    Return absolute atom-index ranges for all folded domains across all chains.

    Returns list of (start_atom_idx, end_atom_idx) (end exclusive).
    """
    ranges = []
    for meta in chain_meta:
        chain_start = meta["start"]
        for (dom_start, dom_end) in meta["folded_domains"]:
            # convert 1-based inclusive domain indices to 0-based absolute
            abs_start = chain_start + dom_start - 1
            abs_end   = chain_start + dom_end       # exclusive
            ranges.append((abs_start, abs_end))
    return ranges
