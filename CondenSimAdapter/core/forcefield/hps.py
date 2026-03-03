"""
HPS-Urry force field (Urry hydrophobicity scale).

Physics:
  - Harmonic CA-CA backbone bonds
  - Ashbaugh-Hatch short-range non-bonded with pair-specific sigma/lambda
    (implemented via Discrete2DFunction for the 20x20 parameter tables)
  - Debye-Hückel electrostatics
  - ENM for MDP folded domains

Reference: Regy et al. Protein Sci. 2021; Dignon et al. PLOS Comput. Biol. 2018.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit

from .base import CGForceField, debye_huckel_params

_DATA = Path(__file__).parent / "data"

# Canonical amino acid ordering for HPS index (0-19)
_AA_ORDER = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
]
_AA_INDEX: Dict[str, int] = {aa: i for i, aa in enumerate(_AA_ORDER)}
_N_AA = 20


class HPSFF(CGForceField):
    """
    HPS-Urry force field.

    The Ashbaugh-Hatch interaction uses pairwise sigma and lambda values
    tabulated for all 20x20 amino acid combinations.
    """

    EPS_AH   = 0.8368    # kJ/mol
    RC_AH    = 2.0       # nm
    RC_DH    = 3.5       # nm
    DH_DIELECTRIC = 80.0
    # HPS-Urry correction: lambda_eff = lambda - delta
    DELTA    = 0.08
    K_BOND   = 8368.0    # kJ/mol/nm^2

    def __init__(self):
        self._sigma_table: np.ndarray = np.zeros((_N_AA, _N_AA))
        self._lambda_table: np.ndarray = np.zeros((_N_AA, _N_AA))
        self._charges: Dict[str, float] = {}
        self._load_params()

    def _load_params(self) -> None:
        """Build symmetric 20x20 sigma and lambda tables from the CSV."""
        fname = _DATA / "HPS_Urry_parameters.csv"
        with open(fname, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                aa1, aa2 = row["atom_type1"], row["atom_type2"]
                sigma  = float(row["sigma"])
                lam    = float(row["lambda"]) - self.DELTA
                i = _AA_INDEX.get(aa1, 0)
                j = _AA_INDEX.get(aa2, 0)
                self._sigma_table[i, j] = sigma
                self._sigma_table[j, i] = sigma
                self._lambda_table[i, j] = lam
                self._lambda_table[j, i] = lam

        # Residue charges (standard protonation at pH 7)
        self._charges = {
            "ARG":  1.0, "LYS":  1.0,
            "ASP": -1.0, "GLU": -1.0,
            "HIS":  0.0,
        }

    def _get_type_index(self, atom: app.topology.Atom) -> int:
        return _AA_INDEX.get(atom.residue.name, 7)   # 7 = GLY fallback

    def _get_charge(self, atom: app.topology.Atom) -> float:
        return self._charges.get(atom.residue.name, 0.0)

    # ------------------------------------------------------------------

    def create_nonbonded_forces(
        self,
        topology: app.Topology,
        chain_meta: List[dict],
        temperature: float,
        ionic: float,
    ) -> List[mm.Force]:
        forces = []

        # Flatten tables into row-major lists for Discrete2DFunction
        sigma_flat  = self._sigma_table.flatten().tolist()
        lambda_flat = self._lambda_table.flatten().tolist()

        # --- Ashbaugh-Hatch via tabulated sigma / lambda ---
        eps = self.EPS_AH
        rc  = self.RC_AH
        expr = (
            f"{eps}*select(step(r-2^(1/6)*sigma_ah),"
            f"4*lam_ah*((sigma_ah/r)^12-(sigma_ah/r)^6-shift),"
            f"4*((sigma_ah/r)^12-(sigma_ah/r)^6-lam_ah*shift)+(1-lam_ah));"
            f"sigma_ah=sigma_table(type1,type2);"
            f"lam_ah=lambda_table(type1,type2);"
            f"shift=(sigma_ah/{rc})^12-(sigma_ah/{rc})^6"
        )
        ah = mm.CustomNonbondedForce(expr)
        ah.addTabulatedFunction(
            "sigma_table",
            mm.Discrete2DFunction(_N_AA, _N_AA, sigma_flat)
        )
        ah.addTabulatedFunction(
            "lambda_table",
            mm.Discrete2DFunction(_N_AA, _N_AA, lambda_flat)
        )
        ah.addPerParticleParameter("type1")   # integer index
        # Note: Discrete2DFunction needs two separate per-particle params
        # OpenMM syntax: sigma_table(type1, type2) references params 0 and 1
        # so we use a single int parameter indexed as (atom_type) twice
        # We use a trick: one per-particle int, function f(i,j) with same index
        ah.addPerParticleParameter("type2")
        ah.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        ah.setCutoffDistance(rc * unit.nanometer)
        ah.setForceGroup(0)

        # --- Debye-Hückel electrostatics ---
        eps_yu, k_yu = debye_huckel_params(temperature, ionic)
        shift_dh = float(np.exp(-k_yu * self.RC_DH) / self.RC_DH)
        dh_expr = (
            f"q1*q2*{eps_yu:.6f}*(exp(-{k_yu:.6f}*r)/r-{shift_dh:.6e})*"
            f"step({self.RC_DH}-r)"
        )
        dh = mm.CustomNonbondedForce(dh_expr)
        dh.addPerParticleParameter("q1")
        dh.addPerParticleParameter("q2")
        dh.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        dh.setCutoffDistance(self.RC_DH * unit.nanometer)
        dh.setForceGroup(1)

        # Add per-particle values
        for atom in topology.atoms():
            idx = self._get_type_index(atom)
            ah.addParticle([idx, idx])
            q = self._get_charge(atom)
            dh.addParticle([q, q])

        bonds = [(b[0].index, b[1].index) for b in topology.bonds()]
        ah.createExclusionsFromBonds(bonds, 1)
        dh.createExclusionsFromBonds(bonds, 1)

        forces.extend([ah, dh])
        return forces
