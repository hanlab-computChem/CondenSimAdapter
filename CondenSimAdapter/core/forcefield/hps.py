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
from typing import Dict, List, Optional, Tuple

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
    RC_AH    = 2.0       # nm (fallback, but actual cutoff uses 4*sigma dynamic)
    RC_DH    = 3.5       # nm
    DH_DIELECTRIC = 80.0
    # HPS-Urry correction: lambda_eff = lambda - delta
    DELTA    = 0.08
    K_BOND   = 8368.0    # kJ/mol/nm^2
    
    # ENM parameters for folded domains (MDP support)
    # Matching CALVADOS standard: k=700, cutoff=0.9, min_seq_sep=3
    ENM_K       = 700.0     # kJ/mol/nm^2
    ENM_CUTOFF  = 0.9       # nm
    ENM_MIN_SEP = 3         # minimum sequence separation

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
        # HIS = 0.5 matches original HPS-Urry implementation
        self._charges = {
            "ARG":  1.0, "LYS":  1.0,
            "ASP": -1.0, "GLU": -1.0,
            "HIS":  0.5,
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
        # Dynamic cutoff: 4*sigma (matches original HPS-Urry implementation)
        eps = self.EPS_AH
        lj_at_cutoff = 4*eps*((1/4)**12 - (1/4)**6)
        expr = (
            f"(f1+f2-offset)*step(4*sigma_ah-r);"
            f"offset=lam_ah*{lj_at_cutoff};"
            f"f1=(lj+(1-lam_ah)*{eps})*step(2^(1/6)*sigma_ah-r);"
            f"f2=lam_ah*lj*step(r-2^(1/6)*sigma_ah);"
            f"lj=4*{eps}*((sigma_ah/r)^12-(sigma_ah/r)^6);"
            f"sigma_ah=sigma_table(atom_type1,atom_type2);"
            f"lam_ah=lambda_table(atom_type1,atom_type2)"
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
        ah.addPerParticleParameter("atom_type")  # integer index
        ah.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        # Dynamic cutoff: 4*max(sigma) to cover all pair interactions
        max_sigma = float(np.max(self._sigma_table))
        ah.setCutoffDistance(4.0 * max_sigma * unit.nanometer)
        ah.setForceGroup(0)

        # --- Debye-Hückel electrostatics ---
        # HPS-Urry uses fixed Debye length of 1.0 nm (matching original implementation)
        # instead of computing from ionic strength
        _LDBY_HPS = 1.0  # nm, fixed Debye length for HPS-Urry
        eps_yu, _ = debye_huckel_params(temperature, ionic)
        k_yu = 1.0 / _LDBY_HPS  # inverse Debye length
        shift_dh = float(np.exp(-k_yu * self.RC_DH) / self.RC_DH)
        dh_expr = (
            f"q1*q2*{eps_yu:.6f}*(exp(-{k_yu:.6f}*r)/r-{shift_dh:.6e})*"
            f"step({self.RC_DH}-r)"
        )
        dh = mm.CustomNonbondedForce(dh_expr)
        dh.addPerParticleParameter("q")
        dh.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        dh.setCutoffDistance(self.RC_DH * unit.nanometer)
        dh.setForceGroup(1)

        # Add per-particle values
        for atom in topology.atoms():
            idx = self._get_type_index(atom)
            ah.addParticle([idx])  # single atom type index
            q = self._get_charge(atom)
            dh.addParticle([q])

        bonds = [(b[0].index, b[1].index) for b in topology.bonds()]
        ah.createExclusionsFromBonds(bonds, 1)
        dh.createExclusionsFromBonds(bonds, 1)

        forces.extend([ah, dh])
        return forces

    # ------------------------------------------------------------------
    # ENM for folded domains (override base to use HPS-specific parameters)
    # ------------------------------------------------------------------

    def build_enm_bonds(
        self,
        positions: np.ndarray,
        chain_meta: List[dict],
        restraint_type: str = "harmonic",
        k: float = None,
        cutoff: float = None,
    ) -> Optional[mm.Force]:
        """
        Elastic Network Model (ENM) for folded domains.
        
        HPS-specific implementation using CALVADOS-standard parameters:
        - k = 700 kJ/mol/nm^2 (matching original CALVADOS/HPS)
        - cutoff = 0.9 nm
        - min_seq_sep = 3
        
        Reference: Original CALVADOS and OpenABC HPS implementation
        """
        # Use HPS-specific defaults if not provided
        if k is None:
            k = self.ENM_K
        if cutoff is None:
            cutoff = self.ENM_CUTOFF
        min_seq_sep = self.ENM_MIN_SEP
        
        if restraint_type == "go":
            expr = "k*(5*(s/r)^12-6*(s/r)^10); s=s; k=k"
            cs = mm.CustomBondForce(expr)
            cs.addPerBondParameter("s")
            cs.addPerBondParameter("k")
        else:
            cs = mm.HarmonicBondForce()
        
        cs.setUsesPeriodicBoundaryConditions(True)
        n_bonds = 0
        
        for meta in chain_meta:
            if not meta["folded_domains"]:
                continue
            chain_start = meta["start"]
            for (dom_s, dom_e) in meta["folded_domains"]:
                a0 = chain_start + dom_s - 1   # absolute, 0-based
                a1 = chain_start + dom_e        # exclusive
                indices = list(range(a0, a1))
                for ii in range(len(indices)):
                    for jj in range(ii + min_seq_sep, len(indices)):
                        gi, gj = indices[ii], indices[jj]
                        d = float(np.linalg.norm(positions[gi] - positions[gj]))
                        if d <= cutoff:
                            if restraint_type == "go":
                                cs.addBond(gi, gj,
                                           [d * unit.nanometer,
                                            k * unit.kilojoule_per_mole])
                            else:
                                cs.addBond(
                                    gi, gj,
                                    d * unit.nanometer,
                                    k * unit.kilojoule_per_mole / unit.nanometer ** 2)
                            n_bonds += 1
        
        return cs if n_bonds > 0 else None
