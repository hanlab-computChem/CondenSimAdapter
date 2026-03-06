"""
Mpipi-Recharged force field.

Physics:
  - Harmonic CA-CA backbone bonds (IDR regions)
  - ENM for MDP folded domains (globular regions)
  - Wang-Frenkel short-range non-bonded (tabulated 21x21 parameters)
  - Yukawa electrostatics (tabulated 21x21 A-matrix)
  - Globular-domain scaling factor (glob_factor)

Reference: Joseph et al. eLife 2021; recharged version -- Bremer et al. Nat. Chem. 2022.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit

from .base import CGForceField

_DATA = Path(__file__).parent / "data"

# Mpipi residue ordering: 20 amino acids + rU (RNA uridine, kept for table
# completeness but excluded from protein-only simulations via a fixed index)
_MPIPI_ORDER = [
    "pM", "pG", "pK", "pT", "pR", "pA", "pD", "pE",
    "pY", "pV", "pL", "pQ", "pW", "pF", "pS", "pH",
    "pN", "pP", "pC", "pI",
    # index 20: rU (RNA) -- unused in protein-only mode
]
_N_MPIPI = 21   # table size (includes rU)

# Map standard three-letter code -> Mpipi index
_AA3_TO_MPIPI: Dict[str, int] = {
    "MET": 0,  "GLY": 1,  "LYS": 2,  "THR": 3,  "ARG": 4,
    "ALA": 5,  "ASP": 6,  "GLU": 7,  "TYR": 8,  "VAL": 9,
    "LEU": 10, "GLN": 11, "TRP": 12, "PHE": 13, "SER": 14,
    "HIS": 15, "ASN": 16, "PRO": 17, "CYS": 18, "ILE": 19,
}

# Mpipi residue masses
_MPIPI_MASS: Dict[str, float] = {
    "MET": 131.20, "GLY": 57.05,  "LYS": 128.20, "THR": 101.10, "ARG": 156.20,
    "ALA": 71.08,  "ASP": 115.10, "GLU": 129.10, "TYR": 163.20, "VAL": 99.07,
    "LEU": 113.20, "GLN": 128.10, "TRP": 186.20, "PHE": 147.20, "SER": 87.08,
    "HIS": 137.10, "ASN": 114.10, "PRO": 97.12,  "CYS": 103.10, "ILE": 113.20,
}

_K_BOND = 8031.0   # kJ/mol/nm^2  (harmonic backbone)
_D_IDR  = 0.381    # nm  (IDR bond length)

# ENM parameters for folded domains
_ENM_K       = 8031.0   # kJ/mol/nm^2
_ENM_CUTOFF  = 0.75     # nm
_ENM_MIN_SEP = 1        # MPIPI uses no sequence separation (neighbors included)


class MpipiFF(CGForceField):
    """Mpipi-Recharged force field for protein condensates."""

    def __init__(self):
        data = np.loadtxt(_DATA / "recharged_params.txt")
        # First 21*21*3 values: WF table (eps, sigma, mu) stored row-major [i,j,k]
        # Last  21*21   values: Yukawa A-matrix
        n3 = _N_MPIPI * _N_MPIPI * 3
        self._wf_params   = data[:n3].tolist()
        self._yukawa_A    = data[n3:].tolist()

    # ------------------------------------------------------------------

    def add_masses(self, system: mm.System, chain_meta: List[dict]) -> None:
        for meta in chain_meta:
            from ..topology import ONE_TO_THREE
            for aa1 in meta["sequence"]:
                three = ONE_TO_THREE.get(aa1, "GLY")
                m = _MPIPI_MASS.get(three, 57.05)
                system.addParticle(m * unit.amu)

    def create_nonbonded_forces(
        self,
        topology: app.Topology,
        chain_meta: List[dict],
        temperature: float,
        ionic: float,
        debye_length: float = 1.0,
    ) -> List[mm.Force]:
        kappa = 1.0 / debye_length   # nm^-1; can be overridden by caller

        # --- Wang-Frenkel short-range ---
        wf_str = (
            "glob_factor * step(rc-r) * epsilon * alpha * ((sigma/r)^(2*mu)-1) * ((rc/r)^(2*mu)-1)^2;"
            "alpha = 2*(3^(2*mu))*((3)/(2*((3^(2*mu))-1)))^3;"
            "rc = 3*sigma;"
            "glob_factor = select(globular1*globular2, 0.7,"
            "              select(globular1+globular2, sqrt(0.7), 1.0));"
            "epsilon = wf_table(index1, index2, 0);"
            "sigma   = wf_table(index1, index2, 1);"
            "mu      = floor(wf_table(index1, index2, 2))"
        )
        wf = mm.CustomNonbondedForce(wf_str)
        wf.addTabulatedFunction(
            "wf_table",
            mm.Discrete3DFunction(_N_MPIPI, _N_MPIPI, 3, self._wf_params)
        )
        wf.addPerParticleParameter("index")
        wf.addPerParticleParameter("globular")
        wf.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        wf.setCutoffDistance(2.5 * unit.nanometer)
        wf.setForceGroup(0)

        # --- Yukawa electrostatics ---
        yu_str = "(A_table(index1, index2)/r) * exp(-kappa*r)"
        yu = mm.CustomNonbondedForce(yu_str)
        yu.addTabulatedFunction(
            "A_table",
            mm.Discrete2DFunction(_N_MPIPI, _N_MPIPI, self._yukawa_A)
        )
        yu.addPerParticleParameter("index")
        yu.addGlobalParameter("kappa", kappa / unit.nanometer)
        yu.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        yu.setCutoffDistance(3.5 * unit.nanometer)
        yu.setForceGroup(1)

        # Identify folded (globular) atoms
        folded_atoms = self._get_folded_set(chain_meta)

        for atom in topology.atoms():
            from ..topology import THREE_TO_ONE
            aa1 = THREE_TO_ONE.get(atom.residue.name, "G")
            from ..topology import ONE_TO_THREE
            three = ONE_TO_THREE.get(aa1, "GLY")
            idx = _AA3_TO_MPIPI.get(three, 1)   # fallback to GLY
            is_glob = 1 if atom.index in folded_atoms else 0
            wf.addParticle([idx, is_glob])
            yu.addParticle([idx])

        bonds = [(b[0].index, b[1].index) for b in topology.bonds()]
        wf.createExclusionsFromBonds(bonds, 1)
        yu.createExclusionsFromBonds(bonds, 1)

        return [wf, yu]

    # ------------------------------------------------------------------

    def build_harmonic_bonds(
        self,
        topology: app.Topology,
        r0: float = _D_IDR,
        k: float = _K_BOND,
    ) -> mm.HarmonicBondForce:
        """MPIPI CA-CA harmonic bonds using MPIPI-specific defaults (k=8031.0, r0=0.381)."""
        return super().build_harmonic_bonds(topology, r0=r0, k=k)

    # ------------------------------------------------------------------

    def build_enm_bonds(
        self,
        positions: np.ndarray,
        chain_meta: List[dict],
        restraint_type: str = "harmonic",
        k: float = None,
        cutoff: float = None,
    ) -> mm.Force:
        """
        Elastic Network Model (ENM) for folded domains.
        
        MPIPI-specific implementation using original MPIPI parameters:
        - k = 8031 kJ/mol/nm^2
        - cutoff = 0.75 nm
        - min_seq_sep = 1 (no sequence separation, matching original behavior)
        
        Reference: Original MPIPI implementation uses KDTree without sequence separation.
        """
        from typing import Optional
        
        # Use MPIPI-specific defaults if not provided
        if k is None:
            k = _ENM_K
        if cutoff is None:
            cutoff = _ENM_CUTOFF
        min_seq_sep = _ENM_MIN_SEP
        
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

    # ------------------------------------------------------------------

    @staticmethod
    def _get_folded_set(chain_meta: List[dict]) -> set:
        """Return the set of absolute atom indices inside folded domains."""
        folded = set()
        for meta in chain_meta:
            cs = meta["start"]
            for (dom_s, dom_e) in meta["folded_domains"]:
                for i in range(cs + dom_s - 1, cs + dom_e):
                    folded.add(i)
        return folded
