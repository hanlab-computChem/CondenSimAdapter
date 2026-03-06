"""
CALVADOS2 / CALVADOS3 force field.

Physics:
  - Harmonic CA-CA backbone bonds (k = 8368 kJ/mol/nm^2)
  - Ashbaugh-Hatch (AH) short-range non-bonded (with ID-based mixing)
  - Yukawa (Debye-Hückel) electrostatics
  - Go-like ENM for folded domains (MDP only)

"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit

from .base import CGForceField, debye_huckel_params

# Parameter files relative to this directory
_DATA = Path(__file__).parent / "data"


class CalvadosFF(CGForceField):
    """
    Protein-only CALVADOS force field (version 2 or 3).

    Parameters are read from the residues CSV bundled in forcefield/data/.
    """

    # CALVADOS default LJ parameters
    EPS_LJ   = 0.8368    # kJ/mol  (0.2 kcal/mol)
    RC_LJ    = 2.0       # nm
    RC_YU    = 4.0       # nm
    K_BOND   = 8368.0    # kJ/mol/nm^2  (CALVADOS standard; matches original implementation)

    def __init__(self, version: int = 2):
        if version not in (2, 3):
            raise ValueError("CALVADOS version must be 2 or 3.")
        self.version = version
        self._params: Dict[str, dict] = {}
        self._load_params()

    def _load_params(self) -> None:
        fname = _DATA / f"residues_CALVADOS{self.version}.csv"
        with open(fname, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                one = row["one"]
                self._params[one] = {
                    "sigma" : float(row["sigmas"]),
                    "lambda": float(row["lambdas"]),
                    "q"     : float(row["q"]),
                    "mass"  : float(row["MW"]),
                    "r0"    : float(row["bondlength"]),
                }

    # ------------------------------------------------------------------
    # CGForceField interface
    # ------------------------------------------------------------------

    def add_masses(self, system: mm.System, chain_meta: List[dict]) -> None:
        for meta in chain_meta:
            for aa in meta["sequence"]:
                m = self._params.get(aa, {}).get("mass", 57.05)
                system.addParticle(m * unit.amu)

    def create_nonbonded_forces(
        self,
        topology: app.Topology,
        chain_meta: List[dict],
        temperature: float,
        ionic: float,
    ) -> List[mm.Force]:
        forces = []

        # --- Ashbaugh-Hatch (with ID parameter for molecule-type mixing) ---
        # ID parameter: protein=1, lipid=0, crowder=-1
        # Uses select(id1+id2, (id1*id2)*0.5*(l1+l2), fixed_lambda) for mixing
        eps = self.EPS_LJ
        rc  = self.RC_LJ
        # fixed_lambda used for non-protein interactions (lipid/lipid, protein/lipid, etc.)
        fixed_lambda = 0.5  # Default for cross-type interactions
        expr = (
            f"{eps}*select(step(r-2^(1/6)*s),"
            f"4*l*((s/r)^12-(s/r)^6-shift),"
            f"4*((s/r)^12-(s/r)^6-l*shift)+(1-l));"
            f"l=select(id1+id2,(id1*id2)*0.5*(l1+l2),{fixed_lambda});"
            f"shift=(s/{rc})^12-(s/{rc})^6;"
            f"s=0.5*(s1+s2)"
        )
        ah = mm.CustomNonbondedForce(expr)
        ah.addPerParticleParameter("s")
        ah.addPerParticleParameter("l")
        ah.addPerParticleParameter("id")
        ah.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        ah.setCutoffDistance(rc * unit.nanometer)
        ah.setForceGroup(0)

        eps_yu, k_yu = debye_huckel_params(temperature, ionic)
        shift_yu = float(np.exp(-k_yu * self.RC_YU) / self.RC_YU)
        yu_expr = f"q*{eps_yu:.6f}*(exp(-{k_yu:.6f}*r)/r-{shift_yu:.6e}); q=q1*q2"
        yu = mm.CustomNonbondedForce(yu_expr)
        yu.addPerParticleParameter("q")
        yu.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        yu.setCutoffDistance(self.RC_YU * unit.nanometer)
        yu.setForceGroup(1)

        # Add per-particle parameters
        # ID parameter: protein=1 (for proper AH lambda mixing)
        for meta in chain_meta:
            for aa in meta["sequence"]:
                p = self._params.get(aa, self._params.get("G", {}))
                ah.addParticle([p["sigma"], p["lambda"], 1])  # id=1 for protein
                yu.addParticle([p["q"]])

        # Exclude bonded 1-2 pairs
        bonds = [(b[0].index, b[1].index) for b in topology.bonds()]
        ah.createExclusionsFromBonds(bonds, 1)
        yu.createExclusionsFromBonds(bonds, 1)

        forces.extend([ah, yu])
        return forces

    def build_harmonic_bonds(
        self,
        topology: app.Topology,
        r0: float = 0.38,
        k: float = None,
    ) -> mm.HarmonicBondForce:
        """Use per-residue bond length from the parameter CSV."""
        k = k or self.K_BOND
        hb = mm.HarmonicBondForce()
        hb.setUsesPeriodicBoundaryConditions(True)
        atoms = list(topology.atoms())
        for bond in topology.bonds():
            i_atom = bond[0]
            j_atom = bond[1]
            aa_i = self._aa_from_atom(i_atom)
            r0_i = self._params.get(aa_i, {}).get("r0", r0)
            hb.addBond(
                i_atom.index, j_atom.index,
                r0_i * unit.nanometer,
                k * unit.kilojoule_per_mole / unit.nanometer ** 2,
            )
        return hb

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _aa_from_atom(atom: app.topology.Atom) -> str:
        from ..topology import THREE_TO_ONE
        return THREE_TO_ONE.get(atom.residue.name, "G")
