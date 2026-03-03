"""
COCOMO2 force field.

Physics:
  - Harmonic CA-CA bonds
  - HarmonicAngleForce (1-2-3 triplets, theta0 = 180 deg)
  - SASA-screened van der Waals (10-5 LJ)
  - SASA-screened electrostatics (Debye-screened, distance-based)
  - Cation-pi and pi-pi interactions (interaction group sub-sets)
  - Go-like ENM for folded domains

Reference: Gallivan et al. J. Chem. Theory Comput. 2023; Feig Lab COCOMO.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit

from .base import CGForceField


# ---------------------------------------------------------------------------
# Per-residue force-field parameters for COCOMO2
# ---------------------------------------------------------------------------

_FF_PARAM: Dict[str, dict] = {
    "ALA": {"mass": 71.079,  "charge": 0.0,  "radius": 0.2845, "epsilon": 0.295, "azero": 2e-4, "surface": 0.796},
    "ARG": {"mass": 157.197, "charge": 1.0,  "radius": 0.3567, "epsilon": 0.176, "azero": 0.0,  "surface": 1.921},
    "ASN": {"mass": 114.104, "charge": 0.0,  "radius": 0.3150, "epsilon": 0.176, "azero": 0.0,  "surface": 1.281},
    "ASP": {"mass": 114.080, "charge":-1.0,  "radius": 0.3114, "epsilon": 0.176, "azero": 0.0,  "surface": 1.162},
    "CYS": {"mass": 103.139, "charge": 0.0,  "radius": 0.3024, "epsilon": 0.295, "azero": 2e-4, "surface": 1.074},
    "GLN": {"mass": 128.131, "charge": 0.0,  "radius": 0.3311, "epsilon": 0.176, "azero": 0.0,  "surface": 1.575},
    "GLU": {"mass": 128.107, "charge":-1.0,  "radius": 0.3279, "epsilon": 0.176, "azero": 0.0,  "surface": 1.462},
    "GLY": {"mass": 57.052,  "charge": 0.0,  "radius": 0.2617, "epsilon": 0.295, "azero": 2e-4, "surface": 0.544},
    "HIS": {"mass": 137.142, "charge": 0.0,  "radius": 0.3338, "epsilon": 0.176, "azero": 0.0,  "surface": 1.634},
    "ILE": {"mass": 113.160, "charge": 0.0,  "radius": 0.3411, "epsilon": 0.295, "azero": 2e-4, "surface": 1.671},
    "LEU": {"mass": 113.160, "charge": 0.0,  "radius": 0.3411, "epsilon": 0.295, "azero": 2e-4, "surface": 1.671},
    "LYS": {"mass": 128.174, "charge": 1.0,  "radius": 0.3478, "epsilon": 0.176, "azero": 0.0,  "surface": 1.771},
    "MET": {"mass": 131.197, "charge": 0.0,  "radius": 0.3459, "epsilon": 0.295, "azero": 2e-4, "surface": 1.724},
    "PHE": {"mass": 147.177, "charge": 0.0,  "radius": 0.3544, "epsilon": 0.295, "azero": 2e-4, "surface": 1.875},
    "PRO": {"mass": 97.116,  "charge": 0.0,  "radius": 0.3146, "epsilon": 0.295, "azero": 2e-4, "surface": 1.254},
    "SER": {"mass": 87.078,  "charge": 0.0,  "radius": 0.2868, "epsilon": 0.176, "azero": 0.0,  "surface": 0.865},
    "THR": {"mass": 101.105, "charge": 0.0,  "radius": 0.3063, "epsilon": 0.176, "azero": 0.0,  "surface": 1.172},
    "TRP": {"mass": 186.213, "charge": 0.0,  "radius": 0.3827, "epsilon": 0.295, "azero": 2e-4, "surface": 2.376},
    "TYR": {"mass": 163.176, "charge": 0.0,  "radius": 0.3648, "epsilon": 0.176, "azero": 0.0,  "surface": 2.011},
    "VAL": {"mass": 99.133,  "charge": 0.0,  "radius": 0.3275, "epsilon": 0.295, "azero": 2e-4, "surface": 1.416},
}

_CATIONPI   = 0.30   # kJ/mol  cation-pi well depth
_PIPI       = 0.10   # kJ/mol  pi-pi well depth
_KBOND      = 4184.0 # kJ/mol/nm^2
_KANGLE_PRO = 4.184  # kJ/mol/rad^2
_L0_PRO     = 0.38   # nm
_THETA0     = np.pi  # 180 deg
_KAPPA      = 1.0    # nm  (Debye length for SASA electrostatics)
_SURF_REF   = 0.7    # COCOMO surface threshold


class CocomoFF(CGForceField):
    """
    COCOMO2 force field for protein condensates.

    SASA values are estimated from the initial CG structure using
    a simple neighbour-count model, consistent with Feig lab COCOMO.
    """

    def add_masses(self, system: mm.System, chain_meta: List[dict]) -> None:
        for meta in chain_meta:
            for aa in meta["sequence"]:
                from ..topology import ONE_TO_THREE
                three = ONE_TO_THREE.get(aa, "GLY")
                m = _FF_PARAM.get(three, _FF_PARAM["GLY"])["mass"]
                system.addParticle(m * unit.amu)

    def create_nonbonded_forces(
        self,
        topology: app.Topology,
        chain_meta: List[dict],
        temperature: float,
        ionic: float,
        positions: Optional[np.ndarray] = None,
    ) -> List[mm.Force]:
        forces = []
        bonds = [(b[0].index, b[1].index) for b in topology.bonds()]

        # Per-atom SASA-based surface screening factors
        sasa_vals = self._compute_sasa(topology, positions) if positions is not None else None
        surface_factors = self._get_surface_factors(topology, sasa_vals)

        # --- 10-5 LJ van der Waals ---
        vdw_expr = (
            "S*4*eps*((sigma/r)^10-(sigma/r)^5);"
            "sigma=0.5*(sigma1+sigma2);"
            "eps=sqrt(epsilon1*epsilon2);"
            "S=sqrt(S1*S2)"
        )
        vdw = mm.CustomNonbondedForce(vdw_expr)
        vdw.addPerParticleParameter("sigma")
        vdw.addPerParticleParameter("epsilon")
        vdw.addPerParticleParameter("S")
        vdw.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        vdw.setCutoffDistance(3.0 * unit.nanometer)

        # --- Electrostatics ---
        elec_expr = (
            "S*(A+Z)/r*exp(-r/K0);"
            "A=A1*A2;"
            "Z=Z1+Z2;"
            "S=sqrt(S1*S2)"
        )
        elec = mm.CustomNonbondedForce(elec_expr)
        elec.addGlobalParameter("K0", _KAPPA * unit.nanometer)
        elec.addPerParticleParameter("A")
        elec.addPerParticleParameter("Z")
        elec.addPerParticleParameter("S")
        elec.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        elec.setCutoffDistance(3.0 * unit.nanometer)

        # Cation-pi force (same form as VdW, different epsilon, interaction group)
        cpi_expr = (
            "S*4*{eps}*((sigma/r)^10-(sigma/r)^5);"
            "sigma=0.5*(sigma1+sigma2);"
            "S=sqrt(S1*S2)"
        ).format(eps=_CATIONPI)
        cpi = mm.CustomNonbondedForce(cpi_expr)
        cpi.addPerParticleParameter("sigma")
        cpi.addPerParticleParameter("S")
        cpi.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        cpi.setCutoffDistance(3.0 * unit.nanometer)

        # Pi-pi force
        pipi_expr = (
            "S*4*{eps}*((sigma/r)^10-(sigma/r)^5);"
            "sigma=0.5*(sigma1+sigma2);"
            "S=sqrt(S1*S2)"
        ).format(eps=_PIPI)
        pipi = mm.CustomNonbondedForce(pipi_expr)
        pipi.addPerParticleParameter("sigma")
        pipi.addPerParticleParameter("S")
        pipi.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        pipi.setCutoffDistance(3.0 * unit.nanometer)

        cation_idx, aromatic_idx = [], []

        for atom_idx, atom in enumerate(topology.atoms()):
            p = _FF_PARAM.get(atom.residue.name, _FF_PARAM["GLY"])
            sf = surface_factors[atom_idx]

            # VdW
            sigma_vdw = p["radius"] * 2 * 2 ** (-1.0 / 6.0)
            vdw.addParticle([sigma_vdw * unit.nanometer,
                             p["epsilon"] * unit.kilojoule_per_mole,
                             sf])

            # Electrostatics
            A_val = (np.sqrt(0.75 * abs(p["charge"])) * np.sign(p["charge"])
                     if p["charge"] != 0 else 0.0)
            elec.addParticle([A_val * unit.nanometer * unit.kilojoule_per_mole,
                              p["azero"] * (unit.nanometer * unit.kilojoule_per_mole) ** 0.5,
                              sf])

            # Cation / aromatic tracking
            sigma_sp = p["radius"] * 2 * 2 ** (-1.0 / 6.0)
            cpi.addParticle([sigma_sp * unit.nanometer, sf])
            pipi.addParticle([sigma_sp * unit.nanometer, sf])

            if atom.residue.name in ("ARG", "LYS"):
                cation_idx.append(atom_idx)
            if atom.residue.name in ("PHE", "TRP", "TYR"):
                aromatic_idx.append(atom_idx)

        vdw.createExclusionsFromBonds(bonds, 1)
        elec.createExclusionsFromBonds(bonds, 1)

        forces.extend([vdw, elec])

        if cation_idx and aromatic_idx:
            cpi.createExclusionsFromBonds(bonds, 1)
            cpi.addInteractionGroup(cation_idx, aromatic_idx)
            forces.append(cpi)

        if aromatic_idx:
            pipi.createExclusionsFromBonds(bonds, 1)
            pipi.addInteractionGroup(aromatic_idx, aromatic_idx)
            forces.append(pipi)

        return forces

    def build_angle_force(self, topology: app.Topology) -> mm.HarmonicAngleForce:
        """1-2-3 harmonic angle force for backbone stiffness (theta0 = 180 deg)."""
        ha = mm.HarmonicAngleForce()
        ha.setUsesPeriodicBoundaryConditions(True)
        # Collect per-chain atom lists
        for chain in topology.chains():
            atoms = list(chain.atoms())
            for i in range(len(atoms) - 2):
                ha.addAngle(
                    atoms[i].index, atoms[i + 1].index, atoms[i + 2].index,
                    _THETA0 * unit.radian,
                    _KANGLE_PRO * unit.kilojoule_per_mole / unit.radian ** 2,
                )
        return ha

    def build_harmonic_bonds(
        self,
        topology: app.Topology,
        r0: float = _L0_PRO,
        k: float = _KBOND,
    ) -> mm.HarmonicBondForce:
        return super().build_harmonic_bonds(topology, r0=r0, k=k)

    # ------------------------------------------------------------------
    # SASA estimation
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_sasa(topology: app.Topology, positions: np.ndarray) -> np.ndarray:
        """
        Estimate per-residue solvent accessible surface area (nm^2) using
        MDAnalysis shrake-rupley, or a uniform default if unavailable.
        """
        try:
            import MDAnalysis as mda
            from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
        except ImportError:
            return None

        try:
            import MDAnalysis as mda
            from MDAnalysis.analysis import shrake_rupley as sr
            # Build a minimal MDA universe from positions
            n = positions.shape[0]
            u = mda.Universe.empty(n, n_residues=n, n_segments=1,
                                   atom_resindex=np.arange(n),
                                   residue_segindex=np.zeros(n, dtype=int))
            u.add_TopologyAttr("name", ["CA"] * n)
            u.add_TopologyAttr("mass", [57.0] * n)
            u.atoms.positions = positions * 10.0   # nm -> Å
            srk = sr.ShrakeRupley(u, probe_radius=0.14, n_sphere_points=100)
            srk.run()
            return srk.results.areas[0]   # Å^2 per residue
        except Exception:
            return None

    @staticmethod
    def _get_surface_factors(
        topology: app.Topology,
        sasa_vals: Optional[np.ndarray],
        surf_ref: float = _SURF_REF,
    ) -> List[float]:
        """
        Convert SASA values to COCOMO surface screening factors.

        surface_factor = min(1, sasa_nm / (surf_ref * residue_native_surface))
        """
        factors = []
        for atom_idx, atom in enumerate(topology.atoms()):
            p = _FF_PARAM.get(atom.residue.name, _FF_PARAM["GLY"])
            native_sasa = p["surface"]   # nm^2
            if sasa_vals is not None:
                s_nm = float(sasa_vals[atom_idx]) / 100.0   # Å^2 -> nm^2
                sf = min(1.0, s_nm / (surf_ref * native_sasa + 1e-9))
            else:
                sf = 0.7   # default surface factor
            factors.append(sf)
        return factors
