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
    "ALA": {"mass": 71.079,  "charge": 0.0,  "radius": 0.2845, "epsilon": 0.29519131, "azero": 0.0002, "surface": 0.796},
    "ARG": {"mass": 157.197, "charge": 1.0,  "radius": 0.3567, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.921},
    "ASN": {"mass": 114.104, "charge": 0.0,  "radius": 0.3150, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.281},
    "ASP": {"mass": 114.080, "charge":-1.0,  "radius": 0.3114, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.162},
    "CYS": {"mass": 103.139, "charge": 0.0,  "radius": 0.3024, "epsilon": 0.29519131, "azero": 0.0002, "surface": 1.074},
    "GLN": {"mass": 128.131, "charge": 0.0,  "radius": 0.3311, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.575},
    "GLU": {"mass": 128.107, "charge":-1.0,  "radius": 0.3279, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.462},
    "GLY": {"mass": 57.052,  "charge": 0.0,  "radius": 0.2617, "epsilon": 0.29519131, "azero": 0.0002, "surface": 0.544},
    "HIS": {"mass": 137.142, "charge": 0.0,  "radius": 0.3338, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.634},
    "HSD": {"mass": 137.142, "charge": 0.0,  "radius": 0.3338, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.634},
    "HSE": {"mass": 137.142, "charge": 0.0,  "radius": 0.3338, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.634},
    "ILE": {"mass": 113.160, "charge": 0.0,  "radius": 0.3360, "epsilon": 0.29519131, "azero": 0.0002, "surface": 1.410},
    "LEU": {"mass": 113.160, "charge": 0.0,  "radius": 0.3363, "epsilon": 0.29519131, "azero": 0.0002, "surface": 1.519},
    "LYS": {"mass": 129.183, "charge": 1.0,  "radius": 0.3439, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.923},
    "MET": {"mass": 131.193, "charge": 0.0,  "radius": 0.3381, "epsilon": 0.29519131, "azero": 0.0002, "surface": 1.620},
    "PHE": {"mass": 147.177, "charge": 0.0,  "radius": 0.3556, "epsilon": 0.29519131, "azero": 0.0002, "surface": 1.869},
    "PRO": {"mass": 98.125,  "charge": 0.0,  "radius": 0.3187, "epsilon": 0.29519131, "azero": 0.0002, "surface": 0.974},
    "SER": {"mass": 87.078,  "charge": 0.0,  "radius": 0.2927, "epsilon": 0.17596101, "azero": 0.0,        "surface": 0.933},
    "THR": {"mass": 101.105, "charge": 0.0,  "radius": 0.3108, "epsilon": 0.17596101, "azero": 0.0,        "surface": 1.128},
    "TRP": {"mass": 186.214, "charge": 0.0,  "radius": 0.3754, "epsilon": 0.29519131, "azero": 0.0002, "surface": 2.227},
    "TYR": {"mass": 163.176, "charge": 0.0,  "radius": 0.3611, "epsilon": 0.29519131, "azero": 0.0002, "surface": 2.018},
    "VAL": {"mass": 99.133,  "charge": 0.0,  "radius": 0.3205, "epsilon": 0.29519131, "azero": 0.0002, "surface": 1.232},
}

_CATIONPI   = 0.30   # kJ/mol  cation-pi well depth
_PIPI       = 0.10   # kJ/mol  pi-pi well depth
_KBOND      = 4184.0 # kJ/mol/nm^2
_KANGLE_PRO = 4.184  # kJ/mol/rad^2
_L0_PRO     = 0.38   # nm
_THETA0     = np.pi  # 180 deg
_KAPPA      = 1.0    # nm  (Debye length for SASA electrostatics)
_SURF_REF   = 0.7    # COCOMO surface threshold
_CUTOFF     = 3.0    # nm  (nonbonded cutoff, matching original COCOMO)
_DEFAULT_SASA = 999.0  # nm^2 (default SASA for exposed residues, matching original COCOMO)


class CocomoFF(CGForceField):
    """
    COCOMO2 force field for protein condensates.

    SASA values are estimated from the initial CG structure using
    a simple neighbour-count model, consistent with Feig lab COCOMO.
    """
    
    # COCOMO-specific ENM parameters (matching original COCOMO implementation)
    ENM_K = 500.0       # kJ/mol/nm^2
    ENM_CUTOFF = 0.9    # nm

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
        # Build atom-to-chain mapping for IDP/MDP distinction
        atom_to_chain = {}
        for meta in chain_meta:
            for i in range(meta["start"], meta["end"]):
                atom_to_chain[i] = meta
        
        # Compute SASA using mdsim for MDP components (matching original COCOMO)
        # Build per-atom SASA array (IDP residues use default value)
        sasa_vals = self._compute_all_sasa(chain_meta)
        surface_factors = self._get_surface_factors(topology, sasa_vals, atom_to_chain)

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
        vdw.setCutoffDistance(_CUTOFF * unit.nanometer)

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
        elec.setCutoffDistance(_CUTOFF * unit.nanometer)

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
        cpi.setCutoffDistance(_CUTOFF * unit.nanometer)

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
        pipi.setCutoffDistance(_CUTOFF * unit.nanometer)

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
    # ENM for folded domains (COCOMO-specific parameters)
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
        
        COCOMO-specific implementation using original COCOMO parameters:
        - k = 500 kJ/mol/nm^2
        - cutoff = 0.9 nm
        
        Reference: Original COCOMO implementation (cocomo_model.py)
        """
        from typing import Optional
        
        # Use COCOMO-specific defaults if not provided
        if k is None:
            k = self.ENM_K
        if cutoff is None:
            cutoff = self.ENM_CUTOFF
        min_seq_sep = 3  # COCOMO: skip i+1 and i+2, start from i+3 (matching original implementation)
        
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
    # SASA estimation
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_sasa(pdb_path: str) -> np.ndarray:
        """
        Compute per-residue SASA using internal Shrake-Rupley implementation.
        
        Uses n_sphere_points=1920 (matching original COCOMO implementation).
        Requires full-atom PDB (not CG CA-only structure).
        
        Args:
            pdb_path: Path to full-atom PDB file
            
        Returns:
            SASA values in nm^2 per residue, or None if calculation fails
        """
        try:
            from ..sasa import calc_sasa_from_pdb
            return calc_sasa_from_pdb(pdb_path, n_sphere_points=1920)
        except Exception:
            return None

    @staticmethod
    def _compute_all_sasa(chain_meta: List[dict]) -> Optional[np.ndarray]:
        """
        Compute SASA values for all atoms across all chains.
        
        For MDP components with pdb_path, use mdsim to compute SASA from full-atom PDB.
        For IDP components or when SASA computation fails, use default value.
        
        Returns:
            Array of SASA values (nm^2) per atom, or None if no SASA computed
        """
        from ..config import ComponentType
        
        all_sasa = []
        has_sasa = False
        
        for meta in chain_meta:
            n_res = len(meta["sequence"])
            pdb_path = meta.get("pdb_path")
            comp_type = meta.get("comp_type")
            
            if comp_type == ComponentType.MDP and pdb_path:
                # Try to compute SASA from full-atom PDB
                sasa = CocomoFF._compute_sasa(pdb_path)
                if sasa is not None and len(sasa) == n_res:
                    all_sasa.extend(sasa)
                    has_sasa = True
                else:
                    # SASA computation failed or length mismatch, use default
                    all_sasa.extend([_DEFAULT_SASA] * n_res)
            else:
                # IDP or no pdb_path, use default SASA
                all_sasa.extend([_DEFAULT_SASA] * n_res)
        
        return np.array(all_sasa) if has_sasa else None

    @staticmethod
    def _get_surface_factors(
        topology: app.Topology,
        sasa_vals: Optional[np.ndarray],
        atom_to_chain: Optional[dict] = None,
        surf_ref: float = _SURF_REF,
    ) -> List[float]:
        """
        Convert SASA values to COCOMO surface screening factors.

        surface_factor = min(1, sasa_nm / (surf_ref * residue_native_surface))
        
        For IDP (disordered): fully exposed, surface_factor = 1.0
        For MDP (folded domains): computed from SASA or defaults to 0.7 (global surfscale)
        """
        from ..config import ComponentType
        
        factors = []
        for atom_idx, atom in enumerate(topology.atoms()):
            p = _FF_PARAM.get(atom.residue.name, _FF_PARAM["GLY"])
            native_sasa = p["surface"]   # nm^2
            
            if sasa_vals is not None:
                # Use provided SASA values (typically for MDP from PDB)
                s_nm = float(sasa_vals[atom_idx]) / 100.0   # Å^2 -> nm^2
                sf = min(1.0, s_nm / native_sasa / surf_ref)
            else:
                # No SASA provided: distinguish IDP vs MDP
                if atom_to_chain and atom_idx in atom_to_chain:
                    comp_type = atom_to_chain[atom_idx].get("comp_type")
                    if comp_type == ComponentType.IDP:
                        # IDP: fully exposed (surface_factor = 1.0, no screening)
                        sf = 1.0
                    else:
                        # MDP: use default surfscale (0.7)
                        sf = surf_ref
                else:
                    # Fallback: use surf_ref as default
                    sf = surf_ref
            factors.append(sf)
        return factors
