"""
Abstract base class for all CG protein force fields.

Shared helpers (harmonic bonds, ENM, droplet confinement) live here so
concrete force-field classes only implement their non-bonded physics.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Optional

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit


class CGForceField(ABC):
    """
    Protocol for a CG force field.

    Subclasses implement `create_nonbonded_forces()` and optionally
    override `add_masses()`.
    """

    # Whether ENM restraint pairs should be added to the exclusion list of all
    # CustomNonbondedForce objects (VdW, electrostatics).
    #
    # - True (default): matches OpenMpipi (topology.addBond + createExclusionsFromBonds)
    #   and CALVADOS (explicit add_exclusions for every restrained pair).
    # - False: matches original COCOMO, which intentionally allows the 10-5 LJ to
    #   act on ENM pairs (small attractive contribution at native distances deepens
    #   the native-contact potential well).
    _exclude_enm_from_nonbonded: bool = True

    @abstractmethod
    def create_nonbonded_forces(
        self,
        topology: app.Topology,
        chain_meta: List[dict],
        temperature: float,
        ionic: float,
    ) -> List[mm.Force]:
        """
        Return a list of OpenMM Force objects for non-bonded interactions.

        Args:
            topology:    OpenMM Topology (CA-only CG system).
            chain_meta:  Per-chain metadata from molecule.build_all_chains().
            temperature: Simulation temperature in K.
            ionic:       Ionic strength in M.
        """

    def add_masses(self, system: mm.System, chain_meta: List[dict]) -> None:
        """Set per-particle masses. Override to change masses."""
        from ..topology import RESIDUE_MASS

        for meta in chain_meta:
            for aa in meta["sequence"]:
                system.addParticle(RESIDUE_MASS.get(aa, 57.05) * unit.amu)

    # ------------------------------------------------------------------
    # Shared bonded-force builders (used by all sub-classes)
    # ------------------------------------------------------------------

    def build_harmonic_bonds(
        self,
        topology: app.Topology,
        r0: float = 0.38,
        k: float = 8368.0,
    ) -> mm.HarmonicBondForce:
        """
        CA-CA harmonic bonds for the backbone.

        Args:
            r0: Equilibrium bond length in nm.
            k:  Spring constant in kJ/mol/nm^2 (default = 8368 ~ HPS/Mpipi).
        """
        hb = mm.HarmonicBondForce()
        hb.setUsesPeriodicBoundaryConditions(True)
        for bond in topology.bonds():
            i, j = bond[0].index, bond[1].index
            hb.addBond(i, j, r0 * unit.nanometer, k * unit.kilojoule_per_mole / unit.nanometer**2)
        return hb

    def build_enm_bonds(
        self,
        positions: np.ndarray,
        chain_meta: List[dict],
        restraint_type: str = "harmonic",
        k: float = 700.0,
        cutoff: float = 0.9,
    ) -> Optional[mm.Force]:
        """
        Elastic Network Model (ENM) for folded domains.

        Adds harmonic bonds between all CA pairs within `cutoff` nm that
        belong to folded-domain segments.

        Args:
            positions:      (N, 3) nm array of initial CA coordinates.
            chain_meta:     Per-chain metadata with 'folded_domains'.
            restraint_type: 'harmonic' or 'go' (Go-like 10-12 potential).
            k:              Spring constant kJ/mol/nm^2 (harmonic) or well depth. Default 700.0 (CALVADOS standard).
            cutoff:         Contact distance cutoff in nm. Default 0.9 nm (CALVADOS standard).

        Returns:
            Force object or None if no folded domains exist.
        """
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
            for dom_s, dom_e in meta["folded_domains"]:
                a0 = chain_start + dom_s - 1  # absolute, 0-based
                a1 = chain_start + dom_e  # exclusive
                indices = list(range(a0, a1))
                for ii in range(len(indices)):
                    for jj in range(ii + 2, len(indices)):
                        gi, gj = indices[ii], indices[jj]
                        d = float(np.linalg.norm(positions[gi] - positions[gj]))
                        if d <= cutoff:
                            if restraint_type == "go":
                                cs.addBond(
                                    gi, gj, [d * unit.nanometer, k * unit.kilojoule_per_mole]
                                )
                            else:
                                cs.addBond(
                                    gi,
                                    gj,
                                    d * unit.nanometer,
                                    k * unit.kilojoule_per_mole / unit.nanometer**2,
                                )
                            n_bonds += 1

        return cs if n_bonds > 0 else None

    def build_droplet_force(
        self,
        topology: app.Topology,
        positions: np.ndarray,
        radius: float,
        k: float = 1.0,
        stride: int = 10,
    ) -> mm.CustomExternalForce:
        """
        Harmonic confinement force to keep chains inside a spherical droplet.

        The force is applied every `stride` residues to save compute.

        Energy: k * max(0, r - radius)^2
        where r = distance from box centre.
        """
        n_atoms = positions.shape[0]
        centre = positions.mean(axis=0)

        expr = (
            "k_drop * step(dist - r_drop) * (dist - r_drop)^2;"
            "dist = sqrt((x - cx)^2 + (y - cy)^2 + (z - cz)^2)"
        )
        force = mm.CustomExternalForce(expr)
        force.addGlobalParameter("k_drop", k * unit.kilojoule_per_mole / unit.nanometer**2)
        force.addGlobalParameter("r_drop", radius * unit.nanometer)
        force.addGlobalParameter("cx", float(centre[0]) * unit.nanometer)
        force.addGlobalParameter("cy", float(centre[1]) * unit.nanometer)
        force.addGlobalParameter("cz", float(centre[2]) * unit.nanometer)

        for i in range(0, n_atoms, stride):
            force.addParticle(i, [])
        return force


# ---------------------------------------------------------------------------
# Shared physics helpers
# ---------------------------------------------------------------------------


def debye_huckel_params(temperature: float, ionic: float):
    """
    Compute Yukawa / Debye-Hückel parameters.

    Returns:
        eps_yu: energy prefactor (kJ/mol · nm)
        k_yu:   inverse Debye length (nm^-1)
    """
    kT = 8.3145 * temperature * 1e-3  # kJ/mol

    def fepsw(T):
        return 5321 / T + 233.76 - 0.9297 * T + 1.417e-3 * T**2 - 8.292e-7 * T**3

    epsw = fepsw(temperature)
    lB = (1.6021766**2 / (4 * np.pi * 8.854188 * epsw)) * 6.02214076e3 / kT
    eps_yu = lB * kT
    k_yu = np.sqrt(8 * np.pi * lB * ionic * 6.02214076 / 10.0)
    return eps_yu, k_yu
