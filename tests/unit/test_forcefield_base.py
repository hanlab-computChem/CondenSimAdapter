"""
Unit tests for core/forcefield/base.py

Tests cover:
- CGForceField abstract interface
- Shared helper methods (harmonic bonds, ENM, droplet)
- Debye-Huckel parameter calculation
"""

from __future__ import annotations

import numpy as np
import pytest

# Skip all tests if OpenMM is not available
pytest.importorskip("openmm")

import openmm as mm
import openmm.app as app
import openmm.unit as unit

from CondenSimAdapter.core.forcefield.base import (
    CGForceField,
    debye_huckel_params,
)


# =============================================================================
# Concrete test implementation of abstract base class
# =============================================================================

class TestForceField(CGForceField):
    """Concrete implementation for testing base class methods."""

    def create_nonbonded_forces(
        self,
        topology: app.Topology,
        chain_meta: list,
        temperature: float,
        ionic: float,
    ) -> list:
        """Return empty list for testing."""
        return []


# =============================================================================
# Debye-Huckel Tests
# =============================================================================

class TestDebyeHuckelParams:
    """Tests for debye_huckel_params function."""

    def test_returns_tuple(self):
        """Test function returns two values."""
        eps_yu, k_yu = debye_huckel_params(300.0, 0.15)
        assert isinstance(eps_yu, (float, np.floating))
        assert isinstance(k_yu, (float, np.floating))

    def test_eps_yu_positive(self):
        """Test energy prefactor is positive."""
        eps_yu, _ = debye_huckel_params(300.0, 0.15)
        assert eps_yu > 0

    def test_k_yu_positive(self):
        """Test inverse Debye length is positive."""
        _, k_yu = debye_huckel_params(300.0, 0.15)
        assert k_yu > 0

    def test_temperature_dependence(self):
        """Test higher temperature increases eps_yu."""
        eps_yu_300, _ = debye_huckel_params(300.0, 0.15)
        eps_yu_350, _ = debye_huckel_params(350.0, 0.15)
        assert eps_yu_350 > eps_yu_300

    def test_ionic_dependence(self):
        """Test higher ionic strength increases k_yu."""
        _, k_yu_low = debye_huckel_params(300.0, 0.05)
        _, k_yu_high = debye_huckel_params(300.0, 0.5)
        assert k_yu_high > k_yu_low


# =============================================================================
# Build Harmonic Bonds Tests
# =============================================================================

class TestBuildHarmonicBonds:
    """Tests for build_harmonic_bonds method."""

    @pytest.fixture
    def simple_topology(self) -> app.Topology:
        """Create a simple 3-atom topology with 2 bonds."""
        top = app.Topology()
        chain = top.addChain()
        res = top.addResidue("ALA", chain)
        elem = app.element.carbon
        
        atom1 = top.addAtom("CA", elem, res)
        atom2 = top.addAtom("CA", elem, res)
        atom3 = top.addAtom("CA", elem, res)
        
        top.addBond(atom1, atom2)
        top.addBond(atom2, atom3)
        
        return top

    def test_returns_harmonic_bond_force(self, simple_topology):
        """Test method returns HarmonicBondForce."""
        ff = TestForceField()
        force = ff.build_harmonic_bonds(simple_topology)
        assert isinstance(force, mm.HarmonicBondForce)

    def test_correct_number_of_bonds(self, simple_topology):
        """Test force contains correct number of bonds."""
        ff = TestForceField()
        force = ff.build_harmonic_bonds(simple_topology)
        assert force.getNumBonds() == 2

    def test_default_parameters(self, simple_topology):
        """Test default bond parameters."""
        ff = TestForceField()
        force = ff.build_harmonic_bonds(simple_topology)
        
        # Get first bond parameters
        atom1, atom2, length, k = force.getBondParameters(0)
        assert length.value_in_unit(unit.nanometer) == pytest.approx(0.38, abs=1e-6)
        assert k.value_in_unit(unit.kilojoule_per_mole / unit.nanometer ** 2) == pytest.approx(8368.0, abs=1e-6)

    def test_custom_parameters(self, simple_topology):
        """Test custom bond parameters."""
        ff = TestForceField()
        force = ff.build_harmonic_bonds(simple_topology, r0=0.40, k=10000.0)
        
        atom1, atom2, length, k = force.getBondParameters(0)
        assert length.value_in_unit(unit.nanometer) == pytest.approx(0.40, abs=1e-6)
        assert k.value_in_unit(unit.kilojoule_per_mole / unit.nanometer ** 2) == pytest.approx(10000.0, abs=1e-6)

    def test_periodic_boundary_conditions(self, simple_topology):
        """Test PBC is enabled."""
        ff = TestForceField()
        force = ff.build_harmonic_bonds(simple_topology)
        assert force.usesPeriodicBoundaryConditions()


# =============================================================================
# Build ENM Bonds Tests
# =============================================================================

class TestBuildENMBonds:
    """Tests for build_enm_bonds method."""

    @pytest.fixture
    def chain_meta_single_domain(self):
        """Chain metadata with single folded domain."""
        # 10 atoms, folded domain from 2-5 (1-based)
        return [{
            "name": "test",
            "start": 0,
            "end": 10,
            "sequence": "AAAAAAAAAA",
            "folded_domains": [(2, 5)],  # atoms 1,2,3,4 in 0-based
        }]

    @pytest.fixture
    def positions_10_atoms(self):
        """10 atoms in a line, 0.38 nm apart."""
        pos = np.zeros((10, 3))
        pos[:, 2] = np.arange(10) * 0.38
        return pos

    def test_no_domains_returns_none(self, positions_10_atoms):
        """Test that chain with no folded domains returns None."""
        ff = TestForceField()
        chain_meta = [{
            "name": "test",
            "start": 0,
            "end": 10,
            "sequence": "AAAAAAAAAA",
            "folded_domains": [],
        }]
        result = ff.build_enm_bonds(positions_10_atoms, chain_meta)
        assert result is None

    def test_harmonic_mode(self, chain_meta_single_domain, positions_10_atoms):
        """Test harmonic ENM mode."""
        ff = TestForceField()
        force = ff.build_enm_bonds(
            positions_10_atoms,
            chain_meta_single_domain,
            restraint_type="harmonic",
            k=700.0,
            cutoff=1.0,
        )
        assert isinstance(force, mm.HarmonicBondForce)
        assert force.getNumBonds() > 0

    def test_go_mode(self, chain_meta_single_domain, positions_10_atoms):
        """Test Go-like ENM mode."""
        ff = TestForceField()
        force = ff.build_enm_bonds(
            positions_10_atoms,
            chain_meta_single_domain,
            restraint_type="go",
            k=1.0,
            cutoff=1.0,
        )
        assert isinstance(force, mm.CustomBondForce)

    def test_respects_cutoff(self, positions_10_atoms):
        """Test that ENM only adds bonds within cutoff."""
        ff = TestForceField()
        # Domain with atoms close together in a cluster
        chain_meta = [{
            "name": "test",
            "start": 0,
            "end": 10,
            "sequence": "AAAAAAAAAA",
            "folded_domains": [(1, 10)],  # All atoms
        }]
        
        # Create positions where some pairs are within cutoff and some are not
        # Place atoms in a compact cluster
        positions_cluster = np.random.rand(10, 3) * 0.5  # Random within 0.5 nm
        
        # Small cutoff: should have few bonds
        force_small = ff.build_enm_bonds(
            positions_cluster, chain_meta, cutoff=0.2
        )
        n_bonds_small = force_small.getNumBonds()
        
        # Large cutoff: should have more bonds
        force_large = ff.build_enm_bonds(
            positions_cluster, chain_meta, cutoff=1.0
        )
        n_bonds_large = force_large.getNumBonds()
        
        assert n_bonds_large >= n_bonds_small


# =============================================================================
# Build Droplet Force Tests
# =============================================================================

class TestBuildDropletForce:
    """Tests for build_droplet_force method."""

    @pytest.fixture
    def simple_topology_5atoms(self) -> app.Topology:
        """Create a simple 5-atom topology."""
        top = app.Topology()
        chain = top.addChain()
        res = top.addResidue("ALA", chain)
        elem = app.element.carbon
        
        for i in range(5):
            top.addAtom(f"CA{i}", elem, res)
        
        return top

    @pytest.fixture
    def positions_5atoms(self):
        """5 atoms centered at origin."""
        return np.array([
            [0, 0, 0],
            [0.5, 0, 0],
            [-0.5, 0, 0],
            [0, 0.5, 0],
            [0, -0.5, 0],
        ], dtype=np.float64)

    def test_returns_custom_external_force(self, simple_topology_5atoms, positions_5atoms):
        """Test method returns CustomExternalForce."""
        ff = TestForceField()
        force = ff.build_droplet_force(simple_topology_5atoms, positions_5atoms, radius=10.0)
        assert isinstance(force, mm.CustomExternalForce)

    def test_correct_number_of_particles(self, simple_topology_5atoms, positions_5atoms):
        """Test force acts on correct number of particles."""
        ff = TestForceField()
        force = ff.build_droplet_force(
            simple_topology_5atoms, positions_5atoms, radius=10.0, stride=1
        )
        # With stride=1, all 5 atoms should have the force
        assert force.getNumParticles() == 5

    def test_stride_works(self, simple_topology_5atoms, positions_5atoms):
        """Test stride parameter reduces number of particles."""
        ff = TestForceField()
        force = ff.build_droplet_force(
            simple_topology_5atoms, positions_5atoms, radius=10.0, stride=2
        )
        # With stride=2, should have particles 0, 2, 4
        assert force.getNumParticles() == 3

    def test_global_parameters_set(self, simple_topology_5atoms, positions_5atoms):
        """Test global parameters are set correctly."""
        ff = TestForceField()
        force = ff.build_droplet_force(
            simple_topology_5atoms, positions_5atoms, radius=5.0, k=2.0
        )
        
        # Check global parameters by index
        k_drop = force.getGlobalParameterDefaultValue(0)  # k_drop is first
        r_drop = force.getGlobalParameterDefaultValue(1)  # r_drop is second
        
        assert k_drop == pytest.approx(2.0, abs=1e-6)
        assert r_drop == pytest.approx(5.0, abs=1e-6)


# =============================================================================
# Add Masses Tests
# =============================================================================

class TestAddMasses:
    """Tests for add_masses method."""

    def test_adds_correct_number_of_masses(self):
        """Test correct number of masses are added to system."""
        ff = TestForceField()
        system = mm.System()
        
        chain_meta = [
            {"name": "A", "sequence": "AAA"},
            {"name": "B", "sequence": "BB"},
        ]
        
        ff.add_masses(system, chain_meta)
        assert system.getNumParticles() == 5  # 3 + 2

    def test_default_mass_values(self):
        """Test default mass values from RESIDUE_MASS table."""
        ff = TestForceField()
        system = mm.System()
        
        chain_meta = [{"name": "test", "sequence": "AG"}]
        
        ff.add_masses(system, chain_meta)
        
        from CondenSimAdapter.core.topology import RESIDUE_MASS
        assert system.getParticleMass(0).value_in_unit(unit.amu) == pytest.approx(RESIDUE_MASS["A"])
        assert system.getParticleMass(1).value_in_unit(unit.amu) == pytest.approx(RESIDUE_MASS["G"])
