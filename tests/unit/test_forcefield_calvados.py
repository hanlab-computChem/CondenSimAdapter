"""
Unit tests for core/forcefield/calvados.py

Tests cover:
- CALVADOS2 and CALVADOS3 parameter loading
- Nonbonded force creation (Ashbaugh-Hatch + Yukawa)
- Per-residue bond lengths
- Version-specific behaviors
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("openmm")

import openmm as mm
import openmm.app as app
import openmm.unit as unit

from CondenSimAdapter.core.forcefield.calvados import CalvadosFF


# =============================================================================
# Initialization Tests
# =============================================================================

class TestCalvadosInitialization:
    """Tests for CalvadosFF initialization."""

    def test_default_version_is_2(self):
        """Test default version is 2."""
        ff = CalvadosFF()
        assert ff.version == 2

    def test_explicit_version_2(self):
        """Test explicit version 2."""
        ff = CalvadosFF(version=2)
        assert ff.version == 2

    def test_explicit_version_3(self):
        """Test explicit version 3."""
        ff = CalvadosFF(version=3)
        assert ff.version == 3

    def test_invalid_version_raises(self):
        """Test invalid version raises ValueError."""
        with pytest.raises(ValueError, match="version must be 2 or 3"):
            CalvadosFF(version=1)
        with pytest.raises(ValueError, match="version must be 2 or 3"):
            CalvadosFF(version=4)

    def test_parameters_loaded(self):
        """Test that residue parameters are loaded."""
        ff = CalvadosFF(version=2)
        assert len(ff._params) == 20  # 20 amino acids
        assert "A" in ff._params
        assert "G" in ff._params
        assert "R" in ff._params

    def test_parameter_keys(self):
        """Test that each residue has expected parameter keys."""
        ff = CalvadosFF(version=2)
        for aa, params in ff._params.items():
            assert "sigma" in params
            assert "lambda" in params
            assert "q" in params
            assert "mass" in params
            assert "r0" in params


# =============================================================================
# Parameter Value Tests
# =============================================================================

class TestCalvadosParameters:
    """Tests for parameter values."""

    def test_sigma_positive(self):
        """Test all sigma values are positive."""
        ff = CalvadosFF(version=2)
        for params in ff._params.values():
            assert params["sigma"] > 0

    def test_lambda_range(self):
        """Test lambda values are in reasonable range."""
        ff = CalvadosFF(version=2)
        for params in ff._params.values():
            # Lambda should be between 0 and 1 for hydrophobicity scale
            assert 0 <= params["lambda"] <= 1

    def test_charge_values(self):
        """Test charge values are reasonable."""
        ff = CalvadosFF(version=2)
        for aa, params in ff._params.items():
            q = params["q"]
            # Charges should be -1, 0, or +1 for standard residues
            assert q in [-1.0, 0.0, 1.0], f"Unexpected charge for {aa}: {q}"

    def test_mass_positive(self):
        """Test all mass values are positive."""
        ff = CalvadosFF(version=2)
        for params in ff._params.values():
            assert params["mass"] > 0

    def test_r0_positive(self):
        """Test all bond length values are positive."""
        ff = CalvadosFF(version=2)
        for params in ff._params.values():
            assert params["r0"] > 0

    def test_version_2_vs_3_differences(self):
        """Test that v2 and v3 have some different parameters."""
        ff2 = CalvadosFF(version=2)
        ff3 = CalvadosFF(version=3)
        
        # At least some parameters should differ between versions
        differences = []
        for aa in ff2._params:
            for key in ["sigma", "lambda", "q", "mass", "r0"]:
                if ff2._params[aa][key] != ff3._params[aa][key]:
                    differences.append((aa, key))
        
        assert len(differences) > 0, "Expected some differences between v2 and v3"


# =============================================================================
# Add Masses Tests
# =============================================================================

class TestAddMasses:
    """Tests for add_masses method."""

    def test_adds_correct_masses(self):
        """Test that correct masses are added based on sequence."""
        ff = CalvadosFF(version=2)
        system = mm.System()
        
        chain_meta = [{"name": "test", "sequence": "AG"}]
        ff.add_masses(system, chain_meta)
        
        assert system.getNumParticles() == 2
        
        # Check masses match parameter values
        mass_a = system.getParticleMass(0).value_in_unit(unit.amu)
        mass_g = system.getParticleMass(1).value_in_unit(unit.amu)
        
        assert mass_a == pytest.approx(ff._params["A"]["mass"], abs=1e-6)
        assert mass_g == pytest.approx(ff._params["G"]["mass"], abs=1e-6)


# =============================================================================
# Nonbonded Forces Tests
# =============================================================================

class TestCreateNonbondedForces:
    """Tests for create_nonbonded_forces method."""

    @pytest.fixture
    def simple_topology(self) -> app.Topology:
        """Create a simple 3-residue topology."""
        top = app.Topology()
        chain = top.addChain()
        
        for res_name in ["ALA", "GLY", "SER"]:
            res = top.addResidue(res_name, chain)
            top.addAtom("CA", app.element.carbon, res)
        
        # Add bonds
        atoms = list(top.atoms())
        top.addBond(atoms[0], atoms[1])
        top.addBond(atoms[1], atoms[2])
        
        return top

    @pytest.fixture
    def chain_meta(self):
        """Simple chain metadata."""
        return [{
            "name": "test",
            "sequence": "AGS",
        }]

    def test_returns_two_forces(self, simple_topology, chain_meta):
        """Test that two forces are returned (AH + Yukawa)."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        assert len(forces) == 2

    def test_first_force_is_ah(self, simple_topology, chain_meta):
        """Test first force is Ashbaugh-Hatch."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        assert isinstance(forces[0], mm.CustomNonbondedForce)

    def test_second_force_is_yukawa(self, simple_topology, chain_meta):
        """Test second force is Yukawa."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        assert isinstance(forces[1], mm.CustomNonbondedForce)

    def test_ah_has_three_per_particle_params(self, simple_topology, chain_meta):
        """Test AH force has sigma, lambda, id parameters."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        ah = forces[0]
        assert ah.getNumPerParticleParameters() == 3

    def test_yukawa_has_one_per_particle_param(self, simple_topology, chain_meta):
        """Test Yukawa force has charge parameter."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        yu = forces[1]
        assert yu.getNumPerParticleParameters() == 1

    def test_correct_number_of_particles(self, simple_topology, chain_meta):
        """Test forces have correct number of particles."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        assert forces[0].getNumParticles() == 3
        assert forces[1].getNumParticles() == 3

    def test_uses_periodic_boundary_conditions(self, simple_topology, chain_meta):
        """Test forces use periodic boundary conditions."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        for force in forces:
            method = force.getNonbondedMethod()
            assert method == mm.CustomNonbondedForce.CutoffPeriodic

    def test_ah_cutoff_is_2nm(self, simple_topology, chain_meta):
        """Test AH force cutoff is 2 nm."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        cutoff = forces[0].getCutoffDistance().value_in_unit(unit.nanometer)
        assert cutoff == pytest.approx(2.0, abs=1e-6)

    def test_yukawa_cutoff_is_4nm(self, simple_topology, chain_meta):
        """Test Yukawa force cutoff is 4 nm."""
        ff = CalvadosFF(version=2)
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        cutoff = forces[1].getCutoffDistance().value_in_unit(unit.nanometer)
        assert cutoff == pytest.approx(4.0, abs=1e-6)


# =============================================================================
# Harmonic Bonds Tests
# =============================================================================

class TestBuildHarmonicBonds:
    """Tests for build_harmonic_bonds method."""

    @pytest.fixture
    def topology_with_different_residues(self) -> app.Topology:
        """Create topology with different residue types."""
        top = app.Topology()
        chain = top.addChain()
        
        residues = ["ALA", "GLY", "ALA"]
        atoms = []
        for res_name in residues:
            res = top.addResidue(res_name, chain)
            atom = top.addAtom("CA", app.element.carbon, res)
            atoms.append(atom)
        
        # Add bonds
        top.addBond(atoms[0], atoms[1])
        top.addBond(atoms[1], atoms[2])
        
        return top

    def test_uses_per_residue_r0(self, topology_with_different_residues):
        """Test that per-residue r0 values are used."""
        ff = CalvadosFF(version=2)
        force = ff.build_harmonic_bonds(topology_with_different_residues)
        
        # Each bond should have r0 from the first residue's parameters
        a1, a2, length1, k1 = force.getBondParameters(0)
        a3, a4, length2, k2 = force.getBondParameters(1)
        
        # ALA bond length should be used (both bonds start with ALA)
        expected_r0 = ff._params["A"]["r0"]
        assert length1.value_in_unit(unit.nanometer) == pytest.approx(expected_r0, abs=1e-6)
        assert length2.value_in_unit(unit.nanometer) == pytest.approx(expected_r0, abs=1e-6)

    def test_default_k_is_8368(self, topology_with_different_residues):
        """Test default spring constant is 8368 kJ/mol/nm^2."""
        ff = CalvadosFF(version=2)
        force = ff.build_harmonic_bonds(topology_with_different_residues)
        
        a1, a2, length, k = force.getBondParameters(0)
        expected_k = 8368.0
        actual_k = k.value_in_unit(unit.kilojoule_per_mole / unit.nanometer ** 2)
        assert actual_k == pytest.approx(expected_k, abs=1e-6)


# =============================================================================
# Constants Tests
# =============================================================================

class TestCalvadosConstants:
    """Tests for class constants."""

    def test_eps_lj_value(self):
        """Test epsilon LJ constant."""
        assert CalvadosFF.EPS_LJ == pytest.approx(0.8368, abs=1e-6)

    def test_rc_lj_value(self):
        """Test LJ cutoff constant."""
        assert CalvadosFF.RC_LJ == pytest.approx(2.0, abs=1e-6)

    def test_rc_yu_value(self):
        """Test Yukawa cutoff constant."""
        assert CalvadosFF.RC_YU == pytest.approx(4.0, abs=1e-6)

    def test_k_bond_value(self):
        """Test bond spring constant."""
        assert CalvadosFF.K_BOND == pytest.approx(8368.0, abs=1e-6)
