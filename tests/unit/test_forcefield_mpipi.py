"""
Unit tests for core/forcefield/mpipi.py

Tests cover:
- Mpipi parameter loading from recharged_params.txt
- Wang-Frenkel + Yukawa force creation
- Globular domain handling
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("openmm")

import openmm as mm
import openmm.app as app
import openmm.unit as unit

from CondenSimAdapter.core.forcefield.mpipi import (
    MpipiFF,
    _MPIPI_ORDER,
    _N_MPIPI,
    _AA3_TO_MPIPI,
    _MPIPI_MASS,
)


class TestMpipiInitialization:
    """Tests for MpipiFF initialization."""

    def test_parameters_loaded(self):
        """Test that WF and Yukawa parameters are loaded."""
        ff = MpipiFF()
        # WF params: N_MPIPI * N_MPIPI * 3 values
        expected_wf_size = _N_MPIPI * _N_MPIPI * 3
        assert len(ff._wf_params) == expected_wf_size
        # Yukawa A: N_MPIPI * N_MPIPI values
        expected_yukawa_size = _N_MPIPI * _N_MPIPI
        assert len(ff._yukawa_A) == expected_yukawa_size


class TestMpipiMappings:
    """Tests for residue mappings."""

    def test_mpi_order_has_20(self):
        """Test MPI order has 20 entries (20 amino acids, rU excluded in recharged version)."""
        assert len(_MPIPI_ORDER) == 20

    def test_aa3_to_mpipi_has_20(self):
        """Test 3-letter to Mpipi index mapping has 20 entries."""
        assert len(_AA3_TO_MPIPI) == 20

    def test_mpipi_mass_has_20(self):
        """Test mass table has 20 entries."""
        assert len(_MPIPI_MASS) == 20

    def test_common_residues_mapped(self):
        """Test common residues are in mapping."""
        for res in ["ALA", "GLY", "ARG", "LYS", "ASP", "GLU"]:
            assert res in _AA3_TO_MPIPI
            assert res in _MPIPI_MASS


class TestMpipiNonbondedForces:
    """Tests for create_nonbonded_forces method."""

    @pytest.fixture
    def simple_topology(self) -> app.Topology:
        """Create a simple topology."""
        top = app.Topology()
        chain = top.addChain()
        
        for res_name in ["ALA", "GLY", "ARG"]:
            res = top.addResidue(res_name, chain)
            top.addAtom("CA", app.element.carbon, res)
        
        atoms = list(top.atoms())
        top.addBond(atoms[0], atoms[1])
        top.addBond(atoms[1], atoms[2])
        
        return top

    @pytest.fixture
    def chain_meta(self):
        return [{"name": "test", "sequence": "AGR", "start": 0, "end": 3, "folded_domains": []}]

    def test_returns_two_forces(self, simple_topology, chain_meta):
        """Test that WF and Yukawa forces are returned."""
        ff = MpipiFF()
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        assert len(forces) == 2

    def test_wf_uses_3d_table(self, simple_topology, chain_meta):
        """Test WF force uses Discrete3DFunction."""
        ff = MpipiFF()
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        wf = forces[0]
        assert wf.getNumTabulatedFunctions() == 1

    def test_yukawa_uses_2d_table(self, simple_topology, chain_meta):
        """Test Yukawa force uses Discrete2DFunction."""
        ff = MpipiFF()
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        yu = forces[1]
        assert yu.getNumTabulatedFunctions() == 1

    def test_wf_has_globular_parameter(self, simple_topology, chain_meta):
        """Test WF force has globular parameter for domain scaling."""
        ff = MpipiFF()
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        wf = forces[0]
        # Should have 2 per-particle params: index and globular
        assert wf.getNumPerParticleParameters() == 2


class TestMpipiFoldedSet:
    """Tests for _get_folded_set method."""

    def test_no_domains_returns_empty(self):
        """Test empty folded domains returns empty set."""
        chain_meta = [{
            "name": "test",
            "start": 0,
            "end": 10,
            "sequence": "AAAAAAAAAA",
            "folded_domains": [],
        }]
        folded = MpipiFF._get_folded_set(chain_meta)
        assert len(folded) == 0

    def test_single_domain(self):
        """Test single folded domain returns correct indices."""
        chain_meta = [{
            "name": "test",
            "start": 0,
            "end": 10,
            "sequence": "AAAAAAAAAA",
            "folded_domains": [(2, 5)],  # 1-based: atoms 1,2,3,4
        }]
        folded = MpipiFF._get_folded_set(chain_meta)
        assert folded == {1, 2, 3, 4}

    def test_multiple_domains(self):
        """Test multiple folded domains."""
        chain_meta = [
            {
                "name": "A",
                "start": 0,
                "end": 10,
                "sequence": "AAAAAAAAAA",
                "folded_domains": [(1, 3)],
            },
            {
                "name": "B",
                "start": 10,
                "end": 20,
                "sequence": "AAAAAAAAAA",
                "folded_domains": [(5, 8)],  # 1-based in chain B
            },
        ]
        folded = MpipiFF._get_folded_set(chain_meta)
        # Chain A: domain at 1-3 -> indices 0,1,2
        # Chain B: domain at 5-8, chain starts at 10 -> indices 14,15,16,17
        assert folded == {0, 1, 2, 14, 15, 16, 17}


class TestMpipiConstants:
    """Tests for Mpipi class constants."""

    def test_k_bond_value(self):
        """Test bond spring constant."""
        from CondenSimAdapter.core.forcefield.mpipi import _K_BOND
        assert _K_BOND == pytest.approx(8031.0, abs=1e-6)

    def test_d_idr_value(self):
        """Test IDR bond length."""
        from CondenSimAdapter.core.forcefield.mpipi import _D_IDR
        assert _D_IDR == pytest.approx(0.381, abs=1e-6)
