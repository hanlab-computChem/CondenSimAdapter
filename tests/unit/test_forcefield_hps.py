"""
Unit tests for core/forcefield/hps.py

Tests cover:
- HPS-Urry parameter loading
- 20x20 sigma/lambda tables
- Nonbonded force creation (AH + Debye-Huckel)
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("openmm")

import openmm.app as app
import openmm.unit as unit

from CondenSimAdapter.core.forcefield.hps import _AA_ORDER, _N_AA, HPSFF


class TestHPSInitialization:
    """Tests for HPSFF initialization."""

    def test_parameters_loaded(self):
        """Test that parameter tables are loaded."""
        ff = HPSFF()
        assert ff._sigma_table.shape == (_N_AA, _N_AA)
        assert ff._lambda_table.shape == (_N_AA, _N_AA)
        assert len(ff._charges) > 0

    def test_tables_symmetric(self):
        """Test that parameter tables are symmetric."""
        ff = HPSFF()
        np.testing.assert_allclose(ff._sigma_table, ff._sigma_table.T)
        np.testing.assert_allclose(ff._lambda_table, ff._lambda_table.T)

    def test_charges_expected_values(self):
        """Test that expected charged residues have correct charges."""
        ff = HPSFF()
        assert ff._charges.get("ARG") == 1.0
        assert ff._charges.get("LYS") == 1.0
        assert ff._charges.get("ASP") == -1.0
        assert ff._charges.get("GLU") == -1.0
        assert ff._charges.get("HIS") == 0.5


class TestHPSSigmaLambdaTables:
    """Tests for sigma and lambda parameter tables."""

    def test_sigma_positive(self):
        """Test all sigma values are positive."""
        ff = HPSFF()
        assert np.all(ff._sigma_table > 0)

    def test_lambda_within_range(self):
        """Test lambda values are in reasonable range."""
        ff = HPSFF()
        # Lambda should be between -delta and 1-delta after correction
        assert np.all(ff._lambda_table >= -0.1)
        assert np.all(ff._lambda_table <= 1.0)


class TestHPSNonbondedForces:
    """Tests for create_nonbonded_forces method."""

    @pytest.fixture
    def simple_topology(self) -> app.Topology:
        """Create a simple 3-residue topology."""
        top = app.Topology()
        chain = top.addChain()

        for res_name in ["ALA", "GLY", "ASP"]:
            res = top.addResidue(res_name, chain)
            top.addAtom("CA", app.element.carbon, res)

        atoms = list(top.atoms())
        top.addBond(atoms[0], atoms[1])
        top.addBond(atoms[1], atoms[2])

        return top

    @pytest.fixture
    def chain_meta(self):
        return [{"name": "test", "sequence": "AGD"}]

    def test_returns_two_forces(self, simple_topology, chain_meta):
        """Test that AH and DH forces are returned."""
        ff = HPSFF()
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        assert len(forces) == 2

    def test_ah_uses_tabulated_functions(self, simple_topology, chain_meta):
        """Test AH force uses Discrete2DFunction for tables."""
        ff = HPSFF()
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        ah = forces[0]
        # Check that tabulated functions exist
        assert ah.getNumTabulatedFunctions() == 2

    def test_dynamic_cutoff(self, simple_topology, chain_meta):
        """Test AH cutoff is based on max sigma."""
        ff = HPSFF()
        forces = ff.create_nonbonded_forces(
            simple_topology, chain_meta, temperature=300.0, ionic=0.15
        )
        ah = forces[0]
        max_sigma = float(np.max(ff._sigma_table))
        expected_cutoff = 4.0 * max_sigma
        actual_cutoff = ah.getCutoffDistance().value_in_unit(unit.nanometer)
        assert actual_cutoff == pytest.approx(expected_cutoff, abs=1e-6)


class TestHPSConstants:
    """Tests for HPS class constants."""

    def test_aa_order_has_20(self):
        """Test AA order list has 20 amino acids."""
        assert len(_AA_ORDER) == 20

    def test_n_aa_is_20(self):
        """Test N_AA constant."""
        assert _N_AA == 20

    def test_eps_ah_value(self):
        """Test epsilon AH constant."""
        assert HPSFF.EPS_AH == pytest.approx(0.8368, abs=1e-6)

    def test_delta_value(self):
        """Test delta correction constant."""
        assert HPSFF.DELTA == pytest.approx(0.08, abs=1e-6)

    def test_k_bond_value(self):
        """Test bond spring constant."""
        assert HPSFF.K_BOND == pytest.approx(8368.0, abs=1e-6)
