"""
Unit tests for core/topology.py

Tests cover:
- OpenMM Topology construction
- Mass table correctness
- Folded domain range calculation
- Box vector construction
"""

from __future__ import annotations

import numpy as np

from CondenSimAdapter.core.topology import (
    ONE_TO_THREE,
    RESIDUE_MASS,
    THREE_TO_ONE,
    _box_vectors,
    build_topology,
    get_folded_atom_ranges,
    get_masses,
)

# =============================================================================
# Residue Constants Tests
# =============================================================================


class TestResidueConstants:
    """Tests for residue conversion tables."""

    def test_three_to_one_coverage(self):
        """Test that all 20 standard amino acids are covered."""
        assert len(THREE_TO_ONE) == 20
        expected = set("ARNDCEQGHILKMFPSTWYV")
        assert set(THREE_TO_ONE.values()) == expected

    def test_one_to_three_inverse(self):
        """Test that ONE_TO_THREE is correct inverse of THREE_TO_ONE."""
        for three, one in THREE_TO_ONE.items():
            assert ONE_TO_THREE[one] == three

    def test_residue_mass_coverage(self):
        """Test that all 20 amino acids have mass entries."""
        for aa in ONE_TO_THREE.keys():
            assert aa in RESIDUE_MASS
            assert RESIDUE_MASS[aa] > 0

    def test_residue_mass_reasonable(self):
        """Test that masses are in reasonable range for amino acids."""
        for aa, mass in RESIDUE_MASS.items():
            assert 50 < mass < 200, f"Mass for {aa} seems unreasonable: {mass}"


# =============================================================================
# Box Vectors Tests
# =============================================================================


class TestBoxVectors:
    """Tests for _box_vectors helper."""

    def test_cubic_box(self):
        """Test cubic box vector construction."""
        box = [10.0, 10.0, 10.0]
        vectors = _box_vectors(box)
        # Should be a Quantity with 3 Vec3 objects
        assert len(vectors.value_in_unit(vectors.unit)) == 3
        v0, v1, v2 = vectors.value_in_unit(vectors.unit)
        np.testing.assert_allclose([v0.x, v0.y, v0.z], [10.0, 0, 0])
        np.testing.assert_allclose([v1.x, v1.y, v1.z], [0, 10.0, 0])
        np.testing.assert_allclose([v2.x, v2.y, v2.z], [0, 0, 10.0])

    def test_rectangular_box(self):
        """Test rectangular box vector construction."""
        box = [20.0, 15.0, 10.0]
        vectors = _box_vectors(box)
        v0, v1, v2 = vectors.value_in_unit(vectors.unit)
        assert v0.x == 20.0
        assert v1.y == 15.0
        assert v2.z == 10.0


# =============================================================================
# Build Topology Tests
# =============================================================================


class TestBuildTopology:
    """Tests for build_topology function."""

    def test_single_chain(self):
        """Test building topology for single chain."""
        chain_meta = [
            {
                "name": "FUS",
                "start": 0,
                "end": 5,
                "sequence": "AGSVL",
                "folded_domains": [],
            }
        ]
        positions = np.zeros((5, 3), dtype=np.float64)
        box = [20.0, 20.0, 20.0]

        top, pos_qty = build_topology(chain_meta, positions, box)

        # Check topology structure
        assert top.getNumChains() == 1
        assert top.getNumResidues() == 5
        assert top.getNumAtoms() == 5
        assert top.getNumBonds() == 4  # 5 residues -> 4 bonds

        # Check positions are wrapped in Quantity
        assert hasattr(pos_qty, "unit")
        assert pos_qty.shape == (5, 3)

    def test_multiple_chains(self):
        """Test building topology for multiple chains."""
        chain_meta = [
            {"name": "A", "start": 0, "end": 3, "sequence": "AAA", "folded_domains": []},
            {"name": "B", "start": 3, "end": 7, "sequence": "BBBB", "folded_domains": []},
        ]
        positions = np.zeros((7, 3), dtype=np.float64)
        box = [20.0, 20.0, 20.0]

        top, _ = build_topology(chain_meta, positions, box)

        assert top.getNumChains() == 2
        assert top.getNumResidues() == 7
        assert top.getNumAtoms() == 7
        # Bonds: chain A has 2, chain B has 3, total 5
        assert top.getNumBonds() == 5

    def test_residue_names(self):
        """Test that residue names are set correctly from sequence."""
        chain_meta = [
            {
                "name": "test",
                "start": 0,
                "end": 3,
                "sequence": "AGS",  # Ala-Gly-Ser
                "folded_domains": [],
            }
        ]
        positions = np.zeros((3, 3))
        box = [10.0, 10.0, 10.0]

        top, _ = build_topology(chain_meta, positions, box)

        residues = list(top.residues())
        assert residues[0].name == "ALA"
        assert residues[1].name == "GLY"
        assert residues[2].name == "SER"

    def test_atom_names(self):
        """Test that all atoms are named CA."""
        chain_meta = [
            {
                "name": "test",
                "start": 0,
                "end": 5,
                "sequence": "AAAAA",
                "folded_domains": [],
            }
        ]
        positions = np.zeros((5, 3))
        box = [10.0, 10.0, 10.0]

        top, _ = build_topology(chain_meta, positions, box)

        for atom in top.atoms():
            assert atom.name == "CA"

    def test_bonds_within_chains_only(self):
        """Test that bonds don't connect different chains."""
        chain_meta = [
            {"name": "A", "start": 0, "end": 2, "sequence": "AA", "folded_domains": []},
            {"name": "B", "start": 2, "end": 4, "sequence": "BB", "folded_domains": []},
        ]
        positions = np.zeros((4, 3))
        box = [10.0, 10.0, 10.0]

        top, _ = build_topology(chain_meta, positions, box)

        # Each chain should have 1 bond
        assert top.getNumBonds() == 2

        # Verify bonds are within chains
        list(top.atoms())
        for bond in top.bonds():
            a1, a2 = bond[0], bond[1]
            # Atoms should be consecutive in the same chain
            assert abs(a1.index - a2.index) == 1

    def test_periodic_box_set(self):
        """Test that periodic box vectors are set correctly."""
        chain_meta = [
            {
                "name": "test",
                "start": 0,
                "end": 3,
                "sequence": "AAA",
                "folded_domains": [],
            }
        ]
        positions = np.zeros((3, 3))
        box = [25.0, 25.0, 25.0]

        top, _ = build_topology(chain_meta, positions, box)

        vectors = top.getPeriodicBoxVectors()
        assert vectors is not None
        v0, v1, v2 = vectors.value_in_unit(vectors.unit)
        assert v0.x == 25.0
        assert v1.y == 25.0
        assert v2.z == 25.0


# =============================================================================
# Get Masses Tests
# =============================================================================


class TestGetMasses:
    """Tests for get_masses function."""

    def test_single_chain_masses(self):
        """Test mass calculation for single chain."""
        chain_meta = [
            {
                "name": "test",
                "sequence": "AG",
            }
        ]
        masses = get_masses(chain_meta)

        assert len(masses) == 2
        assert masses[0] == RESIDUE_MASS["A"]  # Alanine
        assert masses[1] == RESIDUE_MASS["G"]  # Glycine

    def test_multiple_chains_masses(self):
        """Test mass calculation for multiple chains."""
        chain_meta = [
            {"name": "A", "sequence": "A"},
            {"name": "B", "sequence": "GS"},
        ]
        masses = get_masses(chain_meta)

        assert len(masses) == 3
        assert masses[0] == RESIDUE_MASS["A"]
        assert masses[1] == RESIDUE_MASS["G"]
        assert masses[2] == RESIDUE_MASS["S"]

    def test_unknown_residue_fallback(self):
        """Test that unknown residues get default mass."""
        chain_meta = [
            {
                "name": "test",
                "sequence": "AXA",  # X is unknown
            }
        ]
        masses = get_masses(chain_meta)

        assert len(masses) == 3
        assert masses[0] == RESIDUE_MASS["A"]
        assert masses[1] == 57.05  # Default (Glycine) mass
        assert masses[2] == RESIDUE_MASS["A"]

    def test_return_type(self):
        """Test that return is numpy array of float64."""
        chain_meta = [{"name": "test", "sequence": "AAAA"}]
        masses = get_masses(chain_meta)

        assert isinstance(masses, np.ndarray)
        assert masses.dtype == np.float64


# =============================================================================
# Get Folded Atom Ranges Tests
# =============================================================================


class TestGetFoldedAtomRanges:
    """Tests for get_folded_atom_ranges function."""

    def test_single_domain(self):
        """Test single folded domain."""
        chain_meta = [
            {
                "name": "MDP1",
                "start": 0,
                "end": 100,
                "sequence": "A" * 100,
                "folded_domains": [(10, 30)],  # 1-based inclusive
            }
        ]
        ranges = get_folded_atom_ranges(chain_meta)

        assert len(ranges) == 1
        # 10-30 (1-based) -> 9-30 (0-based, exclusive end)
        assert ranges[0] == (9, 30)

    def test_multiple_domains_single_chain(self):
        """Test multiple folded domains in one chain."""
        chain_meta = [
            {
                "name": "MDP1",
                "start": 0,
                "end": 200,
                "sequence": "A" * 200,
                "folded_domains": [(1, 50), (100, 150)],
            }
        ]
        ranges = get_folded_atom_ranges(chain_meta)

        assert len(ranges) == 2
        assert ranges[0] == (0, 50)  # 1-50 -> 0-50
        assert ranges[1] == (99, 150)  # 100-150 -> 99-150

    def test_multiple_chains_with_domains(self):
        """Test folded domains across multiple chains."""
        chain_meta = [
            {
                "name": "A",
                "start": 0,
                "end": 50,
                "sequence": "A" * 50,
                "folded_domains": [(1, 20)],
            },
            {
                "name": "B",
                "start": 50,
                "end": 100,
                "sequence": "A" * 50,
                "folded_domains": [(10, 30)],
            },
        ]
        ranges = get_folded_atom_ranges(chain_meta)

        assert len(ranges) == 2
        # Chain A: domain at 1-20, chain starts at 0 -> absolute 0-20
        assert ranges[0] == (0, 20)
        # Chain B: domain at 10-30, chain starts at 50 -> absolute 59-80
        assert ranges[1] == (59, 80)

    def test_no_folded_domains(self):
        """Test chain with no folded domains."""
        chain_meta = [
            {
                "name": "IDP",
                "start": 0,
                "end": 50,
                "sequence": "A" * 50,
                "folded_domains": [],
            }
        ]
        ranges = get_folded_atom_ranges(chain_meta)

        assert ranges == []

    def test_mixed_chains(self):
        """Test system with both folded and unfolded chains."""
        chain_meta = [
            {
                "name": "MDP",
                "start": 0,
                "end": 100,
                "sequence": "A" * 100,
                "folded_domains": [(1, 50)],
            },
            {
                "name": "IDP",
                "start": 100,
                "end": 150,
                "sequence": "A" * 50,
                "folded_domains": [],
            },
        ]
        ranges = get_folded_atom_ranges(chain_meta)

        assert len(ranges) == 1
        assert ranges[0] == (0, 50)
