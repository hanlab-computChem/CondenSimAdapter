"""
Unit tests for minimize/topology_builder.py

Tests cover:
- PDB2GMX input building
- Topology modification helpers
- Utility functions
"""

from __future__ import annotations

import pytest

from CondenSimAdapter.minimize.topology_builder import (
    _build_pdb2gmx_input,
    count_chains_in_pdb,
    modify_topology_molecule_name,
    verify_files_valid,
)

# =============================================================================
# PDB2GMX Input Building Tests
# =============================================================================


class TestBuildPDB2GMXInput:
    """Tests for _build_pdb2gmx_input function."""

    def test_no_options_returns_none(self):
        """Test with no options returns None."""
        result = _build_pdb2gmx_input(disable_disulfide=False, his_type=None, his_repeat_count=30)
        assert result is None

    def test_disulfide_only(self):
        """Test disulfide disabling only."""
        result = _build_pdb2gmx_input(disable_disulfide=True, his_type=None, his_repeat_count=30)
        assert result is not None
        assert "n\n" in result

    def test_his_type_only(self):
        """Test histidine type only."""
        result = _build_pdb2gmx_input(
            disable_disulfide=False,
            his_type=0,  # HID
            his_repeat_count=5,
        )
        assert result is not None
        assert "0\n" in result
        # Should have 5 repetitions
        assert result.count("0\n") == 5

    def test_both_options(self):
        """Test both disulfide and histidine options."""
        result = _build_pdb2gmx_input(
            disable_disulfide=True,
            his_type=1,  # HIE
            his_repeat_count=3,
        )
        assert result is not None
        assert "n\n" in result
        assert "1\n" in result

    def test_invalid_his_type_raises(self):
        """Test invalid histidine type raises ValueError."""
        with pytest.raises(ValueError, match="his_type must be 0 or 1"):
            _build_pdb2gmx_input(
                disable_disulfide=False,
                his_type=2,  # Invalid
                his_repeat_count=30,
            )


# =============================================================================
# Topology Modification Tests
# =============================================================================


class TestModifyTopologyMoleculeName:
    """Tests for modify_topology_molecule_name function."""

    def test_modifies_molecule_name(self, tmp_path):
        """Test molecule name is modified."""
        topology_content = """\
; This is a topology
#include "forcefield.itp"

[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
1 ALA CA 1 ALA CA 1 0.0 12.01
"""
        top_path = tmp_path / "topol.top"
        top_path.write_text(topology_content)

        modify_topology_molecule_name(top_path, "NewName")

        result = top_path.read_text()
        assert "NewName" in result
        assert "Protein" not in result

    def test_preserves_other_content(self, tmp_path):
        """Test other topology content is preserved."""
        topology_content = """\
[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
1 ALA CA 1 ALA CA 1 0.0 12.01

[ bonds ]
1 2 1 0.1 1000
"""
        top_path = tmp_path / "topol.top"
        top_path.write_text(topology_content)

        modify_topology_molecule_name(top_path, "Test")

        result = top_path.read_text()
        assert "[ atoms ]" in result
        assert "[ bonds ]" in result
        assert "Test" in result


# =============================================================================
# Chain Counting Tests
# =============================================================================


class TestCountChainsInPDB:
    """Tests for count_chains_in_pdb function."""

    def test_single_chain(self, tmp_path):
        """Test single chain PDB."""
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00  0.00           C
END
"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)

        count = count_chains_in_pdb(pdb_path)
        assert count == 1

    def test_multiple_chains(self, tmp_path):
        """Test multiple chain PDB."""
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA B   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA  ALA C   1       2.000   0.000   0.000  1.00  0.00           C
END
"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)

        count = count_chains_in_pdb(pdb_path)
        assert count == 3

    def test_hetatm_included(self, tmp_path):
        """Test HETATM records are included."""
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O   HOH B   1       1.000   0.000   0.000  1.00  0.00           O
END
"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)

        count = count_chains_in_pdb(pdb_path)
        assert count == 2

    def test_empty_chain_id(self, tmp_path):
        """Test empty chain ID is handled."""
        pdb_content = """\
ATOM      1  CA  ALA     1       0.000   0.000   0.000  1.00  0.00           C
END
"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)

        count = count_chains_in_pdb(pdb_path)
        # Empty chain ID should not be counted
        assert count == 0


# =============================================================================
# File Verification Tests
# =============================================================================


class TestVerifyFilesValid:
    """Tests for verify_files_valid function."""

    def test_missing_topology_file(self, tmp_path):
        """Test missing topology file returns error."""
        top_path = tmp_path / "missing.top"
        gro_path = tmp_path / "test.gro"
        gro_path.write_text("test")

        success, message = verify_files_valid(str(top_path), str(gro_path))
        assert success is False
        assert "not found" in message.lower()

    def test_missing_gro_file(self, tmp_path):
        """Test missing GRO file returns error."""
        top_path = tmp_path / "test.top"
        gro_path = tmp_path / "missing.gro"
        top_path.write_text("test")

        success, message = verify_files_valid(str(top_path), str(gro_path))
        assert success is False
        assert "not found" in message.lower()

    def test_topology_missing_atoms_section(self, tmp_path):
        """Test topology without [ atoms ] section fails."""
        top_path = tmp_path / "test.top"
        gro_path = tmp_path / "test.gro"

        top_path.write_text("""\
[ moleculetype ]
Test 3

[ molecules ]
Test 1
""")
        # Create minimal valid GRO
        gro_path.write_text("""\
Test
2
    1ALA    CA    1   0.000   0.000   0.000
    1ALA     C    2   0.100   0.000   0.000
   0.10000   0.10000   0.10000
""")

        success, message = verify_files_valid(str(top_path), str(gro_path))
        assert success is False
        assert "missing [ atoms ]" in message

    def test_topology_missing_molecules_section(self, tmp_path):
        """Test topology without [ molecules ] section fails."""
        top_path = tmp_path / "test.top"
        gro_path = tmp_path / "test.gro"

        top_path.write_text("""\
[ moleculetype ]
Test 3

[ atoms ]
1 CA 1 ALA CA 1 0.0 12.01
""")
        # Create minimal valid GRO
        gro_path.write_text("""\
Test
1
    1ALA    CA    1   0.000   0.000   0.000
   0.10000   0.10000   0.10000
""")

        success, message = verify_files_valid(str(top_path), str(gro_path))
        assert success is False
        assert "missing [ molecules ]" in message
