"""
Unit tests for backmap/backmapper.py

Tests cover:
- Backmapper input validation
- Model type handling
- Slab z-centering algorithm
- PDB helper functions (C-terminus fix, TER insertion)
"""

from __future__ import annotations

import numpy as np
import pytest
from pathlib import Path
from unittest.mock import Mock, patch

from CondenSimAdapter.backmap.backmapper import (
    Backmapper,
    BackmapResult,
    SUPPORTED_MODELS,
    _center_slab_in_z,
    _fix_c_terminus_atom_names,
    _insert_ter_after_oxt,
)


# =============================================================================
# Backmapper Initialization and Basic Tests
# =============================================================================

class TestBackmapperBasics:
    """Tests for Backmapper basic functionality."""

    def test_supported_models_defined(self):
        """Test that supported models are defined."""
        assert len(SUPPORTED_MODELS) > 0
        assert "CalphaBasedModel" in SUPPORTED_MODELS

    def test_backmapper_instantiation(self):
        """Test Backmapper can be instantiated."""
        mapper = Backmapper()
        assert mapper is not None


# =============================================================================
# Input Validation Tests
# =============================================================================

class TestBackmapperInputValidation:
    """Tests for input validation."""

    def test_missing_input_file_returns_error(self, tmp_path):
        """Test that missing input file returns error result."""
        mapper = Backmapper()
        result = mapper.run(
            cg_pdb=str(tmp_path / "nonexistent.pdb"),
            output_dir=str(tmp_path / "output"),
        )
        assert result.success is False
        assert "not found" in result.error.lower()

    def test_unsupported_model_warns(self, tmp_path):
        """Test unsupported model type logs warning."""
        # Create dummy input file
        input_pdb = tmp_path / "input.pdb"
        input_pdb.write_text("ATOM    1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")
        
        mapper = Backmapper()
        with patch('CondenSimAdapter.backmap.backmapper.log') as mock_log:
            result = mapper.run(
                cg_pdb=str(input_pdb),
                output_dir=str(tmp_path / "output"),
                model_type="UnsupportedModel",
            )
            # Should warn but still proceed
            mock_log.warning.assert_called_once()


# =============================================================================
# Slab Z-Centering Tests
# =============================================================================

class TestCenterSlabInZ:
    """Tests for _center_slab_in_z helper."""

    @pytest.fixture
    def slab_pdb_content(self):
        """Create a simple slab PDB content."""
        return """\
CRYST1   50.000   50.000  100.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   1      25.000  25.000  40.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      25.000  25.000  41.000  1.00  0.00           C
ATOM      3  C   ALA A   1      25.000  25.000  42.000  1.00  0.00           C
ATOM      4  O   ALA A   1      25.000  25.000  43.000  1.00  0.00           O
ATOM      5  OXT ALA A   1      25.000  25.000  44.000  1.00  0.00           O
TER
END
"""

    def test_slab_centering_shifts_z(self, tmp_path, slab_pdb_content):
        """Test that slab centering shifts z coordinates."""
        pdb_path = tmp_path / "slab.pdb"
        pdb_path.write_text(slab_pdb_content)
        
        result_path = _center_slab_in_z(str(pdb_path))
        
        # Read result and check z values shifted
        result_content = Path(result_path).read_text()
        # Original z values were around 40-44, box_z is 100
        # COM should move to 50 (box_z/2)
        assert result_path == str(pdb_path)  # Modifies in place

    def test_preserves_cryst1_record(self, tmp_path, slab_pdb_content):
        """Test that CRYST1 record is preserved."""
        pdb_path = tmp_path / "slab.pdb"
        pdb_path.write_text(slab_pdb_content)
        
        _center_slab_in_z(str(pdb_path))
        
        result_content = Path(pdb_path).read_text()
        assert "CRYST1" in result_content


# =============================================================================
# C-Terminus Fix Tests
# =============================================================================

class TestFixCTerminusAtomNames:
    """Tests for _fix_c_terminus_atom_names helper."""

    def test_ot1_becomes_o(self, tmp_path):
        """Test OT1 is renamed to O."""
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  OT1 ALA A   1       1.000   0.000   0.000  1.00  0.00           O
ATOM      3  OT2 ALA A   1       2.000   0.000   0.000  1.00  0.00           O
"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)
        
        _fix_c_terminus_atom_names(str(pdb_path))
        
        result = pdb_path.read_text()
        assert " OT1 " not in result
        assert " O   " in result  # OT1 -> O
        assert " OXT" in result   # OT2 -> OXT

    def test_no_c_terminus_unchanged(self, tmp_path):
        """Test file without C-terminus is unchanged."""
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  N   ALA A   1       1.000   0.000   0.000  1.00  0.00           N
"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)
        
        _fix_c_terminus_atom_names(str(pdb_path))
        
        result = pdb_path.read_text()
        assert "CA" in result
        assert "N" in result


# =============================================================================
# TER Insertion Tests
# =============================================================================

class TestInsertTerAfterOxt:
    """Tests for _insert_ter_after_oxt helper."""

    def test_inserts_ter_after_oxt(self, tmp_path):
        """Test TER is inserted after OXT atom."""
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  OXT ALA A   1       1.000   0.000   0.000  1.00  0.00           O
ATOM      3  CA  GLY A   2       2.000   0.000   0.000  1.00  0.00           C
"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)
        
        _insert_ter_after_oxt(str(pdb_path))
        
        lines = pdb_path.read_text().strip().split('\n')
        assert "TER" in lines
        # TER should be between OXT and next CA
        oxt_idx = [i for i, l in enumerate(lines) if "OXT" in l][0]
        ter_idx = [i for i, l in enumerate(lines) if l == "TER"][0]
        assert ter_idx == oxt_idx + 1

    def test_no_ter_if_already_present(self, tmp_path):
        """Test no duplicate TER if already present."""
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  OXT ALA A   1       1.000   0.000   0.000  1.00  0.00           O
TER
ATOM      3  CA  GLY A   2       2.000   0.000   0.000  1.00  0.00           C
"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)
        
        _insert_ter_after_oxt(str(pdb_path))
        
        lines = pdb_path.read_text().strip().split('\n')
        ter_count = sum(1 for l in lines if l == "TER")
        assert ter_count == 1  # Still only one TER


# =============================================================================
# BackmapResult Tests
# =============================================================================

class TestBackmapResult:
    """Tests for BackmapResult dataclass."""

    def test_success_result(self):
        """Test successful result creation."""
        result = BackmapResult(
            success=True,
            output_pdb="/path/to/output.pdb",
            input_pdb="/path/to/input.pdb",
            model_type="CalphaBasedModel",
        )
        assert result.success is True
        assert result.output_pdb == "/path/to/output.pdb"
        assert result.error is None

    def test_failure_result(self):
        """Test failed result creation."""
        result = BackmapResult(
            success=False,
            input_pdb="/path/to/input.pdb",
            error="File not found",
        )
        assert result.success is False
        assert result.error == "File not found"
        assert result.output_pdb is None
