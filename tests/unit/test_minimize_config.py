"""
Unit tests for minimize/config_loader.py and minimizer.py config

Tests cover:
- MinimizeConfig defaults and validation
- Config loading from YAML
"""

from __future__ import annotations

import pytest
from pathlib import Path
import yaml

from CondenSimAdapter.minimize.minimizer import MinimizeConfig, MinimizeResult


# =============================================================================
# MinimizeConfig Tests
# =============================================================================

class TestMinimizeConfigDefaults:
    """Tests for MinimizeConfig default values."""

    def test_default_forcefield_type(self):
        """Test default force field type."""
        config = MinimizeConfig()
        assert config.forcefield_type == "1-a99SBdisp"

    def test_default_gb_model(self):
        """Test default GB model."""
        config = MinimizeConfig()
        assert config.gb_model == "GBn2"

    def test_default_softcore_lambdas(self):
        """Test default softcore lambda schedule."""
        config = MinimizeConfig()
        assert config.softcore_lambdas == [0.75, 0.85, 0.95]

    def test_default_platform(self):
        """Test default platform."""
        config = MinimizeConfig()
        assert config.platform == "CUDA"

    def test_default_gpu_id(self):
        """Test default GPU ID."""
        config = MinimizeConfig()
        assert config.gpu_id == 0

    def test_default_tolerance(self):
        """Test default minimization tolerance."""
        config = MinimizeConfig()
        assert config.tolerance == 100.0

    def test_default_max_iterations(self):
        """Test default max iterations."""
        config = MinimizeConfig()
        assert config.max_iterations == 5000

    def test_default_nonbonded_cutoff(self):
        """Test default nonbonded cutoff."""
        config = MinimizeConfig()
        assert config.nonbonded_cutoff == 2.0

    def test_default_solvate(self):
        """Test default solvation flag."""
        config = MinimizeConfig()
        assert config.solvate is False

    def test_default_ion_concentration(self):
        """Test default ion concentration."""
        config = MinimizeConfig()
        assert config.ion_concentration == 0.15

    def test_default_box_resize(self):
        """Test default box resize flag."""
        config = MinimizeConfig()
        assert config.box_resize is False


class TestMinimizeConfigCustomization:
    """Tests for MinimizeConfig customization."""

    def test_custom_forcefield(self):
        """Test custom force field type."""
        config = MinimizeConfig(forcefield_type="8-amber14sb")
        assert config.forcefield_type == "8-amber14sb"

    def test_custom_gb_model(self):
        """Test custom GB model."""
        config = MinimizeConfig(gb_model="OBC2")
        assert config.gb_model == "OBC2"

    def test_custom_softcore_schedule(self):
        """Test custom softcore lambda schedule."""
        config = MinimizeConfig(softcore_lambdas=[0.5, 0.8, 0.9, 1.0])
        assert config.softcore_lambdas == [0.5, 0.8, 0.9, 1.0]

    def test_custom_platform(self):
        """Test custom platform."""
        config = MinimizeConfig(platform="CPU")
        assert config.platform == "CPU"

    def test_custom_tolerance(self):
        """Test custom tolerance."""
        config = MinimizeConfig(tolerance=50.0)
        assert config.tolerance == 50.0

    def test_custom_max_iterations(self):
        """Test custom max iterations."""
        config = MinimizeConfig(max_iterations=10000)
        assert config.max_iterations == 10000

    def test_enable_solvation(self):
        """Test enabling explicit solvation."""
        config = MinimizeConfig(solvate=True, ion_concentration=0.2)
        assert config.solvate is True
        assert config.ion_concentration == 0.2

    def test_droplet_box_settings(self):
        """Test droplet box settings."""
        config = MinimizeConfig(
            droplet_box_type="dodecahedron",
            droplet_distance=3.0,
        )
        assert config.droplet_box_type == "dodecahedron"
        assert config.droplet_distance == 3.0

    def test_box_resize_settings(self):
        """Test box resize settings."""
        config = MinimizeConfig(
            box_resize=True,
            box_resize_dims=[20.0, 20.0, 20.0],
        )
        assert config.box_resize is True
        assert config.box_resize_dims == [20.0, 20.0, 20.0]


# =============================================================================
# MinimizeResult Tests
# =============================================================================

class TestMinimizeResult:
    """Tests for MinimizeResult dataclass."""

    def test_success_result(self):
        """Test successful result creation."""
        result = MinimizeResult(
            success=True,
            output_pdb="/path/to/minimized.pdb",
            input_pdb="/path/to/input.pdb",
            output_dir="/path/to/output",
        )
        assert result.success is True
        assert result.output_pdb == "/path/to/minimized.pdb"
        assert result.errors == []

    def test_failure_result(self):
        """Test failed result with errors."""
        result = MinimizeResult(
            success=False,
            output_pdb="",
            input_pdb="/path/to/input.pdb",
            output_dir="/path/to/output",
            errors=["Error 1", "Error 2"],
        )
        assert result.success is False
        assert len(result.errors) == 2
        assert "Error 1" in result.errors

    def test_intermediate_files(self):
        """Test result with intermediate files."""
        result = MinimizeResult(
            success=True,
            output_pdb="/path/to/minimized.pdb",
            input_pdb="/path/to/input.pdb",
            output_dir="/path/to/output",
            intermediate_files={
                "topology": "/path/to/topol.top",
                "structure": "/path/to/structure.gro",
            },
        )
        assert result.intermediate_files["topology"] == "/path/to/topol.top"
        assert result.intermediate_files["structure"] == "/path/to/structure.gro"
