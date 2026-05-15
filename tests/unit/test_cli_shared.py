"""
Unit tests for cli/shared.py

Tests cover:
- Force field validation
- Component pattern parsing
- Box parsing
"""

from __future__ import annotations

import click
import pytest

from CondenSimAdapter.cli.shared import (
    CG_FORCE_FIELDS,
    GEOMETRY_DEFAULTS,
    MINIMIZE_FORCE_FIELDS,
    parse_box,
    parse_component_pattern,
    validate_cg_force_field,
    validate_minimize_force_field,
)

# =============================================================================
# Force Field Constants Tests
# =============================================================================


class TestForceFieldConstants:
    """Tests for force field constants."""

    def test_cg_force_fields_not_empty(self):
        """Test CG force fields list is not empty."""
        assert len(CG_FORCE_FIELDS) > 0

    def test_cg_force_fields_contains_calvados(self):
        """Test CG force fields contains calvados."""
        assert "calvados" in CG_FORCE_FIELDS

    def test_minimize_force_fields_not_empty(self):
        """Test minimize force fields list is not empty."""
        assert len(MINIMIZE_FORCE_FIELDS) > 0

    def test_geometry_defaults_has_slab(self):
        """Test geometry defaults has slab."""
        assert "slab" in GEOMETRY_DEFAULTS

    def test_geometry_defaults_has_droplet(self):
        """Test geometry defaults has droplet."""
        assert "droplet" in GEOMETRY_DEFAULTS


# =============================================================================
# CG Force Field Validation Tests
# =============================================================================


class TestValidateCGForceField:
    """Tests for validate_cg_force_field function."""

    def test_valid_calvados(self):
        """Test valid calvados."""
        result = validate_cg_force_field(None, None, "calvados")
        assert result == "calvados"

    def test_valid_hps_urry(self):
        """Test valid hps_urry (alias for hps)."""
        result = validate_cg_force_field(None, None, "hps_urry")
        assert result == "hps"

    def test_valid_cocomo(self):
        """Test valid cocomo."""
        result = validate_cg_force_field(None, None, "cocomo")
        assert result == "cocomo"

    def test_valid_mpipi_recharged(self):
        """Test valid mpipi_recharged (alias for mpipi)."""
        result = validate_cg_force_field(None, None, "mpipi_recharged")
        assert result == "mpipi"

    def test_invalid_force_field_raises(self):
        """Test invalid force field raises BadParameter."""
        with pytest.raises(click.BadParameter):
            validate_cg_force_field(None, None, "invalid_ff")

    def test_case_insensitive(self):
        """Test force field validation is case insensitive."""
        result = validate_cg_force_field(None, None, "CALVADOS")
        assert result == "calvados"


# =============================================================================
# Minimize Force Field Validation Tests
# =============================================================================


class TestValidateMinimizeForceField:
    """Tests for validate_minimize_force_field function."""

    def test_valid_amber_ff(self):
        """Test valid AMBER force field."""
        result = validate_minimize_force_field(None, None, "1-a99SBdisp")
        assert result == "1-a99SBdisp"

    def test_valid_charmm_ff(self):
        """Test valid CHARMM force field."""
        result = validate_minimize_force_field(None, None, "9-charmm36m")
        assert result == "9-charmm36m"

    def test_invalid_force_field_raises(self):
        """Test invalid force field raises BadParameter."""
        with pytest.raises(click.BadParameter):
            validate_minimize_force_field(None, None, "invalid_ff")


# =============================================================================
# Component Pattern Parsing Tests
# =============================================================================


class TestParseComponentPattern:
    """Tests for parse_component_pattern function."""

    def test_single_component(self):
        """Test single component pattern."""
        result = parse_component_pattern("I", nmol=1)
        assert len(result) == 1
        assert result[0]["type"] == "IDP"

    def test_multiple_idps(self):
        """Test multiple IDP pattern."""
        result = parse_component_pattern("III", nmol=2)
        assert len(result) == 3
        for comp in result:
            assert comp["type"] == "IDP"
            assert comp["nmol"] == 2

    def test_mixed_pattern(self):
        """Test mixed IDP/MDP pattern."""
        result = parse_component_pattern("IIMII", nmol=1)
        assert len(result) == 5
        assert result[0]["type"] == "IDP"
        assert result[1]["type"] == "IDP"
        assert result[2]["type"] == "MDP"
        assert result[3]["type"] == "IDP"
        assert result[4]["type"] == "IDP"


# =============================================================================
# Box Parsing Tests
# =============================================================================


class TestParseBox:
    """Tests for parse_box function (accepts tuple from nargs=3)."""

    def test_grid_default_when_none(self):
        """Test grid default box when None provided."""
        result = parse_box(None, "grid", GEOMETRY_DEFAULTS["grid"]["box"])
        assert result == GEOMETRY_DEFAULTS["grid"]["box"]

    def test_slab_default_when_none(self):
        """Test slab default box when None provided."""
        result = parse_box(None, "slab", GEOMETRY_DEFAULTS["slab"]["box"])
        assert result == GEOMETRY_DEFAULTS["slab"]["box"]

    def test_droplet_default_when_none(self):
        """Test droplet default box when None provided."""
        result = parse_box(None, "droplet", GEOMETRY_DEFAULTS["droplet"]["box"])
        assert result == GEOMETRY_DEFAULTS["droplet"]["box"]

    def test_custom_box_tuple(self):
        """Test custom box with tuple."""
        result = parse_box((20.0, 30.0, 40.0), "grid", [10.0, 10.0, 10.0])
        assert result == [20.0, 30.0, 40.0]

    def test_droplet_uses_x_for_all(self):
        """Test droplet uses x value for all dimensions."""
        result = parse_box((15.0, 20.0, 30.0), "droplet", [10.0, 10.0, 10.0])
        assert result == [15.0, 15.0, 15.0]  # x value used for all
