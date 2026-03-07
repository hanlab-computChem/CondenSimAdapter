"""
Unit tests for core/config.py

Tests cover:
- Component creation and sequence resolution
- CGConfig loading from dict/YAML
- Force field auto-resolution (calvados -> 2/3)
- SimulationParams defaults and overrides
"""

from __future__ import annotations

import pytest
import yaml
from pathlib import Path

from CondenSimAdapter.core.config import (
    Component,
    ComponentType,
    CGConfig,
    SimulationParams,
    TopologyType,
    _read_fasta,
    _parse_fdomains,
)


# =============================================================================
# Component Tests
# =============================================================================

class TestComponent:
    """Tests for Component dataclass."""

    def test_component_idp_basic(self):
        """Test basic IDP component creation."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="GSMASAS",
            nmol=5,
        )
        assert comp.name == "FUS"
        assert comp.comp_type == ComponentType.IDP
        assert comp.get_sequence() == "GSMASAS"
        assert comp.nmol == 5
        assert comp.folded_domains == []

    def test_component_from_dict_idp(self):
        """Test Component creation from dictionary."""
        d = {
            "name": "ProteinA",
            "type": "idp",
            "nmol": 3,
            "sequence": "AAAA",
        }
        comp = Component.from_dict(d)
        assert comp.name == "ProteinA"
        assert comp.comp_type == ComponentType.IDP
        assert comp.nmol == 3
        assert comp.get_sequence() == "AAAA"

    def test_component_from_dict_mdp(self):
        """Test Component creation for MDP with folded domains."""
        d = {
            "name": "MDP1",
            "type": "mdp",
            "nmol": 1,
            "sequence": "AAAAAAAAAA",
            "fdomains": [[1, 5], [8, 10]],
        }
        comp = Component.from_dict(d)
        assert comp.comp_type == ComponentType.MDP
        assert comp.folded_domains == [(1, 5), (8, 10)]

    def test_component_get_sequence_priority(self, tmp_path: Path):
        """Test sequence resolution priority: sequence > fasta > pdb."""
        # Write test files
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">Test\nSEQFROMFASTA\n")
        
        # Direct sequence should be used first
        comp = Component(
            name="Test",
            comp_type=ComponentType.IDP,
            sequence="DIRECT",
            fasta_path=str(fasta),
        )
        assert comp.get_sequence() == "DIRECT"
        
        # Without sequence, fasta should be used
        comp2 = Component(
            name="Test",
            comp_type=ComponentType.IDP,
            fasta_path=str(fasta),
        )
        assert comp2.get_sequence() == "SEQFROMFASTA"

    def test_component_no_sequence_raises(self):
        """Test that Component without any sequence source raises ValueError."""
        comp = Component(name="Empty", comp_type=ComponentType.IDP)
        with pytest.raises(ValueError, match="no sequence"):
            comp.get_sequence()


# =============================================================================
# Helper Function Tests
# =============================================================================

class TestReadFasta:
    """Tests for _read_fasta helper."""

    def test_read_fasta_single_sequence(self, tmp_path: Path):
        """Test reading a single sequence from FASTA."""
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">Protein1\nAAAA\nBBBB\n>Protein2\nCCCC\n")
        
        result = _read_fasta(str(fasta), "Protein1")
        assert result == "AAAABBBB"
        
        result2 = _read_fasta(str(fasta), "Protein2")
        assert result2 == "CCCC"

    def test_read_fasta_first_sequence_if_no_name(self, tmp_path: Path):
        """Test reading first sequence when no name specified."""
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">AnyName\nMYSEQUENCE\n")
        
        result = _read_fasta(str(fasta))
        assert result == "MYSEQUENCE"

    def test_read_fasta_not_found_raises(self, tmp_path: Path):
        """Test that missing sequence name raises ValueError."""
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">Other\nSEQ\n")
        
        with pytest.raises(ValueError, match="not found"):
            _read_fasta(str(fasta), "Missing")


class TestParseFdomains:
    """Tests for _parse_fdomains helper."""

    def test_parse_list_of_pairs(self):
        """Test parsing list of [start, end] pairs."""
        raw = [[1, 50], [80, 130]]
        result = _parse_fdomains(raw)
        assert result == [(1, 50), (80, 130)]

    def test_parse_string_ranges(self):
        """Test parsing string format like '1-50, 80-130'."""
        raw = "1-50, 80-130"
        result = _parse_fdomains(raw)
        assert result == [(1, 50), (80, 130)]

    def test_parse_string_with_semicolon(self):
        """Test parsing string with semicolon separators."""
        raw = "1-50; 80-130"
        result = _parse_fdomains(raw)
        assert result == [(1, 50), (80, 130)]

    def test_parse_dict_format(self):
        """Test parsing dict format with component name."""
        raw = {"CompName": [[1, 50], [80, 130]]}
        result = _parse_fdomains(raw)
        assert result == [(1, 50), (80, 130)]

    def test_parse_empty(self):
        """Test parsing empty/None input."""
        assert _parse_fdomains([]) == []
        assert _parse_fdomains(None) == []
        assert _parse_fdomains("") == []


# =============================================================================
# SimulationParams Tests
# =============================================================================

class TestSimulationParams:
    """Tests for SimulationParams dataclass."""

    def test_default_values(self):
        """Test default parameter values."""
        params = SimulationParams()
        assert params.steps == 100_000_000
        assert params.dt == 0.01
        assert params.wfreq == 10_000
        assert params.log_freq == 1_000_000
        assert params.friction == 0.01
        assert params.platform == "CUDA"
        assert params.gpu_id == 0

    def test_from_dict_overrides(self):
        """Test that from_dict correctly overrides defaults."""
        d = {"steps": 1000, "dt": 0.02, "platform": "CPU", "gpu_id": 1}
        params = SimulationParams.from_dict(d)
        assert params.steps == 1000
        assert params.dt == 0.02
        assert params.wfreq == 10_000  # default
        assert params.platform == "CPU"
        assert params.gpu_id == 1


# =============================================================================
# CGConfig Tests
# =============================================================================

class TestCGConfig:
    """Tests for CGConfig dataclass."""

    def test_basic_config(self):
        """Test basic CGConfig creation."""
        comp = Component(name="A", comp_type=ComponentType.IDP, sequence="AAAA")
        config = CGConfig(
            system_name="test_system",
            force_field="calvados2",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.CUBIC,
        )
        assert config.system_name == "test_system"
        assert config.n_molecules == 1
        assert config.box == [20.0, 20.0, 20.0]

    def test_n_molecules_multi(self):
        """Test n_molecules counts all component copies."""
        comp1 = Component(name="A", comp_type=ComponentType.IDP, sequence="AAA", nmol=5)
        comp2 = Component(name="B", comp_type=ComponentType.IDP, sequence="BBB", nmol=3)
        config = CGConfig(
            system_name="multi",
            force_field="hps",
            components=[comp1, comp2],
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
        )
        assert config.n_molecules == 8

    def test_resolved_force_field_calvados_auto_idp(self):
        """Test that 'calvados' resolves to calvados2 for all-IDP system."""
        comp = Component(name="A", comp_type=ComponentType.IDP, sequence="AAAA")
        config = CGConfig(
            system_name="test",
            force_field="calvados",
            components=[comp],
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
        )
        assert config.resolved_force_field == "calvados2"

    def test_resolved_force_field_calvados_auto_mdp(self):
        """Test that 'calvados' resolves to calvados3 when MDP present."""
        comp = Component(name="A", comp_type=ComponentType.MDP, sequence="AAAA")
        config = CGConfig(
            system_name="test",
            force_field="calvados",
            components=[comp],
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
        )
        assert config.resolved_force_field == "calvados3"

    def test_resolved_force_field_explicit(self):
        """Test that explicit force field names are preserved."""
        for ff in ["calvados2", "calvados3", "hps", "cocomo", "mpipi"]:
            config = CGConfig(
                system_name="test",
                force_field=ff,
                components=[Component(name="A", comp_type=ComponentType.IDP, sequence="AAA")],
                box=[10.0, 10.0, 10.0],
                topology=TopologyType.CUBIC,
            )
            assert config.resolved_force_field == ff

    def test_get_component(self):
        """Test retrieving component by name."""
        comp_a = Component(name="ProteinA", comp_type=ComponentType.IDP, sequence="AAA")
        comp_b = Component(name="ProteinB", comp_type=ComponentType.IDP, sequence="BBB")
        config = CGConfig(
            system_name="test",
            force_field="hps",
            components=[comp_a, comp_b],
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
        )
        assert config.get_component("ProteinA") == comp_a
        assert config.get_component("ProteinB") == comp_b
        assert config.get_component("Missing") is None

    def test_from_dict_basic(self):
        """Test CGConfig creation from dictionary."""
        d = {
            "system_name": "from_dict_test",
            "force_field": "hps",
            "components": [{"name": "A", "type": "idp", "sequence": "AAAA"}],
            "box": [15.0, 15.0, 15.0],
            "topology": "cubic",
            "temperature": 310.0,
            "ionic_strength": 0.1,
        }
        config = CGConfig.from_dict(d)
        assert config.system_name == "from_dict_test"
        assert config.force_field == "hps"
        assert config.temperature == 310.0
        assert config.ionic_strength == 0.1
        assert config.topology == TopologyType.CUBIC

    def test_from_dict_legacy_keys(self):
        """Test CGConfig accepts legacy key names."""
        d = {
            "sysname": "legacy_test",
            "ff": "cocomo",
            "components": [{"name": "A", "type": "idp", "sequence": "AAA"}],
            "box": [10.0, 10.0, 10.0],
            "topol": "slab",
            "temp": 298.0,
            "ionic": 0.15,
        }
        config = CGConfig.from_dict(d)
        assert config.system_name == "legacy_test"
        assert config.force_field == "cocomo"
        assert config.temperature == 298.0
        assert config.topology == TopologyType.SLAB

    def test_from_dict_grid_to_cubic(self):
        """Test that 'grid' topology is converted to 'cubic'."""
        d = {
            "system_name": "test",
            "components": [{"name": "A", "type": "idp", "sequence": "AAA"}],
            "box": [10.0, 10.0, 10.0],
            "topol": "grid",
        }
        config = CGConfig.from_dict(d)
        assert config.topology == TopologyType.CUBIC

    def test_from_yaml(self, tmp_path: Path):
        """Test CGConfig loading from YAML file."""
        config_data = {
            "system_name": "yaml_test",
            "force_field": "mpipi",
            "components": [{"name": "FUS", "type": "idp", "nmol": 10, "sequence": "GSMASAS"}],
            "box": [20.0, 20.0, 20.0],
            "topology": "droplet",
            "droplet_radius": 8.0,
            "droplet_k": 2.0,
        }
        yaml_path = tmp_path / "config.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(config_data, f)
        
        config = CGConfig.from_yaml(str(yaml_path))
        assert config.system_name == "yaml_test"
        assert config.force_field == "mpipi"
        assert config.topology == TopologyType.DROPLET
        assert config.droplet_radius == 8.0
        assert config.droplet_k == 2.0

    def test_slab_specific_defaults(self):
        """Test slab-specific default values (0.6 * Lz)."""
        config = CGConfig(
            system_name="test",
            force_field="calvados",
            components=[Component(name="A", comp_type=ComponentType.IDP, sequence="AAA")],
            box=[50.0, 50.0, 200.0],
            topology=TopologyType.SLAB,
        )
        assert config.slab_width is None  # default is None (will use 0.6 * Lz at build time)

    def test_droplet_specific_defaults(self):
        """Test droplet-specific default values."""
        config = CGConfig(
            system_name="test",
            force_field="hps",
            components=[Component(name="A", comp_type=ComponentType.IDP, sequence="AAA")],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.DROPLET,
        )
        assert config.droplet_radius is None  # will use default calculation
        assert config.droplet_k == 1.0  # default
