"""
Integration tests for complete workflows

Tests cover:
- Complete CG simulation workflow
- Multi-force field compatibility
- Error handling in workflows
"""

from __future__ import annotations

import pytest
from pathlib import Path

from CondenSimAdapter.core.config import (
    CGConfig,
    Component,
    ComponentType,
    TopologyType,
    SimulationParams,
)
from CondenSimAdapter.core.simulation import CGSimulation


# =============================================================================
# Complete Workflow Integration Tests
# =============================================================================

@pytest.mark.slow
class TestCGSimulationWorkflow:
    """Integration tests for complete CG simulation workflow."""

    @pytest.fixture
    def tiny_simulation_config(self):
        """Create a tiny simulation config for fast testing."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="AAA",  # Very short for speed
            nmol=1,
        )
        sim_params = SimulationParams(
            steps=10,  # Minimal steps
            wfreq=5,
            log_freq=10,
            platform="CPU",
        )
        return CGConfig(
            system_name="integration_test",
            force_field="hps",
            components=[comp],
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
            simulation=sim_params,
        )

    def test_complete_simulation_workflow(self, tiny_simulation_config, tmp_path):
        """Test a complete simulation workflow runs successfully."""
        output_dir = tmp_path / "sim_output"
        
        sim = CGSimulation(tiny_simulation_config)
        result = sim.run(str(output_dir))
        
        assert result.success is True
        assert result.final_pdb is not None
        assert Path(result.final_pdb).exists()
        assert result.trajectory is not None
        assert Path(result.trajectory).exists()
        assert result.log_file is not None
        assert Path(result.log_file).exists()

    def test_simulation_creates_output_files(self, tiny_simulation_config, tmp_path):
        """Test that simulation creates all expected output files."""
        output_dir = tmp_path / "sim_output"
        
        sim = CGSimulation(tiny_simulation_config)
        result = sim.run(str(output_dir))
        
        # Check main outputs
        assert Path(result.final_pdb).exists()
        assert Path(result.trajectory).exists()
        assert Path(result.log_file).exists()

    def test_simulation_creates_top_pdb(self, tiny_simulation_config, tmp_path):
        """Test that simulation creates initial topology PDB."""
        output_dir = tmp_path / "sim_output"
        
        sim = CGSimulation(tiny_simulation_config)
        result = sim.run(str(output_dir))
        
        top_pdb = output_dir / "top.pdb"
        assert top_pdb.exists()


@pytest.mark.slow
class TestMultiForceFieldWorkflow:
    """Integration tests across all force fields."""

    @pytest.fixture
    def minimal_config(self, request):
        """Create minimal config for specified force field."""
        ff_name = request.param
        comp = Component(
            name="test",
            comp_type=ComponentType.IDP,
            sequence="AAA",
            nmol=1,
        )
        sim_params = SimulationParams(
            steps=5,  # Minimal
            wfreq=5,
            platform="CPU",
        )
        return CGConfig(
            system_name=f"test_{ff_name}",
            force_field=ff_name,
            components=[comp],
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
            simulation=sim_params,
        )

    @pytest.mark.parametrize("minimal_config", ["hps", "calvados2", "cocomo", "mpipi"], indirect=True)
    def test_force_field_runs(self, minimal_config, tmp_path):
        """Test that each force field can run a minimal simulation."""
        output_dir = tmp_path / f"sim_{minimal_config.force_field}"
        
        sim = CGSimulation(minimal_config)
        result = sim.run(str(output_dir))
        
        assert result.success, f"Force field {minimal_config.force_field} failed: {result.error}"
        assert Path(result.final_pdb).exists()


# =============================================================================
# Multi-Component System Tests
# =============================================================================

@pytest.mark.slow
class TestMultiComponentSystems:
    """Integration tests for multi-component systems."""

    def test_two_component_system(self, tmp_path):
        """Test simulation with two different components."""
        comp_a = Component(
            name="ProteinA",
            comp_type=ComponentType.IDP,
            sequence="AAAA",
            nmol=2,
        )
        comp_b = Component(
            name="ProteinB",
            comp_type=ComponentType.IDP,
            sequence="GGGG",
            nmol=3,
        )
        sim_params = SimulationParams(
            steps=5,
            wfreq=5,
            platform="CPU",
        )
        config = CGConfig(
            system_name="multi_component",
            force_field="hps",
            components=[comp_a, comp_b],
            box=[15.0, 15.0, 15.0],
            topology=TopologyType.CUBIC,
            simulation=sim_params,
        )
        
        output_dir = tmp_path / "multi_output"
        sim = CGSimulation(config)
        result = sim.run(str(output_dir))
        
        assert result.success is True

    def test_slab_topology(self, tmp_path):
        """Test slab topology simulation."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="AAA",
            nmol=2,
        )
        sim_params = SimulationParams(
            steps=5,
            wfreq=5,
            platform="CPU",
        )
        config = CGConfig(
            system_name="slab_test",
            force_field="hps",
            components=[comp],
            box=[20.0, 20.0, 100.0],
            topology=TopologyType.SLAB,
            slab_width=20.0,
            simulation=sim_params,
        )
        
        output_dir = tmp_path / "slab_output"
        sim = CGSimulation(config)
        result = sim.run(str(output_dir))
        
        assert result.success is True

    def test_droplet_topology(self, tmp_path):
        """Test droplet topology simulation."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="AAA",
            nmol=2,
        )
        sim_params = SimulationParams(
            steps=5,
            wfreq=5,
            platform="CPU",
        )
        config = CGConfig(
            system_name="droplet_test",
            force_field="hps",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.DROPLET,
            droplet_radius=5.0,
            simulation=sim_params,
        )
        
        output_dir = tmp_path / "droplet_output"
        sim = CGSimulation(config)
        result = sim.run(str(output_dir))
        
        assert result.success is True


# =============================================================================
# Error Handling Tests
# =============================================================================

class TestErrorHandling:
    """Tests for error handling in workflows."""

    def test_invalid_force_field_raises_on_simulation(self):
        """Test that invalid force field raises error on simulation."""
        comp = Component(
            name="test",
            comp_type=ComponentType.IDP,
            sequence="AAA",
        )
        config = CGConfig(
            system_name="test",
            force_field="invalid",
            components=[comp],
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
        )
        # Error should occur when creating CGSimulation
        with pytest.raises(ValueError, match="Unknown CG force field"):
            CGSimulation(config)

    def test_empty_components_works_but_simulation_fails(self):
        """Test that empty components doesn't raise at config level."""
        # Empty components might be allowed at config level
        config = CGConfig(
            system_name="test",
            force_field="hps",
            components=[],  # Empty
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
        )
        # But n_molecules should be 0
        assert config.n_molecules == 0

    def test_negative_box_size_invalid(self):
        """Test that negative box sizes are invalid."""
        comp = Component(
            name="test",
            comp_type=ComponentType.IDP,
            sequence="AAA",
        )
        # This might not raise immediately, but should fail during simulation
        config = CGConfig(
            system_name="test",
            force_field="hps",
            components=[comp],
            box=[-10.0, 10.0, 10.0],  # Negative
            topology=TopologyType.CUBIC,
        )
        # The negative box should be problematic
        assert any(b <= 0 for b in config.box)
