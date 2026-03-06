"""
Unit tests for core/simulation.py

Tests cover:
- CGSimulation initialization
- Platform selection logic
- System creation (force assembly)
- PDB I/O helpers
"""

from __future__ import annotations

import numpy as np
import pytest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

pytest.importorskip("openmm")

import openmm as mm
import openmm.app as app
import openmm.unit as unit

from CondenSimAdapter.core.simulation import (
    CGSimulation,
    _resolve_platform,
    _save_pdb,
    _save_final_pdb,
)
from CondenSimAdapter.core.config import (
    CGConfig,
    Component,
    ComponentType,
    TopologyType,
    SimulationParams,
)


# =============================================================================
# CGSimulation Initialization Tests
# =============================================================================

class TestCGSimulationInit:
    """Tests for CGSimulation initialization."""

    @pytest.fixture
    def simple_config(self):
        """Create a simple CGConfig for testing."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="GSMASAS",
            nmol=1,
        )
        return CGConfig(
            system_name="test",
            force_field="calvados2",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.CUBIC,
        )

    def test_init_creates_forcefield(self, simple_config):
        """Test that initialization creates forcefield."""
        sim = CGSimulation(simple_config)
        assert sim.config == simple_config
        assert sim.ff is not None
        assert sim._result is None

    def test_init_uses_resolved_force_field(self, simple_config):
        """Test that resolved force field is used."""
        # Change to bare 'calvados' to test resolution
        simple_config.force_field = "calvados"
        sim = CGSimulation(simple_config)
        # Should resolve to calvados2 for all-IDP system
        from CondenSimAdapter.core.forcefield.calvados import CalvadosFF
        assert isinstance(sim.ff, CalvadosFF)


# =============================================================================
# Platform Resolution Tests
# =============================================================================

class TestResolvePlatform:
    """Tests for _resolve_platform function."""

    def test_cpu_fallback(self):
        """Test that CPU is returned when CUDA is not available."""
        plat, props, name = _resolve_platform("CUDA", 0)
        # Should fall back to CPU if CUDA not available
        assert name in ["CPU", "Reference", "CUDA:0"]

    def test_explicit_cpu(self):
        """Test explicit CPU platform request."""
        plat, props, name = _resolve_platform("CPU", 0)
        assert name == "CPU"

    def test_returns_platform_object(self):
        """Test function returns OpenMM Platform object."""
        plat, props, name = _resolve_platform("CPU", 0)
        assert isinstance(plat, mm.Platform)

    def test_gpu_id_in_properties(self):
        """Test that GPU ID is in properties for CUDA."""
        # This will likely fall back, but test the structure
        plat, props, name = _resolve_platform("CUDA", 1)
        if name.startswith("CUDA"):
            assert props.get("DeviceIndex") == "1"


# =============================================================================
# System Creation Tests
# =============================================================================

class TestCreateSystem:
    """Tests for _create_system method."""

    @pytest.fixture
    def simple_config(self):
        """Create a simple CGConfig."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="AAAA",
            nmol=1,
        )
        return CGConfig(
            system_name="test",
            force_field="hps",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.CUBIC,
        )

    @pytest.fixture
    def chain_meta(self):
        """Simple chain metadata."""
        return [{
            "name": "FUS",
            "start": 0,
            "end": 4,
            "sequence": "AAAA",
            "folded_domains": [],
        }]

    @pytest.fixture
    def positions(self):
        """Simple positions."""
        return np.array([
            [0, 0, 0],
            [0.38, 0, 0],
            [0.76, 0, 0],
            [1.14, 0, 0],
        ], dtype=np.float64)

    @pytest.fixture
    def topology(self, chain_meta, positions):
        """Create OpenMM topology."""
        from CondenSimAdapter.core.topology import build_topology
        top, _ = build_topology(chain_meta, positions, [20.0, 20.0, 20.0])
        return top

    def test_system_has_correct_num_particles(self, simple_config, topology, chain_meta, positions):
        """Test system has correct number of particles."""
        sim = CGSimulation(simple_config)
        system = sim._create_system(topology, chain_meta, positions, simple_config)
        assert system.getNumParticles() == 4

    def test_system_has_periodic_box(self, simple_config, topology, chain_meta, positions):
        """Test system has periodic box vectors."""
        sim = CGSimulation(simple_config)
        system = sim._create_system(topology, chain_meta, positions, simple_config)
        box = system.getDefaultPeriodicBoxVectors()
        assert box is not None
        assert box[0][0].value_in_unit(unit.nanometer) == 20.0

    def test_system_has_forces(self, simple_config, topology, chain_meta, positions):
        """Test system has forces added."""
        sim = CGSimulation(simple_config)
        system = sim._create_system(topology, chain_meta, positions, simple_config)
        # Should have at least: bonds, nonbonded (2), CM motion remover
        assert system.getNumForces() >= 3

    def test_hps_no_angle_force(self, simple_config, topology, chain_meta, positions):
        """Test HPS doesn't add angle force."""
        sim = CGSimulation(simple_config)
        system = sim._create_system(topology, chain_meta, positions, simple_config)
        # Check no HarmonicAngleForce
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            assert not isinstance(force, mm.HarmonicAngleForce)

    def test_cocomo_adds_angle_force(self, topology, chain_meta, positions):
        """Test COCOMO adds angle force."""
        comp = Component(name="FUS", comp_type=ComponentType.IDP, sequence="AAAA", nmol=1)
        config = CGConfig(
            system_name="test",
            force_field="cocomo",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.CUBIC,
        )
        sim = CGSimulation(config)
        system = sim._create_system(topology, chain_meta, positions, config)
        # Check for HarmonicAngleForce
        has_angle = False
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), mm.HarmonicAngleForce):
                has_angle = True
                break
        assert has_angle

    def test_enm_added_for_folded_domains(self, topology, positions):
        """Test ENM is added when folded domains exist."""
        chain_meta_with_domains = [{
            "name": "MDP1",
            "start": 0,
            "end": 4,
            "sequence": "AAAA",
            "folded_domains": [(1, 4)],
        }]
        comp = Component(name="MDP1", comp_type=ComponentType.MDP, sequence="AAAA", nmol=1)
        config = CGConfig(
            system_name="test",
            force_field="calvados3",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.CUBIC,
        )
        sim = CGSimulation(config)
        system = sim._create_system(topology, chain_meta_with_domains, positions, config)
        # Should have ENM force
        force_types = [type(system.getForce(i)).__name__ for i in range(system.getNumForces())]
        assert "HarmonicBondForce" in force_types or "CustomBondForce" in force_types

    def test_droplet_adds_confinement_force(self, topology, chain_meta, positions):
        """Test droplet topology adds confinement force."""
        comp = Component(name="FUS", comp_type=ComponentType.IDP, sequence="AAAA", nmol=1)
        config = CGConfig(
            system_name="test",
            force_field="hps",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.DROPLET,
            droplet_radius=5.0,
        )
        sim = CGSimulation(config)
        system = sim._create_system(topology, chain_meta, positions, config)
        # Check for CustomExternalForce (droplet confinement)
        has_confinement = False
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), mm.CustomExternalForce):
                has_confinement = True
                break
        assert has_confinement

    def test_cm_motion_remover_added(self, simple_config, topology, chain_meta, positions):
        """Test CMMotionRemover is added."""
        sim = CGSimulation(simple_config)
        system = sim._create_system(topology, chain_meta, positions, simple_config)
        has_cm = False
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), mm.CMMotionRemover):
                has_cm = True
                break
        assert has_cm


# =============================================================================
# Run Method Tests (Mocked)
# =============================================================================

class TestRunMethod:
    """Tests for run method with mocked OpenMM."""

    @pytest.fixture
    def simple_config(self):
        """Create a simple CGConfig with very short simulation."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="AAAA",
            nmol=1,
        )
        sim_params = SimulationParams(
            steps=100,
            wfreq=50,
            log_freq=100,
            platform="CPU",
        )
        return CGConfig(
            system_name="test",
            force_field="hps",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.CUBIC,
            simulation=sim_params,
        )

    def test_raises_if_output_exists(self, simple_config, tmp_path):
        """Test that run raises FileExistsError if output exists."""
        output_dir = tmp_path / "existing_output"
        output_dir.mkdir()
        (output_dir / "dummy.txt").write_text("dummy")
        
        sim = CGSimulation(simple_config)
        with pytest.raises(FileExistsError, match="already exists"):
            sim.run(str(output_dir), overwrite=False)

    def test_creates_output_directory(self, simple_config, tmp_path):
        """Test that run creates output directory."""
        output_dir = tmp_path / "new_output"
        
        sim = CGSimulation(simple_config)
        # Use mock to avoid actual simulation
        with patch.object(sim, '_run_pipeline') as mock_pipeline:
            mock_pipeline.return_value = Mock(
                success=True,
                output_dir=str(output_dir),
                final_pdb="final.pdb",
                trajectory="traj.xtc",
                log_file="sim.log",
            )
            sim.run(str(output_dir))
        
        assert output_dir.exists()

    def test_returns_result_object(self, simple_config, tmp_path):
        """Test that run returns SimulationResult."""
        output_dir = tmp_path / "output"
        
        sim = CGSimulation(simple_config)
        with patch.object(sim, '_run_pipeline') as mock_pipeline:
            mock_result = Mock(
                success=True,
                output_dir=str(output_dir),
                final_pdb="final.pdb",
                trajectory="traj.xtc",
                log_file="sim.log",
                elapsed_seconds=1.0,
            )
            mock_pipeline.return_value = mock_result
            result = sim.run(str(output_dir))
        
        assert result.success is True
        assert result.output_dir == str(output_dir)

    def test_error_handling(self, simple_config, tmp_path):
        """Test that errors are caught and returned in result."""
        output_dir = tmp_path / "output"
        
        sim = CGSimulation(simple_config)
        with patch.object(sim, '_run_pipeline') as mock_pipeline:
            mock_pipeline.side_effect = RuntimeError("Test error")
            result = sim.run(str(output_dir))
        
        assert result.success is False
        assert "Test error" in result.error
        assert result.elapsed_seconds > 0

    def test_overwrite_flag(self, simple_config, tmp_path):
        """Test that overwrite flag allows reusing directory."""
        output_dir = tmp_path / "existing_output"
        output_dir.mkdir()
        
        sim = CGSimulation(simple_config)
        with patch.object(sim, '_run_pipeline') as mock_pipeline:
            mock_pipeline.return_value = Mock(
                success=True,
                output_dir=str(output_dir),
                final_pdb="final.pdb",
                trajectory="traj.xtc",
                log_file="sim.log",
            )
            # Should not raise with overwrite=True
            result = sim.run(str(output_dir), overwrite=True)
            assert result.success is True


# =============================================================================
# PDB I/O Tests
# =============================================================================

class TestSavePDB:
    """Tests for _save_pdb helper."""

    @pytest.fixture
    def simple_topology(self):
        """Create a simple topology."""
        top = app.Topology()
        chain = top.addChain()
        res = top.addResidue("ALA", chain)
        top.addAtom("CA", app.element.carbon, res)
        return top

    @pytest.fixture
    def positions(self):
        """Simple positions with units."""
        return np.array([[1.0, 2.0, 3.0]]) * unit.nanometer

    def test_creates_file(self, simple_topology, positions, tmp_path):
        """Test that file is created."""
        output_file = tmp_path / "output.pdb"
        _save_pdb(simple_topology, positions, str(output_file), [10.0, 10.0, 10.0])
        assert output_file.exists()

    def test_file_has_atoms(self, simple_topology, positions, tmp_path):
        """Test that file contains ATOM records."""
        output_file = tmp_path / "output.pdb"
        _save_pdb(simple_topology, positions, str(output_file), [10.0, 10.0, 10.0])
        content = output_file.read_text()
        assert "ATOM" in content
        assert "CA" in content
        assert "ALA" in content


class TestSaveFinalPDB:
    """Tests for _save_final_pdb helper."""

    def test_creates_file(self, tmp_path):
        """Test that final PDB is created."""
        # Create a mock simulation with proper position format
        mock_sim = Mock()
        mock_state = Mock()
        # Return positions as numpy array with units (matching OpenMM format)
        positions = np.array([[1.0, 2.0, 3.0]]) * unit.nanometer
        mock_state.getPositions.return_value = positions
        mock_sim.context.getState.return_value = mock_state
        
        # Create simple topology
        top = app.Topology()
        chain = top.addChain()
        res = top.addResidue("ALA", chain)
        top.addAtom("CA", app.element.carbon, res)
        
        output_file = tmp_path / "final.pdb"
        _save_final_pdb(mock_sim, top, [10.0, 10.0, 10.0], str(output_file))
        
        assert output_file.exists()


# =============================================================================
# Integration Smoke Tests
# =============================================================================

@pytest.mark.slow
class TestSimulationIntegration:
    """Integration tests that actually run short simulations."""

    @pytest.fixture
    def tiny_config(self):
        """Create a tiny config for fast testing."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="AAA",  # Very short
            nmol=1,
        )
        sim_params = SimulationParams(
            steps=10,  # Very short
            wfreq=5,
            log_freq=10,
            platform="CPU",
        )
        return CGConfig(
            system_name="tiny_test",
            force_field="hps",
            components=[comp],
            box=[10.0, 10.0, 10.0],
            topology=TopologyType.CUBIC,
            simulation=sim_params,
        )

    def test_full_simulation_runs(self, tiny_config, tmp_path):
        """Test that a full simulation runs successfully."""
        output_dir = tmp_path / "sim_output"
        
        sim = CGSimulation(tiny_config)
        result = sim.run(str(output_dir))
        
        assert result.success is True
        assert result.final_pdb is not None
        assert Path(result.final_pdb).exists()
        assert result.trajectory is not None
        assert Path(result.trajectory).exists()

    def test_simulation_outputs_log(self, tiny_config, tmp_path):
        """Test that simulation generates log file."""
        output_dir = tmp_path / "sim_output"
        
        sim = CGSimulation(tiny_config)
        result = sim.run(str(output_dir))
        
        assert result.log_file is not None
        assert Path(result.log_file).exists()

    def test_multiple_force_fields_run(self, tmp_path):
        """Test that all force fields can run a short simulation."""
        for ff_name in ["calvados2", "hps", "cocomo", "mpipi"]:
            comp = Component(
                name="test",
                comp_type=ComponentType.IDP,
                sequence="AAA",
                nmol=1,
            )
            sim_params = SimulationParams(
                steps=5,
                wfreq=5,
                platform="CPU",
            )
            config = CGConfig(
                system_name=f"test_{ff_name}",
                force_field=ff_name,
                components=[comp],
                box=[10.0, 10.0, 10.0],
                topology=TopologyType.CUBIC,
                simulation=sim_params,
            )
            
            output_dir = tmp_path / f"sim_{ff_name}"
            sim = CGSimulation(config)
            result = sim.run(str(output_dir))
            
            assert result.success, f"Force field {ff_name} failed: {result.error}"
