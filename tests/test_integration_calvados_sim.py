from pathlib import Path

import pytest
from CondenSimAdapter.src.cg import (
    CGSimulationConfig,
    CGSimulator,
    ComputePlatform,
    SimulationParams,
)


@pytest.mark.slow
def test_calvados_runs_short_simulation(tmp_path: Path) -> None:
    pytest.importorskip("openmm")

    repo_root = Path(__file__).resolve().parents[1]
    config_path = repo_root / "tests/data/config.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))

    # Keep the same inputs, but shrink for a quick smoke run.
    config.simulation = SimulationParams(
        steps=10,
        wfreq=2,
        platform=ComputePlatform.CPU,
        verbose=False,
    )
    if config.components:
        config.components[0].nmol = 1

    simulator = CGSimulator(config)
    simulator.setup(str(tmp_path / "calvados_output"))
    result = simulator.run_calvados(gpu_id=0)

    assert result.success, f"Calvados failed: {result.errors}"
    assert result.structure is not None
    assert Path(result.structure).exists()
