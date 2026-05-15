from pathlib import Path

from CondenSimAdapter.src.cg import (
    CGComponent,
    CGSimulationConfig,
    CGSimulator,
    ComponentType,
    TopologyType,
)


def test_cg_simulator_setup_marks_ready(tmp_path: Path) -> None:
    fasta_path = tmp_path / "seq.fasta"
    fasta_path.write_text(">A\nAAAA\n", encoding="utf-8")

    component = CGComponent(
        name="ProteinA",
        type=ComponentType.IDP,
        nmol=1,
        ffasta=str(fasta_path),
    )

    config = CGSimulationConfig(
        system_name="cg_smoke",
        components=[component],
        topol=TopologyType.SLAB,
    )

    simulator = CGSimulator(config)
    output_dir = tmp_path / "cg_output"
    result = simulator.setup(str(output_dir))

    assert simulator.is_setup is True
    assert Path(result["output_dir"]).resolve() == output_dir.resolve()
