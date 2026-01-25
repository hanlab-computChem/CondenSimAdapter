from pathlib import Path

from CondenSimAdapter.src.backmap import BackmapSimulator


def test_cg2all_converts_sample_pdb(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    input_pdb = repo_root / "tests/data/final.pdb"

    cg_dir = tmp_path / "cg_output"
    cg_dir.mkdir()
    (cg_dir / "final.pdb").write_bytes(input_pdb.read_bytes())

    output_dir = tmp_path / "backmap_out"
    simulator = BackmapSimulator()
    result = simulator.run(str(cg_dir), output_dir=str(output_dir))

    assert result.success, f"Backmap failed: {result.errors}"
    assert result.output_pdb
    assert Path(result.output_pdb).exists()
    assert Path(result.output_pdb).stat().st_size > 0
