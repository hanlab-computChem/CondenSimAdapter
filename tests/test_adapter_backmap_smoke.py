from pathlib import Path

from CondenSimAdapter.src.backmap import BackmapSimulator, SourceType


def test_backmap_detects_ms2_cg_dir(tmp_path: Path) -> None:
    cg_dir = tmp_path / "sample_CG"
    cg_dir.mkdir()
    (cg_dir / "final.pdb").write_text("HEADER    CG OUTPUT\n", encoding="utf-8")

    sim = BackmapSimulator()
    assert sim.detect_source_type(str(cg_dir)) == SourceType.MS2_CG


def test_backmap_detects_user_pdb_file(tmp_path: Path) -> None:
    pdb_path = tmp_path / "input.pdb"
    pdb_path.write_text("HEADER    USER PDB\n", encoding="utf-8")

    sim = BackmapSimulator()
    assert sim.detect_source_type(str(pdb_path)) == SourceType.USER_PROVIDED


def test_backmap_model_checkpoint_exists() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    model_path = repo_root / "CondenSimAdapter/extern/ms2_cg2all/model/Martini3-FIX.ckpt"
    assert model_path.exists()
