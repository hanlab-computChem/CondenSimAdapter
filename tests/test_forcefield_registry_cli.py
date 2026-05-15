from __future__ import annotations

import importlib
import shutil
from pathlib import Path

import CondenSimAdapter.src.minimize as minimize_module
from click.testing import CliRunner

from CondenSimAdapter.cli.commands_refactored.forcefield_command import forcefield_command
from CondenSimAdapter.cli.shared import validate_minimize_force_field
from CondenSimAdapter.forcefield import registry as registry_module
from CondenSimAdapter.forcefield.registry import ForceFieldRegistry


def _sample_ff_dir() -> Path:
    repo_root = Path(__file__).resolve().parents[1]
    new_ff_dir = repo_root / "new_ff" / "amber03.ff"
    if new_ff_dir.exists():
        return new_ff_dir
    return repo_root / "CondenSimAdapter" / "forcefield" / "amber03wsc.ff"


def _prepare_test_ff_dir(tmp_path) -> Path:
    src = _sample_ff_dir()
    dst = tmp_path / "test_custom.ff"
    shutil.copytree(src, dst)
    return dst


def _prepare_isolated_registry(tmp_path, monkeypatch) -> ForceFieldRegistry:
    monkeypatch.setattr(registry_module, "CUSTOM_FORCEFIELD_DIR", tmp_path / "custom")
    monkeypatch.setattr(
        registry_module, "CUSTOM_FORCEFIELD_INDEX", tmp_path / "user_forcefields.json"
    )
    return ForceFieldRegistry()


def test_registry_register_remove_and_minimize_recognition(tmp_path, monkeypatch) -> None:
    registry = _prepare_isolated_registry(tmp_path, monkeypatch)

    ff = registry.register_custom_force_field(
        ff_dir=str(_prepare_test_ff_dir(tmp_path)),
        family="AMBER",
        water_model="tip3p",
        solvate_cs="spc216",
        gbsa_mapping="AMBER99SB-ILDN",
        description="integration smoke ff",
    )
    assert ff.name == "a1"
    assert registry.get_force_field("a1") is not None
    assert registry.get_force_field_path("a1") is not None

    # Make minimize config validation use this isolated registry.
    monkeypatch.setattr(minimize_module, "FORCE_FIELD_REGISTRY", registry)
    config = minimize_module.MinimizeConfig(forcefield_type="a1")
    assert config.forcefield_type == "a1"

    removed = registry.remove_custom_force_field("a1")
    assert removed.name == "a1"
    assert registry.get_force_field("a1") is None


def test_registry_infers_gbsa_mapping_from_family(tmp_path, monkeypatch) -> None:
    registry = _prepare_isolated_registry(tmp_path, monkeypatch)
    ff = registry.register_custom_force_field(
        ff_dir=str(_prepare_test_ff_dir(tmp_path)),
        family="CHARMM",
        water_model="tip3p",
        solvate_cs="spc216",
        description="family mapping test",
    )
    assert ff.gbsa_mapping == "CHARMM36"


def test_builtin_forcefield_cannot_be_removed(tmp_path, monkeypatch) -> None:
    registry = _prepare_isolated_registry(tmp_path, monkeypatch)
    try:
        registry.remove_custom_force_field("1-a99SBdisp")
        assert False, "Expected built-in remove to fail"
    except ValueError as exc:
        assert "cannot be removed" in str(exc)


def test_forcefield_cli_add_list_remove_and_validate_minimize(tmp_path, monkeypatch) -> None:
    registry = _prepare_isolated_registry(tmp_path, monkeypatch)

    # Patch CLI/shared to use isolated registry instance.
    ff_cmd_module = importlib.import_module(
        "CondenSimAdapter.cli.commands_refactored.forcefield_command"
    )
    from CondenSimAdapter.cli import shared as shared_module

    monkeypatch.setattr(ff_cmd_module, "REGISTRY", registry)
    monkeypatch.setattr(shared_module, "REGISTRY", registry)

    runner = CliRunner()
    add_result = runner.invoke(
        forcefield_command,
        [
            "add",
            "--ff-dir",
            str(_prepare_test_ff_dir(tmp_path)),
            "--family",
            "AMBER",
            "--water-model",
            "tip3p",
            "--solvate-cs",
            "spc216",
        ],
    )
    assert add_result.exit_code == 0, add_result.output
    assert "a1" in add_result.output
    ff = registry.get_force_field("a1")
    assert ff is not None
    assert ff.gbsa_mapping == "AMBER99SB-ILDN"

    list_result = runner.invoke(forcefield_command, ["list"])
    assert list_result.exit_code == 0, list_result.output
    assert "a1" in list_result.output

    # Custom id is accepted by minimize force field validator.
    assert validate_minimize_force_field(None, None, "a1") == "a1"

    remove_result = runner.invoke(forcefield_command, ["remove", "a1"])
    assert remove_result.exit_code == 0, remove_result.output
    assert "Removed custom force field: a1" in remove_result.output
