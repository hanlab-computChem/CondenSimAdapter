"""Smoke tests for the top-level CLI entrypoint."""

from __future__ import annotations

from click.testing import CliRunner

from CondenSimAdapter.cli import main


def test_top_level_help_does_not_require_heavy_runtime_dependencies():
    result = CliRunner().invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "CORE COMMANDS" in result.output
    assert "models" in result.output


def test_command_help_paths_load_lazily():
    runner = CliRunner()

    for command in ("backmap", "forcefield", "models"):
        result = runner.invoke(main, [command, "--help"])
        assert result.exit_code == 0, result.output
        assert "Usage:" in result.output


def test_models_status_does_not_download_models():
    result = CliRunner().invoke(main, ["models", "status"])

    assert result.exit_code == 0
    assert "Available Models" in result.output
