#!/usr/bin/env python3
"""
CLI command for managing neural network models.

Usage:
    adapter models list        # List available models
    adapter models download    # Download all models
    adapter models download CalphaBasedModel  # Download specific model
    adapter models status      # Show download status
"""

from __future__ import annotations

import click

from ...backmap.cg2all.model import (
    MODEL_URLS,
    ensure_models,
    print_model_status,
)


@click.group("models")
def models_command():
    """Manage neural network models for backmapping."""
    pass


@models_command.command("list")
def list_cmd():
    """List all available models and their download status."""
    print_model_status()


@models_command.command("download")
@click.argument("model_names", nargs=-1, required=False)
def download_cmd(model_names):
    """Download model files.

    If no model names are specified, downloads all models.

    Examples:
        adapter models download                    # Download all
        adapter models download CalphaBasedModel   # Download specific
    """
    if not model_names:
        # Download all
        print("Downloading all models...")
        ensure_models()
    else:
        # Validate model names
        available = set(MODEL_URLS.keys())
        for name in model_names:
            if name not in available and not name.endswith(".ckpt"):
                name = name + ".ckpt"
            if name not in available:
                click.echo(f"Error: Unknown model '{name}'", err=True)
                click.echo(f"Available: {', '.join(available)}", err=True)
                raise click.Abort()

        print(f"Downloading: {', '.join(model_names)}")
        ensure_models(list(model_names))


@models_command.command("status")
def status_cmd():
    """Show detailed status of model files."""
    print_model_status()
