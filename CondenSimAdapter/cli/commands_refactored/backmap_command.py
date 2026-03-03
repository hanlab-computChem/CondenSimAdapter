#!/usr/bin/env python3
"""Backmap command: CG -> all-atom."""

import sys
from pathlib import Path
from typing import Optional

import click

from ...backmap.backmapper import SUPPORTED_MODELS


@click.command('backmap', context_settings={'help_option_names': ['-h', '--help']})
@click.option('--input', '-i', 'input_path', type=click.Path(exists=True), required=True,
              help='Input CG PDB file (final.pdb from adapter cg).')
@click.option('--input-file', '-f', type=click.Path(exists=True), default=None,
              help='Simulation config YAML (used to detect topology type for slab centring).')
@click.option('--output', '-o', 'output_dir', type=click.Path(), default=None,
              help='Output directory (default: <system_name>_backmap/).')
@click.option('--model-type', '-m', type=str, default='CalphaBasedModel',
              show_default=True,
              help=f'cg2all model type. Options: {", ".join(SUPPORTED_MODELS)}.')
@click.option('--device', '-d', type=str, default='cpu', show_default=True,
              help='Compute device: cpu | cuda | cuda:0 etc.')
@click.option('--verbose', '-v', is_flag=True, default=False)
def backmap_command(
    input_path: str,
    input_file: Optional[str],
    output_dir: Optional[str],
    model_type: str,
    device: str,
    verbose: bool,
):
    """\b
    Convert a CG structure to all-atom via the cg2all neural network.

    \b
    Examples:
        adapter backmap -i FUS_LC_CG/final.pdb
        adapter backmap -i FUS_LC_CG/final.pdb -m Martini3
        adapter backmap -i FUS_LC_CG/final.pdb -f FUS_LC.yaml -d cuda
    """
    try:
        from ...backmap.backmapper import Backmapper
        from ...minimize.config_loader import load_config_from_yaml
    except ImportError as e:
        click.echo(f"Import error: {e}", err=True)
        sys.exit(1)

    click.echo(f"\n{'=' * 60}\nBackmapping\n{'=' * 60}")
    click.echo(f"  Input:      {input_path}")
    click.echo(f"  Model type: {model_type}")
    click.echo(f"  Device:     {device}")

    # Detect topology type for optional slab centring
    topology_type = None
    system_name = Path(input_path).stem
    if input_file:
        try:
            name, _ = load_config_from_yaml(input_file)
            system_name = name
            import yaml
            with open(input_file) as f:
                raw = yaml.safe_load(f)
            topology_type = raw.get("topol") or raw.get("topology")
        except Exception as e:
            if verbose:
                click.echo(f"  Warning: could not parse config: {e}")

    out = output_dir or f"{system_name}_backmap"

    backmapper = Backmapper()
    result = backmapper.run(
        cg_pdb=input_path,
        output_dir=out,
        model_type=model_type,
        device=device,
        topology_type=topology_type,
    )

    if result.success:
        click.echo(f"\n  Completed. Output PDB: {result.output_pdb}")
    else:
        click.echo(f"\n  Backmapping failed: {result.error}", err=True)
        if verbose and result.error:
            click.echo(result.error, err=True)
        sys.exit(1)
