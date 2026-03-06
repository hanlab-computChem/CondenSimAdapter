#!/usr/bin/env python3
"""CG simulation command."""

import sys
from pathlib import Path
from typing import Optional

import click

from ..shared import CG_FORCE_FIELDS, normalize_cg_ff, validate_cg_force_field


@click.command('cg', context_settings={'help_option_names': ['-h', '--help']})
@click.option('--input-file', '-f', type=click.Path(exists=True), required=True,
              help='Input configuration YAML file.')
@click.option('--force-field', '-ff', type=str, default=None,
              callback=validate_cg_force_field,
              help=f'CG force field. Options: {", ".join(CG_FORCE_FIELDS)} (default: from YAML).')
@click.option('--output-dir', '-o', type=click.Path(), default=None,
              help='Output directory (default: <system_name>_CG/).')
@click.option('--overwrite', '-w', is_flag=True, default=False,
              help='Overwrite existing output directory.')
@click.option('--gpu-id', '-g', type=int, default=0,
              help='GPU device ID (default: 0).')
@click.option('--verbose', '-v', is_flag=True, default=False,
              help='Enable verbose output.')
def cg_command(
    input_file: str,
    force_field: Optional[str],
    output_dir: Optional[str],
    overwrite: bool,
    gpu_id: int,
    verbose: bool,
):
    """\b
    Run a coarse-grained simulation.

    \b
    Examples:
        adapter cg -f FUS_LC.yaml
        adapter cg -f FUS_LC.yaml -ff mpipi
        adapter cg -f FUS_LC.yaml -g 1
    """
    try:
        from ...core.config import CGConfig
        from ...core.simulation import CGSimulation
    except ImportError as e:
        click.echo(f"Import error: {e}", err=True)
        sys.exit(1)

    click.echo(f"\n{'=' * 60}\nCG Simulation\n{'=' * 60}")

    # Load config
    try:
        config = CGConfig.from_yaml(input_file)
    except Exception as e:
        click.echo(f"Failed to load config: {e}", err=True)
        if verbose:
            import traceback; traceback.print_exc()
        sys.exit(1)

    # Override force field if specified on CLI
    if force_field:
        config.force_field = normalize_cg_ff(force_field)
    else:
        config.force_field = normalize_cg_ff(config.force_field)

    click.echo(f"  System:      {config.system_name}")
    click.echo(f"  Force field: {config.force_field}")
    click.echo(f"  Requested:   {config.simulation.platform} (GPU {gpu_id})")

    # Determine output directory
    out = output_dir or f"{config.system_name}_CG"

    try:
        sim = CGSimulation(config)
        result = sim.run(out, gpu_id=gpu_id, overwrite=overwrite)
    except Exception as e:
        click.echo(f"\nSimulation error: {e}", err=True)
        if verbose:
            import traceback; traceback.print_exc()
        sys.exit(1)

    if result.success:
        click.echo(f"\n  Completed. Output: {result.output_dir}")
        click.echo(f"  Final PDB: {result.final_pdb}")
    else:
        click.echo(f"\n  Simulation failed: {result.error}", err=True)
        sys.exit(1)
