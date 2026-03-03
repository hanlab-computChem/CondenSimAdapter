#!/usr/bin/env python3
"""Minimize command: AA softcore energy minimization."""

import sys
from pathlib import Path
from typing import List, Optional

import click
from click_option_group import optgroup

from ..shared import get_minimize_force_fields, validate_minimize_force_field


@click.command('minimize', context_settings={'help_option_names': ['-h', '--help']})
@click.option('--input', '-i', 'input_pdb', type=click.Path(exists=True), required=True,
              help='Input all-atom PDB (output from adapter backmap).')
@click.option('--input-file', '-f', type=click.Path(exists=True), default=None,
              help='Simulation config YAML (for component definitions).')
@click.option('--output', '-o', 'output_dir', type=click.Path(), default=None,
              help='Output directory (default: <system_name>_minimize/).')
@click.option('--force-field', '-ff', type=str, default='1-a99SBdisp',
              callback=validate_minimize_force_field,
              help='All-atom force field (e.g. 1-a99SBdisp, 7-amber99sb-ildn).')
@click.option('--gpu-id', '-g', type=int, default=0, show_default=True,
              help='GPU device ID.')
@click.option('--platform', type=str, default='CUDA', show_default=True,
              help='OpenMM platform: CUDA | OpenCL | CPU.')
@click.option('--tolerance', type=float, default=100.0, show_default=True,
              help='Energy tolerance kJ/(mol·nm).')
@click.option('--solvate', is_flag=True, default=False,
              help='Add explicit TIP3P solvent.')
@click.option('--no-disulfide', is_flag=True, default=False,
              help='Disable automatic disulfide bond detection.')
@click.option('--verbose', '-v', is_flag=True, default=False)
def minimize_command(
    input_pdb: str,
    input_file: Optional[str],
    output_dir: Optional[str],
    force_field: str,
    gpu_id: int,
    platform: str,
    tolerance: float,
    solvate: bool,
    no_disulfide: bool,
    verbose: bool,
):
    """\b
    Run three-stage softcore OpenMM energy minimization.

    \b
    Examples:
        adapter minimize -i FUS_LC_backmap/backmapped.pdb -f FUS_LC.yaml
        adapter minimize -i output.pdb -f config.yaml -ff 7-amber99sb-ildn
        adapter minimize -i output.pdb -f config.yaml --solvate
    """
    try:
        from ...minimize.minimizer import MinimizeSimulator, MinimizeConfig
        from ...minimize.config_loader import load_config_from_yaml, get_system_name
        from ...core.config import Component
    except ImportError as e:
        click.echo(f"Import error: {e}", err=True)
        sys.exit(1)

    click.echo(f"\n{'=' * 60}\nEnergy Minimization\n{'=' * 60}")
    click.echo(f"  Input PDB:   {input_pdb}")
    click.echo(f"  Force field: {force_field}")
    click.echo(f"  Platform:    {platform} (GPU {gpu_id})")

    # Load components from YAML
    system_name = Path(input_pdb).stem
    components = []
    if input_file:
        try:
            system_name, raw_comps = load_config_from_yaml(input_file)
            components = [Component.from_dict(c) for c in raw_comps]
            click.echo(f"  System:      {system_name} ({len(components)} component types)")
        except Exception as e:
            click.echo(f"  Warning: could not load config: {e}")

    out = output_dir or f"{system_name}_minimize"

    cfg = MinimizeConfig(
        forcefield_type=force_field,
        platform=platform,
        gpu_id=gpu_id,
        tolerance=tolerance,
        solvate=solvate,
        disable_disulfide=no_disulfide,
    )

    sim = MinimizeSimulator(cfg, components, system_name)
    try:
        result = sim.run(input_pdb, output_dir=out)
    except Exception as e:
        click.echo(f"\nMinimization error: {e}", err=True)
        if verbose:
            import traceback; traceback.print_exc()
        sys.exit(1)

    if result.success:
        click.echo(f"\n  Completed. Output PDB: {result.output_pdb}")
    else:
        click.echo(f"\n  Minimization failed.", err=True)
        for err in result.errors:
            click.echo(f"    {err}", err=True)
        sys.exit(1)
