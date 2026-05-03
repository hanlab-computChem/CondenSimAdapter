#!/usr/bin/env python3
"""Minimize command: all-atom energy minimization."""

import sys
from pathlib import Path
from typing import Optional

import click

from ..shared import validate_minimize_force_field


def _ff_help() -> str:
    """Build a compact force-field help string listing all available FFs by number."""
    try:
        from ...forcefield.registry import BUILTIN_FORCE_FIELDS
        names = [ff.name for ff in BUILTIN_FORCE_FIELDS]
        ff_list = "  ".join(names)
        return (
            f"Force field (name or number). Available:\n  {ff_list}\n"
            "  Input the number (e.g. -ff 1) or full name (e.g. -ff 1-a99SBdisp)."
        )
    except Exception:
        return "Force field name or number (e.g. -ff 1 or -ff 1-a99SBdisp)."


@click.command('minimize', context_settings={'help_option_names': ['-h', '--help']})
@click.option('--input', '-i', 'input_pdb', type=click.Path(exists=True), required=True,
              help='Input all-atom PDB (from backmap step).')
@click.option('--input-file', '-f', type=click.Path(exists=True), default=None,
              help='Simulation config YAML (to load component definitions).')
@click.option('--output', '-o', 'output_dir', type=click.Path(), default=None,
              help='Output directory (default: {system_name}_minimize).')
@click.option('--force-field', '-ff', 'force_field', type=str, default='1',
              callback=validate_minimize_force_field,
              help=_ff_help())
@click.option('--platform', '-p', type=str, default='CUDA', show_default=True,
              help='OpenMM platform: CUDA | CPU | OpenCL.')
@click.option('--gpu-id', '-g', type=int, default=0, show_default=True,
              help='GPU device index.')
@click.option('--solvate', is_flag=True, default=False,
              help='Add explicit water and ions after minimization.')
@click.option('--salt-conc', type=float, default=0.15, show_default=True,
              help='Ion concentration in M (only used with --solvate).')
@click.option('--level', '-l', type=click.Choice(['high', 'medium', 'low']),
              default='medium', show_default=True,
              help='Softcore optimization level.')
@click.option('--his-type', 'his_type', type=click.Choice(['0', '1']), default='1',
              show_default=True,
              help='Histidine protonation for pdb2gmx -his.  0 = HID (delta, neutral)  1 = HIE (epsilon, neutral).')
@click.option('--no-disulfide', 'disable_disulfide', is_flag=True, default=False,
              help='Pass -ss to pdb2gmx to disable disulfide-bond detection.')
@click.option('--box-type', '-bt', type=click.Choice(['dodecahedron', 'cubic', 'octahedron']),
              default=None, help='Box shape for droplet solvation (gmx editconf -bt).')
@click.option('--box-distance', '-dd', type=float, default=2.0, show_default=True,
              help='Water shell thickness in nm (gmx editconf -d).')
@click.option('--verbose', '-v', is_flag=True, default=False)
def minimize_command(
    input_pdb: str,
    input_file: Optional[str],
    output_dir: Optional[str],
    force_field: str,
    platform: str,
    gpu_id: int,
    solvate: bool,
    salt_conc: float,
    level: str,
    his_type: str,
    disable_disulfide: bool,
    box_type: Optional[str],
    box_distance: float,
    verbose: bool,
):
    """\b
    All-atom energy minimization using softcore potential.

    \b
    Examples:
        adapter minimize -f FUS_LC.yaml -i FUS_LC_backmap/backmapped.pdb -ff 1
        adapter minimize -f FUS_LC.yaml -i FUS_LC_backmap/backmapped.pdb -ff 1 --solvate
        adapter minimize -f FUS_LC.yaml -i FUS_LC_backmap/backmapped.pdb -ff 2 --solvate -bt dodecahedron -dd 2
    """
    try:
        from ...minimize.minimizer import MinimizeConfig, MinimizeSimulator
        from ...minimize.config_loader import load_config_from_yaml
        from ...core.config import Component
    except ImportError as e:
        click.echo(f"Import error: {e}", err=True)
        sys.exit(1)

    his_type_int = int(his_type)
    his_label = {0: "HID (delta, neutral)", 1: "HIE (epsilon, neutral)"}.get(his_type_int, str(his_type_int))

    click.echo(f"\n{'=' * 60}\nEnergy Minimization\n{'=' * 60}")
    click.echo(f"  Input PDB:   {input_pdb}")
    click.echo(f"  Force field: {force_field}")
    click.echo(f"  Platform:    {platform} (GPU {gpu_id})")
    click.echo(f"  HIS type:    {his_type_int} — {his_label}")
    if disable_disulfide:
        click.echo(f"  Disulfide:   disabled (-ss)")

    # Load components from YAML
    components = []
    system_name = Path(input_pdb).stem
    if input_file:
        try:
            name, raw_comps = load_config_from_yaml(input_file)
            system_name = name
            components = [Component.from_dict(c) for c in raw_comps]
            click.echo(f"  System:      {system_name} ({len(components)} component types)")
        except Exception as e:
            click.echo(f"  Warning: could not load config: {e}")

    out = output_dir or f"{system_name}_minimize"

    config = MinimizeConfig(
        forcefield_type=force_field,
        platform=platform,
        gpu_id=gpu_id,
        solvate=solvate,
        ion_concentration=salt_conc,
        his_type=his_type_int,
        disable_disulfide=disable_disulfide,
        droplet_box_type=box_type,
        droplet_distance=box_distance,
    )

    sim = MinimizeSimulator(config=config, components=components, system_name=system_name)

    try:
        result = sim.run(input_pdb, output_dir=out)
    except Exception as e:
        click.echo(f"\nMinimization error: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    if result.success:
        click.echo(f"\n  Completed. Output PDB: {result.output_pdb}")
    else:
        click.echo(f"\nMinimization failed")
        for err in result.errors:
            click.echo(err)
        sys.exit(1)
