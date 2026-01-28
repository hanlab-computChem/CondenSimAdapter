#!/usr/bin/env python3
"""
Init Command

Initializes a new CG simulation configuration template.
"""

import sys
from pathlib import Path
from typing import Optional

import click
import click_option_group as cog

from ..shared import (
    CG_FORCE_FIELDS,
    GEOMETRY_DEFAULTS,
    parse_component_pattern,
    parse_box,
    generate_yaml_with_comments,
)


@click.command('init', context_settings={'help_option_names': ['-h', '--help'], 'allow_extra_args': True})
@cog.optgroup.group(name="Component Options", help="Force field, component type, and molecule count settings")
@cog.optgroup.option(
    '--ff', '-ff',
    type=click.Choice(CG_FORCE_FIELDS),
    default='calvados',
    help='Force field (default: calvados)'
)
@cog.optgroup.option(
    '--type',
    type=click.Choice(['idp', 'mdp', 'mixed']),
    default='idp',
    help='Component type (default: idp). Overridden by --components if provided.'
)
@cog.optgroup.option(
    '--components', '-c',
    type=str,
    default=None,
    help='Component pattern (e.g., "IIIMII" = 3 IDP + 1 MDP + 2 IDP). Overrides --type.'
)
@cog.optgroup.option(
    '--nmol', '-n',
    type=int,
    default=10,
    help='Number of molecules per component (default: 10). For mixed, use "-n A B".'
)
@cog.optgroup.group(name="Topology Options", help="Geometry type and box dimensions")
@cog.optgroup.option(
    '--topol', '-tp',
    type=click.Choice(['grid', 'slab', 'droplet']),
    default='grid',
    help='Topology/geometry type (default: grid)'
)
@cog.optgroup.option(
    '--box', '-b',
    nargs=3,
    type=float,
    default=None,
    help='Box dimensions (e.g., "-b 20 20 20" for grid/slab)'
)
@cog.optgroup.option(
    '--radius', '-r',
    type=float,
    default=None,
    help='Droplet radius in nm (only needed for droplet geometry, overrides --box)'
)
@cog.optgroup.group(name="Simulation Options", help="Simulation parameters: time, temperature, and ionic strength")
@cog.optgroup.option(
    '--time', '-t',
    type=float,
    default=1000.0,
    help='Simulation time in nanoseconds (default: 1000 ns)'
)
@cog.optgroup.option(
    '--temperature', '-T',
    type=float,
    default=310.0,
    help='Temperature in Kelvin (default: 310 K)'
)
@cog.optgroup.option(
    '--ionic', '-I',
    type=float,
    default=0.15,
    help='Ionic strength in Molar (default: 0.15 M)'
)
@cog.optgroup.group(name="Output Options", help="Output file and system name settings")
@cog.optgroup.option(
    '--name',
    type=str,
    default='my_simulation',
    help='System name and output filename (default: my_simulation)'
)
@click.pass_context
def init_command(ctx: click.Context, name: str, ff: str, type: str, topol: str, nmol: int,
                 time: float, components: Optional[str], box: Optional[tuple],
                 radius: Optional[float], temperature: float, ionic: float):
    """\b
    Initialize a new CG simulation configuration template.

    \b
    Examples:
        adapter init --name my_project                       # Custom system name
        adapter init --topol slab                            # Slab geometry [10, 10, 40]
        adapter init --topol slab -b 10 10 40                # Slab with custom box
        adapter init --topol droplet -r 15                   # Droplet with radius 15 nm (box: [30, 30, 30])
        adapter init -c IIIMII --nmol 10                     # 3 IDP + 1 MDP + 2 IDP
        adapter init --time 5000                             # 5000 ns simulation
        adapter init --temperature 293 --ionic 0.2           # Custom parameters
    """
    # Parse extra args for mixed nmol (supports: -n A B)
    extra_args = list(ctx.args)
    if extra_args:
        if components:
            click.echo("Error: Unexpected extra arguments.", err=True)
            sys.exit(1)
        if type == 'mixed' and len(extra_args) == 1:
            try:
                nmol_mixed = (nmol, int(extra_args[0]))
                extra_args = []
            except ValueError:
                click.echo(f"Error: Invalid nmol value '{extra_args[0]}'.", err=True)
                sys.exit(1)
        else:
            click.echo(f"Error: Got unexpected extra argument(s): {' '.join(extra_args)}", err=True)
            sys.exit(1)
    else:
        nmol_mixed = None

    # Get geometry defaults
    geom_defaults = GEOMETRY_DEFAULTS.get(topol, GEOMETRY_DEFAULTS['grid'])
    
    # Handle radius parameter for droplet geometry
    if radius is not None:
        if topol != 'droplet':
            click.echo(f"Warning: --radius is only applicable for droplet topology. Ignoring.", err=True)
            box_values = parse_box(box, topol, geom_defaults['box'])
        else:
            # For droplet: box = [2*r, 2*r, 2*r]
            box_values = [2 * radius, 2 * radius, 2 * radius]
    else:
        # Parse box (use provided or default based on geometry)
        try:
            box_values = parse_box(box, topol, geom_defaults['box'])
        except ValueError as e:
            click.echo(f"Error: {e}", err=True)
            sys.exit(1)
    
    # Calculate steps from time (1 ns = 100,000 steps, 1 step = 10 fs)
    steps = int(time * 100000)
    
    # Build component list
    component_list = []
    if components:
        # Use component pattern if provided (overrides --type)
        try:
            component_list = parse_component_pattern(components, nmol)
            component_note = f"Pattern: {components} ({len(component_list)} components)"
        except ValueError as e:
            click.echo(f"Error: {e}", err=True)
            sys.exit(1)
    else:
        # Use --type (backward compatibility)
        if type == 'idp':
            component_list.append({
                'name': 'protein_A',
                'type': 'IDP',
                'nmol': nmol,
                'ffasta': 'input/protein_A.fasta',
            })
            component_note = "IDP (requires FASTA file)"
        elif type == 'mdp':
            component_list.append({
                'name': 'protein_A',
                'type': 'MDP',
                'nmol': nmol,
                'fpdb': 'input/protein_A.pdb',
                'restraint': True,
                'restraint_type': 'harmonic',
                'charge_termini': 'both',
            })
            component_note = "MDP (requires PDB file, add fdomains in YAML)"
        else:  # mixed
            if nmol_mixed:
                nmol_idp, nmol_mdp = nmol_mixed
            else:
                nmol_idp = nmol
                nmol_mdp = nmol
            component_list.extend([
                {
                    'name': 'protein_A',
                    'type': 'IDP',
                    'nmol': nmol_idp,
                    'ffasta': 'input/protein_A.fasta',
                },
                {
                    'name': 'protein_B',
                    'type': 'MDP',
                    'nmol': nmol_mdp,
                    'fpdb': 'input/protein_B.pdb',
                    'restraint': True,
                    'restraint_type': 'harmonic',
                    'charge_termini': 'both',
                },
            ])
            component_note = "Mixed IDP + MDP"
    
    # Create config (without platform in simulation params)
    from ...src import CGSimulationConfig, CGComponent, SimulationParams, TopologyType
    config = CGSimulationConfig(
        system_name=name,
        box=box_values,
        temperature=temperature,
        ionic=ionic,
        topol=TopologyType(topol),
        components=[CGComponent.from_dict(c) for c in component_list],
        simulation=SimulationParams(
            steps=steps,
            wfreq=5000,
            verbose=True
        )
    )
    
    # Output path (current directory)
    config_file = Path(f'{name}.yaml')
    
    # Ensure output directory exists
    config_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate YAML with enhanced comments
    yaml_content = generate_yaml_with_comments(config, topol, geom_defaults['description'], time, component_list, ff, radius)
    
    # Write to file
    with open(config_file, 'w') as f:
        f.write(yaml_content)
    
    click.echo(f"\n{'=' * 60}")
    click.echo(f"Configuration template created successfully!")
    click.echo(f"{'=' * 60}")
    click.echo(f"\n  File: {config_file}")
    click.echo(f"  System: {name}")
    click.echo(f"  Force field: {ff}")
    click.echo(f"  Topology: {topol} - {geom_defaults['description']}")
    click.echo(f"  Box: [{', '.join(f'{v:.1f}' for v in box_values)}] nm")
    click.echo(f"  Temperature: {temperature} K")
    click.echo(f"  Ionic: {ionic} M")
    click.echo(f"  Time: {time} ns ({steps:,} steps)")
    click.echo(f"  Components: {len(component_list)} ({component_note})")
    if type == 'mixed' and nmol_mixed:
        click.echo(f"  Molecules per component: IDP={nmol_mixed[0]}, MDP={nmol_mixed[1]}")
    else:
        click.echo(f"  Molecules per component: {nmol}")
    
    click.echo(f"\n  Next steps:")
    click.echo(f"    1. Edit {config_file} (especially fdomains for MDP components)")
    click.echo(f"    2. Add your input files (FASTA/PDB) to 'input/' directory")
    click.echo(f"    3. Run: adapter cg -f {config_file} [options]")
    click.echo(f"       (use 'adapter cg -h' to see available options)")
    click.echo()
