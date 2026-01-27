#!/usr/bin/env python3
"""
Init Command

Initializes a new CG simulation configuration template.
"""

import sys
from pathlib import Path
from typing import Optional

import click

from ..shared import (
    CG_FORCE_FIELDS,
    GEOMETRY_DEFAULTS,
    parse_component_pattern,
    parse_box,
    generate_yaml_with_comments,
)


@click.command('init', context_settings={'help_option_names': ['-h', '--help']})
@click.argument('name', type=str, required=False)
@click.option(
    '--ff', '-ff',
    type=click.Choice(CG_FORCE_FIELDS),
    default='calvados',
    help=f'Force field. Available: {", ".join(CG_FORCE_FIELDS)} (default: calvados)'
)
@click.option(
    '--type',
    type=click.Choice(['idp', 'mdp', 'mixed']),
    default='idp',
    help='Component type (default: idp). Overridden by --components if provided.'
)
@click.option(
    '--topol', '-tp',
    type=click.Choice(['grid', 'slab', 'droplet']),
    default='grid',
    help='Topology/geometry type (default: grid)'
)
@click.option(
    '--output', '-o',
    type=click.Path(),
    default='.',
    help='Output directory (default: current directory)'
)
@click.option(
    '--nmol', '-n',
    type=int,
    default=10,
    help='Number of molecules per component (default: 10)'
)
@click.option(
    '--time', '-t',
    type=float,
    default=1000.0,
    help='Simulation time in nanoseconds (default: 1000 ns)'
)
@click.option(
    '--components', '-c',
    type=str,
    default=None,
    help='Component pattern (e.g., "IIIMII" = 3 IDP + 1 MDP + 2 IDP). Overrides --type.'
)
@click.option(
    '--box', '-b',
    nargs=3,
    type=float,
    default=None,
    help='Box dimensions (e.g., "-b 20 20 20" for grid, "-b 15 15 15" for droplet radius)'
)
@click.option(
    '--temperature', '-T',
    type=float,
    default=310.0,
    help='Temperature in Kelvin (default: 310 K)'
)
@click.option(
    '--ionic', '-I',
    type=float,
    default=0.15,
    help='Ionic strength in Molar (default: 0.15 M)'
)
def init_command(name: str, ff: str, type: str, topol: str, output: str, nmol: int, 
                 time: float, components: Optional[str], box: Optional[tuple], 
                 temperature: float, ionic: float):
    """\b
    Initialize a new CG simulation configuration template.

    \b
    NAME: System name (optional, defaults to 'my_simulation').

    \b
    Examples:
        adapter init                                    # Grid geometry, 1000ns
        adapter init my_project                         # Custom name
        adapter init --topol slab                       # Slab geometry [10, 10, 40]
        adapter init --topol slab -b 10 10 40           # Slab with custom box
        adapter init --topol droplet -b 15 15 15        # Droplet with radius 15 nm
        adapter init -c IIIMII --nmol 10                # 3 IDP + 1 MDP + 2 IDP
        adapter init --time 5000                        # 5000 ns simulation
        adapter init --temperature 293 --ionic 0.2      # Custom parameters
    """
    # Set default name if not provided
    if not name:
        name = 'my_simulation'
    
    # Get geometry defaults
    geom_defaults = GEOMETRY_DEFAULTS.get(topol, GEOMETRY_DEFAULTS['grid'])
    
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
            component_list.extend([
                {
                    'name': 'protein_A',
                    'type': 'IDP',
                    'nmol': nmol,
                    'ffasta': 'input/protein_A.fasta',
                },
                {
                    'name': 'protein_B',
                    'type': 'MDP',
                    'nmol': nmol,
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
    
    # Output path
    output_path = Path(output)
    config_file = output_path / f'{name}.yaml'
    
    # Ensure output directory exists
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate YAML with enhanced comments
    yaml_content = generate_yaml_with_comments(config, topol, geom_defaults['description'], time, component_list, ff)
    
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
    click.echo(f"  Molecules per component: {nmol}")
    
    click.echo(f"\n  Next steps:")
    click.echo(f"    1. Edit {config_file} (especially fdomains for MDP components)")
    click.echo(f"    2. Add your input files (FASTA/PDB) to 'input/' directory")
    click.echo(f"    3. Run: adapter cg -f {config_file} [options]")
    click.echo(f"       (use 'adapter cg -h' to see available options)")
    click.echo()
