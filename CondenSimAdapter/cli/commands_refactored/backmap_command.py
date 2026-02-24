#!/usr/bin/env python3
"""
Backmap Command

Backmaps CG structure to all-atom representation.
"""

import sys
from pathlib import Path

import click


@click.command('backmap', context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--input', '-i',
    type=click.Path(exists=True),
    required=False,
    default=None,
    help='Input: CG output directory (default: {system_name}_cg) or PDB file (user provided)',
)
@click.option(
    '--input-file', '-f',
    type=click.Path(exists=True),
    help='Configuration YAML (same as CG stage, same flag as adapter cg)',
)
@click.option(
    '--output', '-o',
    type=click.Path(),
    help='Output directory (default: {system_name}_backmap)',
)
@click.option(
    '--model-type', '-m',
    type=click.Choice(['ResidueBasedModel', 'CalphaBasedModel', 'martini2', 'martini3', 'auto']),
    help='CG model type (overrides config, default: auto-detect)',
)
def backmap_command(input, input_file, output, model_type):
    """\b
    Backmap CG structure to all-atom representation.

    \b
    Uses the same config.yaml as CG stage, with optional backmap section.

    \b
    Examples:
        adapter backmap -f config.yaml                  # Auto-find {system_name}_cg
        adapter backmap -i TDP43_CG -f config.yaml      # With explicit CG directory
        adapter backmap -i my_structure.pdb -f config.yaml  # User provided PDB
        adapter backmap -i TDP43_CG -f config.yaml -o ./results  # Custom output
    """
    # Map CLI model names to library model names
    model_name_map = {
        'martini2': 'Martini',
        'martini3': 'Martini3',
    }
    # Keep original for display
    model_type_display = model_type
    if model_type in model_name_map:
        model_type = model_name_map[model_type]

    from ...src.backmap import BackmapSimulator, BackmapConfig
    from ...src import CGSimulationConfig
    from ...src.pdb2gmx_utils import load_config_from_yaml
    
    # If input_file is not provided, try to find config.yaml in current directory
    if input_file is None:
        import os
        if os.path.exists('config.yaml'):
            input_file = 'config.yaml'
        elif os.path.exists('system.yaml'):
            input_file = 'system.yaml'
    
    if input_file is None:
        click.echo(f"Error: No configuration file found. Please provide -f config.yaml", err=True)
        sys.exit(1)
    
    # Load system_name from config to determine default input/output paths
    try:
        system_name, _ = load_config_from_yaml(input_file)
    except Exception as e:
        click.echo(f"Error: Failed to load configuration file: {e}", err=True)
        sys.exit(1)
    
    # If input is not provided, default to {system_name}_CG
    if input is None:
        input = f"{system_name}_CG"
    
    if not Path(input).exists():
        click.echo(f"Error: Input path '{input}' does not exist.", err=True)
        click.echo(f"  Provide -i flag or ensure '{input}' directory exists.", err=True)
        sys.exit(1)
    
    # If output is not provided, default to {system_name}_backmap
    if output is None:
        output = f"{system_name}_backmap"
    from ...src.backmap import BackmapSimulator, BackmapConfig
    from ...src import CGSimulationConfig
    
    click.echo(f"\n{'=' * 60}")
    click.echo(f"Backmap CG to All-Atom")
    click.echo(f"{'=' * 60}")
    
    click.echo(f"\n  Input: {input}")
    if input_file:
        click.echo(f"  Config: {input_file}")
    if output:
        click.echo(f"  Output: {output}")
    click.echo(f"  Device: cpu (fixed)")
    if model_type_display:
        click.echo(f"  Model: {model_type_display}")
    
    # Create backmap config (CLI parameters take priority)
    backmap_config = BackmapConfig()
    if model_type and model_type != 'auto':
        backmap_config.model_type = model_type
    if output:
        backmap_config.output_dir = output
    
    # Create simulator
    simulator = BackmapSimulator(backmap_config=backmap_config)
    
    # Execute backmap (output_dir parameter overrides backmap_config.output_dir)
    click.echo(f"\n[1/2] Preparing input...")
    try:
        result = simulator.run(input_path=input, config_path=input_file, output_dir=output)
        
        if result.success:
            click.echo(f"  ✓ Backmap completed")
            click.echo(f"  Input PDB: {result.input_pdb}")
            click.echo(f"  Output PDB: {result.output_pdb}")
            click.echo(f"  Model type: {result.model_type}")
            click.echo(f"\n  ✓ Success!")
        else:
            click.echo(f"  ✗ Backmap failed:")
            for error in result.errors:
                click.echo(f"    - {error}")
            sys.exit(1)
            
    except Exception as e:
        click.echo(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    click.echo()
