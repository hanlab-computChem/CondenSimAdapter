#!/usr/bin/env python3
"""
Backmap Command

Backmaps CG structure to all-atom representation.
"""

import sys

import click


@click.command('backmap', context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--input', '-i',
    type=click.Path(exists=True),
    required=True,
    help='Input: CG output directory (adapter cg) or PDB file (user provided)',
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
    type=click.Choice(['ResidueBasedModel', 'CalphaBasedModel', 'auto']),
    help='CG model type (overrides config, default: auto-detect)',
)
def backmap_command(input, input_file, output, model_type):
    """Backmap CG structure to all-atom representation.
    
    Uses the same config.yaml as CG stage, with optional backmap section.
    
    Examples:
        # adapter cg output (with explicit config.yaml)
        adapter backmap -i TDP43_CG -f config.yaml
        
        # adapter cg output (auto-find config.yaml in pwd)
        adapter backmap -i TDP43_CG
        
        # user provided PDB (requires config.yaml)
        adapter backmap -i my_structure.pdb -f config.yaml
        
        # Custom output
        adapter backmap -i TDP43_CG -f config.yaml -o ./results
    """
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
    if model_type:
        click.echo(f"  Model: {model_type}")
    
    # 创建 backmap config（CLI 参数优先）
    backmap_config = BackmapConfig()
    if model_type and model_type != 'auto':
        backmap_config.model_type = model_type
    if output:
        backmap_config.output_dir = output
    
    # 创建 simulator
    simulator = BackmapSimulator(backmap_config=backmap_config)
    
    # 执行 backmap（output_dir 参数会覆盖 backmap_config.output_dir）
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
