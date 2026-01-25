#!/usr/bin/env python3
"""
CG Command

Runs coarse-grained simulation.
"""

import os
import sys
from pathlib import Path
from typing import Optional

import click

from ..shared import (
    CG_FORCE_FIELDS,
    FORCE_FIELD_TO_RUNNER,
    validate_cg_force_field,
)


@click.command('cg', context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--input-file', '-f',
    type=click.Path(exists=True),
    help='Input configuration file (YAML)',
)
@click.option(
    '--force-field', '-ff',
    type=str,
    default=None,
    callback=validate_cg_force_field,
    help=f'Force field to use. Available: {", ".join(CG_FORCE_FIELDS)} (default: read from YAML)',
)
@click.option(
    '--output-dir', '-o',
    type=click.Path(),
    help='Output directory (default: same as input file directory)',
)
@click.option(
    '--dry-run', '-d',
    is_flag=True,
    default=False,
    help='Generate config files without running simulation',
)
@click.option(
    '--overwrite', '-w',
    is_flag=True,
    default=False,
    help='Overwrite existing output directory',
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    default=False,
    help='Enable verbose output',
)
@click.option(
    '--gpu-id', '-g',
    type=int,
    default=0,
    help='GPU device ID (default: 0)',
)
@click.option(
    '--continue-from', '-c',
    type=click.Path(exists=True),
    default=None,
    help='Continue simulation from a PDB coordinate file (Calvados only)',
)
def cg_command(
    input_file: str,
    force_field: str,
    output_dir: Optional[str],
    dry_run: bool,
    overwrite: bool,
    verbose: bool,
    gpu_id: int,
    continue_from: str,
):
    """
    Run coarse-grained simulation.
    
    INPUT_FILE is the configuration YAML file.
    
    Examples:
        adapter cg -f config.yaml                 # Run with calvados (GPU 0)
        adapter cg -f config.yaml -g 1            # Run on GPU 1
        adapter cg -f config.yaml -ff mpipi_recharged  # Run with Mpipi-Recharged
        adapter cg -f config.yaml -ff cocomo      # Run with COCOMO
        adapter cg -f config.yaml --dry-run       # Generate config only
        adapter cg -f config.yaml -o ./results    # Custom output directory
        adapter cg -f config.yaml -c coords.pdb   # Continue from PDB (Calvados only)
    """
    # Handle input file
    if input_file is None:
        click.echo(click.style("Error: No input file specified.", fg='red'))
        click.echo("Use 'adapter cg -f config.yaml' or 'adapter --help' for usage.")
        sys.exit(1)
    
    input_path = Path(input_file).resolve()
    
    # Determine output directory
    if output_dir is None:
        output_dir = str(input_path.parent)
    else:
        output_dir = str(Path(output_dir).resolve())

    click.echo(f"\n{'=' * 60}")
    click.echo(f"CG Simulation")
    click.echo(f"{'=' * 60}")

    click.echo(f"\n  Input file: {input_path}")
    click.echo(f"  GPU ID: {gpu_id}")
    click.echo(f"  Mode: {'Dry run' if dry_run else 'Full simulation'}")
    if continue_from:
        click.echo(f"  Continue from: {continue_from}")

    # Load configuration first to get system_name
    click.echo(f"\n[1/4] Loading configuration...")
    try:
        from ...src import CGSimulationConfig
        config = CGSimulationConfig.from_yaml(str(input_path))
        click.echo(f"  ✓ Loaded: {config.system_name}")
        click.echo(f"  ✓ Components: {len(config.components)}")
        click.echo(f"  ✓ Total molecules: {config.total_molecules()}")
    except Exception as e:
        click.echo(f"  ✗ Failed to load configuration: {e}")
        sys.exit(1)
    
    # Handle force field: read from YAML, override if specified in CLI
    yaml_force_field = config.force_field or 'calvados'
    if force_field is None:
        # Use force field from YAML
        force_field = yaml_force_field
        click.echo(f"  ✓ Force field (from YAML): {force_field}")
    else:
        # Force field specified in CLI
        if force_field != yaml_force_field:
            click.echo(f"  ⚠ WARNING: Force field mismatch!")
            click.echo(f"    YAML: {yaml_force_field}")
            click.echo(f"    CLI:  {force_field}")
            click.echo(f"  → Using CLI force field: {force_field}")
            click.echo(f"  → Updating YAML configuration...")
            # Update config and save back to YAML
            config.force_field = force_field
            try:
                import yaml
                with open(str(input_path), 'r') as f:
                    yaml_data = yaml.safe_load(f)
                yaml_data['force_field'] = force_field
                with open(str(input_path), 'w') as f:
                    yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
                click.echo(f"  ✓ YAML updated with force field: {force_field}")
            except Exception as e:
                click.echo(f"  ⚠ Failed to update YAML: {e}")
        else:
            click.echo(f"  ✓ Force field: {force_field}")

    # Calculate the final output directory with _CG suffix
    final_output_dir = os.path.join(output_dir, f"{config.system_name}_CG")
    click.echo(f"  Output dir: {final_output_dir}")
    
    # Validate configuration
    click.echo(f"\n[2/4] Validating configuration...")
    errors = config.validate()
    if errors:
        click.echo(f"  ✗ Validation failed:")
        for error in errors:
            click.echo(f"    - {error}")
        sys.exit(1)
    click.echo(f"  ✓ Configuration valid")
    
    # Validate --continue-from flag
    if continue_from and force_field != 'calvados':
        click.echo(f"  ✗ --continue-from is only supported with Calvados force field.")
        click.echo(f"    Current force field: {force_field}")
        sys.exit(1)

    # Setup simulator with final output directory
    click.echo(f"\n[3/4] Setting up simulation...")
    try:
        from ...src import CGSimulator

        sim = CGSimulator(config)

        # Setup output directory (with _CG suffix)
        sim.setup(final_output_dir, overwrite=overwrite)
        click.echo(f"  ✓ Output directory: {final_output_dir}")

    except FileExistsError:
        click.echo(f"  ✗ Output directory exists: {final_output_dir}")
        click.echo(f"    Use --overwrite to replace.")
        sys.exit(1)
    except Exception as e:
        click.echo(f"  ✗ Setup failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Run simulation or just generate config
    if dry_run:
        click.echo(f"\n[4/4] Dry run - generating configuration files...")

        if force_field == 'calvados':
            # 使用 CGSimulator 的目录逻辑生成配置
            try:
                from ...src.calvados_wrapper import CalvadosWrapper

                # 准备 raw 目录（与 run_calvados 相同的目录结构）
                raw_dir = os.path.join(final_output_dir, 'raw')

                # 写入配置文件到 raw 目录
                wrapper = CalvadosWrapper(config)
                wrapper._write_to_dir(raw_dir)

                click.echo(f"  ✓ Configuration files generated")
                click.echo(f"  Output: {final_output_dir}")
                click.echo(f"  - config.yaml")
                click.echo(f"  - components.yaml")
                click.echo(f"\n  To run simulation, remove --dry-run flag.")
            except Exception as e:
                click.echo(f"  ✗ Failed to generate config: {e}")
                if verbose:
                    import traceback
                    traceback.print_exc()
                sys.exit(1)
        else:
            # 其他力场：先说明需要用 CALVADOS 生成结构
            click.echo(f"  Force field: {force_field}")
            click.echo(f"\n  For non-CALVADOS force fields:")
            click.echo(f"    1. CALVADOS will first generate initial structure")
            click.echo(f"    2. Then the selected force field will run simulation")
            click.echo(f"\n  To run full simulation, remove --dry-run flag.")
    else:
        click.echo(f"\n[4/4] Running {force_field.upper()} simulation...")

        # Call appropriate runner based on force field
        # Map CLI force field name to internal runner method name
        runner_name = FORCE_FIELD_TO_RUNNER.get(force_field, force_field)
        runner_method = f'run_{runner_name}'
        if not hasattr(sim, runner_method):
            click.echo(f"  ✗ Force field '{force_field}' not yet implemented.")
            sys.exit(1)

        try:
            # Call runner method with gpu_id and continue_from parameters
            # For mpipi_recharged, defaults to use_gmx_insert=True, gmx_radius=0.35
            # continue_from is only supported for calvados
            if force_field == 'calvados':
                result = getattr(sim, runner_method)(gpu_id=gpu_id, continue_from=continue_from)
            else:
                result = getattr(sim, runner_method)(gpu_id=gpu_id)

            # Use the result output directory
            actual_output = result.output_dir if result.output_dir else final_output_dir

            if result.success:
                click.echo(f"  ✓ Simulation completed")
                click.echo(f"  Output: {actual_output}")

                if result.trajectory and os.path.exists(result.trajectory):
                    size_mb = os.path.getsize(result.trajectory) / 1024 / 1024
                    click.echo(f"  Trajectory: {os.path.basename(result.trajectory)} ({size_mb:.1f} MB)")
            else:
                click.echo(f"  ✗ Simulation failed:")
                for error in result.errors:
                    click.echo(f"    - {error}")
                sys.exit(1)

        except Exception as e:
            click.echo(f"  ✗ Simulation error: {e}")
            if verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
    
    click.echo()
