#!/usr/bin/env python3
"""
PACE-Opt Command

Optimizes backmapped structure with PACE force field.
"""

import sys

import click


@click.command("pace-opt", context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--input",
    "-i",
    type=click.Path(exists=True),
    required=True,
    help="Input: backmap output directory (adapter backmap) or PDB file (user provided)",
)
@click.option(
    "--input-file",
    "-f",
    type=click.Path(exists=True),
    help="Configuration YAML (same as CG stage, same flag as adapter cg)",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    help="Output directory (default: {system_name}_pace_opt)",
)
@click.option(
    "--skip-gromacs",
    is_flag=True,
    default=False,
    help="Skip GROMACS optimization step",
)
@click.option(
    "--device",
    "-d",
    type=click.Choice(["cpu", "cuda", "opencl"]),
    default="cuda",
    help="Device for OpenMM (default: cuda)",
)
@click.option(
    "--box-resize",
    type=str,
    help='Box resize dimensions as comma-separated values (e.g., "23.0,17.0,47.0")',
)
@click.option(
    "-l",
    "--level",
    type=click.Choice(["high", "medium", "low"]),
    default="medium",
    help="Optimization level: high (7 steps), medium (5 steps), low (3 steps). Default: medium",
)
@click.option(
    "--forcefield-type",
    "-ff-type",
    type=click.Choice(["pace-new", "pace-asm"]),
    default="pace-new",
    help="PACE force field type: pace-new (standard) or pace-asm (condensate variant). Default: pace-new",
)
def pace_opt_command(
    input, input_file, output, skip_gromacs, device, box_resize, level, forcefield_type
):
    """Optimize backmapped structure with PACE force field.

    Uses multi-step optimization: Gaussian repulsion → Softcore (multiple steps) → Standard force field.

    Force field types:
        pace-new: Standard PACE force field (default) - for general protein structures
        pace-asm: Alternative PACE variant for condensate structures (TODO: not yet implemented)

    Optimization levels:
        high:   7 steps (0.2 → 0.15 → 0.1 → 0.075 → 0.05 → 0.025 → 0.01)
        medium: 5 steps (0.2 → 0.1 → 0.05 → 0.025 → 0.01) - default
        low:    3 steps (0.15 → 0.05 → 0.025)

    Examples:
        # From adapter backmap output (using default PACE-NEW)
        adapter pace-opt -i TDP43_backmap -f config.yaml

        # From user provided PDB
        adapter pace-opt -i my_structure.pdb -f config.yaml

        # With box resize
        adapter pace-opt -i TDP43_backmap -f config.yaml --box-resize "23.0,17.0,47.0"

        # Use PACE-ASM for condensate structures (not yet implemented)
        adapter pace-opt -i condensate_backmap -f config.yaml --forcefield-type pace-asm

        # Use faster optimization (medium mode)
        adapter pace-opt -i TDP43_backmap -f config.yaml --level medium

        # Skip GROMACS optimization
        adapter pace-opt -i TDP43_backmap -f config.yaml --skip-gromacs

        # Use CPU instead of CUDA
        adapter pace-opt -i TDP43_backmap -f config.yaml -d cpu
    """
    from ...src import CGSimulationConfig
    from ...src.pace_opt import (
        PaceOptConfig,
        PaceOptSimulator,
    )

    click.echo(f"\n{'=' * 60}")
    click.echo("PACE Force Field Optimization")
    click.echo(f"{'=' * 60}")

    click.echo(f"\n  Input: {input}")
    if input_file:
        click.echo(f"  Config: {input_file}")
    if output:
        click.echo(f"  Output: {output}")
    click.echo(f"  Device: {device.upper()}")
    click.echo(f"  Force field type: {forcefield_type}")
    if box_resize:
        click.echo(f"  Box resize: {box_resize}")
    if skip_gromacs:
        click.echo("  GROMACS optimization: Skipped")

    # Load configuration
    config = None
    if input_file:
        config = CGSimulationConfig.from_yaml(input_file)

    # Parse box resize if provided
    box_resize_dims = None
    if box_resize:
        try:
            box_resize_dims = [float(x.strip()) for x in box_resize.split(",")]
            if len(box_resize_dims) != 3:
                raise ValueError("Box resize must have 3 dimensions")
        except ValueError as e:
            click.echo(f"  ✗ Invalid box resize format: {e}")
            sys.exit(1)

    # Create pace_opt_config
    pace_opt_config = PaceOptConfig(forcefield_type=forcefield_type)
    pace_opt_config.box_resize_enabled = box_resize_dims is not None
    pace_opt_config.box_resize_dimensions = box_resize_dims
    pace_opt_config.gromacs_enabled = not skip_gromacs
    pace_opt_config.platform = device.upper()
    pace_opt_config.set_optimization_mode(level)

    # Create simulator
    simulator = PaceOptSimulator(config=config, pace_opt_config=pace_opt_config)

    # Run optimization
    try:
        result = simulator.run(input_path=input, config_path=input_file, output_dir=output)

        if result.success:
            click.echo("\n  ✓ PACE optimization completed")
            click.echo(f"  Input PDB: {result.input_pdb}")
            click.echo(f"  Output PDB: {result.output_pdb}")
            if result.intermediate_files:
                click.echo("  Intermediate files:")
                for f in result.intermediate_files:
                    click.echo(f"    - {f}")
            click.echo("\n  ✓ Success!")
        else:
            click.echo("  ✗ PACE optimization failed:")
            for error in result.errors:
                click.echo(f"    - {error}")
            sys.exit(1)

    except Exception as e:
        click.echo(f"  ✗ Error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)

    click.echo()
