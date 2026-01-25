#!/usr/bin/env python3
"""
Minimize Command

Energy minimization using AMBER/CHARMM force fields.
"""

import sys

import click


@click.command('minimize', context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--input', '-i',
    type=click.Path(exists=True),
    required=True,
    help='Input: backmap output directory (adapter backmap) or PDB file (user provided)',
)
@click.option(
    '--input-file', '-f',
    type=click.Path(exists=True),
    required=True,
    help='Configuration YAML (same as CG stage, same flag as adapter cg)',
)
@click.option(
    '--output', '-o',
    type=click.Path(),
    help='Output directory (default: {system_name}_minimize)',
)
@click.option(
    '--force-field', '-ff',
    type=str,
    default='1-a99SBdisp',
    callback=None,  # Will be set in function
    help='Force field for minimization. Use number (1-9) or full name.\n'
         '  1-a99SBdisp (recommended)\n'
         '  2-amber03wsc (recommended)\n'
         '  3-amber99sbws-stqp (recommended)\n'
         '  4-amber99sbws-stq\n'
         '  5-des-amber\n'
         '  6-des-amber-sf1.0\n'
         '  7-amber99sb-ildn\n'
         '  8-amber14sb\n'
         '  9-charmm36m\n'
         'Default: 1 (a99SBdisp)',
)
@click.option(
    '--device', '-d',
    type=click.Choice(['cpu', 'cuda', 'opencl']),
    default='cuda',
    help='Device for OpenMM (default: cuda)',
)
@click.option(
    '-l', '--level',
    type=click.Choice(['high', 'medium', 'low']),
    default='medium',
    help='Optimization level: high (7 steps), medium (5 steps), low (3 steps). Default: medium',
)
@click.option(
    '--tolerance',
    type=float,
    default=100.0,
    help='Minimization tolerance in kJ/(mol·nm). Default: 100.0',
)
@click.option(
    '--iter',
    type=int,
    default=5000,
    help='Maximum minimization iterations per step. Default: 5000',
)
@click.option(
    '--salt-conc',
    type=float,
    default=0.15,
    help='Salt/ion concentration in M (for both implicit and explicit solvent). Default: 0.15',
)
@click.option(
    '--cutoff',
    type=float,
    default=2.0,
    help='Nonbonded cutoff in nm. Default: 2.0',
)
@click.option(
    '--gpu-id', '-g',
    type=int,
    default=0,
    help='GPU device index (only used when device is cuda or opencl, default: 0)',
)
@click.option(
    '--solvate',
    is_flag=True,
    default=False,
    help='Enable explicit solvation (adds water box and ions)',
)
@click.option(
    '--no-disulfide',
    is_flag=True,
    default=False,
    help='Disable disulfide detection in pdb2gmx (-ss). Distance-based SS '
         'assignment can give different bond counts across identical monomers '
         'in a condensate, which may cause build failures (e.g., atom-count '
         'mismatch). Use this flag if you hit that issue.',
)
@click.option(
    '--his-type',
    type=click.Choice(['0', '1']),
    default=None,
    help='Use pdb2gmx -his to set all histidines: 0=HID, 1=HIE. '
         'The CLI pre-fills (total_nmol * 30) entries to avoid interaction.',
)
def minimize_command(input, input_file, output, force_field, device, gpu_id, level, tolerance, iter, salt_conc, cutoff, solvate, no_disulfide, his_type):
    """Energy minimization using AMBER/CHARMM force fields.

    Workflow:
    1. Implicit solvent optimization (GBn2 model) → minimize_final.pdb
    2. (Optional) If --solvate enabled:
       - Add water box and ions to minimize_final.pdb
       - Output: minimize_final_solvated.gro + topol.top

    Uses gromacs pdb2gmx for topology generation and multi-step OpenMM
    minimization (Gaussian -> Softcore -> Standard) with GBn2 implicit solvent.

    Available force fields:
        amber99sb-ildn:    AMBER99SB-ILDN with tip3p (default)
        amber99sbws-STQp:  AMBER99SB-WS (STQp) with tip4p2005s
        amber99sbws-stq:   AMBER99SB-WS (stq) with tip4p2005s
        amber03wsc:        AMBER03wsc with tip4p2005s
        amber14sb:         AMBER14SB with tip3p
        a99SBdisp:         a99SB-disp with custom water
        des-amber:         DES-AMBER with tip4pd
        des-amber-SF1.0:   DES-AMBER SF1.0 with tip4pd
        charmm36:          CHARMM36-jul2021 with tip3p

    Optimization levels:
        high:   4 steps (lambda: 0.65 → 0.75 → 0.85 → 0.95)
        medium: 3 steps (lambda: 0.75 → 0.85 → 0.95) - default
        low:    2 steps (lambda: 0.85 → 0.95)

    Examples:
        # From adapter backmap output (using a99SBdisp with GBn2)
        adapter minimize -i TDP43_backmap -f config.yaml

        # From user provided PDB
        adapter minimize -i my_structure.pdb -f config.yaml

        # With CHARMM36 force field
        adapter minimize -i TDP43_backmap -f config.yaml --force-field 9

        # Use faster optimization (low mode)
        adapter minimize -i TDP43_backmap -f config.yaml --level low

        # Use CPU instead of CUDA
        adapter minimize -i TDP43_backmap -f config.yaml -d cpu
        
        # Use specific GPU (GPU 1)
        adapter minimize -i TDP43_backmap -f config.yaml --gpu-id 1
        
        # With explicit solvation (water box + ions)
        adapter minimize -i TDP43_backmap -f config.yaml --solvate
        
        # With explicit solvation and custom salt concentration
        adapter minimize -i TDP43_backmap -f config.yaml --solvate --salt-conc 0.2
    """
    from ...src.minimize import MinimizeSimulator, MinimizeConfig
    from ..shared import validate_minimize_force_field, REGISTRY

    # Set the callback for force field validation
    for param in minimize_command.params:
        if param.name == 'force_field':
            param.callback = validate_minimize_force_field
            break

    click.echo(f"\n{'=' * 60}")
    click.echo(f"Energy Minimization (AMBER/CHARMM)")
    click.echo(f"{'=' * 60}")

    click.echo(f"\n  Input: {input}")
    click.echo(f"  Config: {input_file}")
    if output:
        click.echo(f"  Output: {output}")
    click.echo(f"  Force field: {force_field}")
    ff_info = REGISTRY.get_force_field(force_field)
    if ff_info:
        click.echo(f"    Family: {ff_info.family}")
        click.echo(f"    Water model: {ff_info.water_model}")
        click.echo(f"    GBSA mapping: {ff_info.family}")  # Show family (AMBER/CHARMM) instead of technical mapping name
    click.echo(f"  GB model: GBn2 (implicit solvent)")
    click.echo(f"  Device: {device.upper()}" + (f" (GPU {gpu_id})" if device.upper() != 'CPU' else ""))
    click.echo(f"  Optimization level: {level}")
    click.echo(f"  Tolerance: {tolerance}")
    click.echo(f"  Iterations: {iter}")
    click.echo(f"  Salt conc: {salt_conc} M")
    click.echo(f"  Cutoff: {cutoff} nm")
    if solvate:
        click.echo(f"  Solvate: Enabled")
    if no_disulfide:
        click.echo(f"  Disulfide detection: Disabled (pdb2gmx -ss)")
    if his_type is not None:
        his_type_int = int(his_type)
        click.echo(f"  Histidine type: {'HID' if his_type_int == 0 else 'HIE'} (pdb2gmx -his)")

    # Create minimize_config (GB model fixed to GBn2)
    minimize_config = MinimizeConfig(
        forcefield_type=force_field,
        gb_model='GBn2',  # Fixed to GBn2
        platform=device.upper(),
        gpu_id=gpu_id,
        tolerance=tolerance,
        max_iterations=iter,
        salt_conc=salt_conc,
        nonbonded_cutoff=cutoff,
        solvate_enabled=solvate,
        ion_concentration=salt_conc,
        disable_disulfide=no_disulfide,
        his_type=int(his_type) if his_type is not None else None
    )
    minimize_config.set_optimization_mode(level)

    # Create simulator with components from YAML
    simulator = MinimizeSimulator.from_yaml(input_file, minimize_config=minimize_config)

    # Run minimization
    try:
        result = simulator.run(input_pdb=input, output_dir=output)

        if result.success:
            click.echo(f"\n  Minimization completed")
            click.echo(f"  Input PDB: {result.input_pdb}")
            click.echo(f"  Output PDB: {result.output_pdb}")
            if result.step_info:
                click.echo(f"  Optimization: {result.step_info}")
            click.echo(f"\n  Success!")
        else:
            click.echo(f"  Minimization failed:")
            for error in result.errors:
                click.echo(f"    - {error}")
            sys.exit(1)

    except Exception as e:
        click.echo(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    click.echo()
