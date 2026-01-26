#!/usr/bin/env python3
"""
Minimization Worker Script

Runs multi-step energy minimization using OpenMM in a subprocess.
Outputs only generated PDB files for real-time feedback.
"""

import os
import sys
import gc
import shutil
from pathlib import Path

try:
    import openmm.unit as unit
    import openmm as mm
    from openmm.app import GromacsGroFile, Simulation, PDBFile
    from openmm import LangevinIntegrator, Platform
except ImportError:
    import simtk.unit as unit
    import simtk.openmm as mm
    from simtk.openmm.app import GromacsGroFile, Simulation, PDBFile
    from simtk.openmm import LangevinIntegrator, Platform

from openmm import openmm as omm

# Import from top_to_softcore_system
sys.path.insert(0, str(Path(__file__).parent))
from top_to_softcore_system import (
    GromacsTopFileWithSoftcore,
    NONBONDED_GAUSSIAN,
    NONBONDED_SOFTCORE,
    NONBONDED_STANDARD,
    IMPLICIT_GBN2,
    IMPLICIT_OBC2
)

from device_utils import resolve_openmm_platform

# Import force field utilities from openmm.app
from openmm.app import forcefield as ff


def get_lambda_values(level: str):
    """Get lambda values for softcore optimization.
    
    Args:
        level: Optimization level (high, medium, low)
    
    Returns:
        List of lambda values (direct, not 1.0 - alpha)
    """
    configs = {
        'high': [0.65, 0.75, 0.85, 0.95],
        'medium': [0.75, 0.85, 0.95],
        'low': [0.85, 0.95],
    }
    level = level.lower()
    if level not in configs:
        raise ValueError(f"Unknown optimization level: {level}. Valid: {list(configs.keys())}")
    return configs[level]


def run_minimization(
    input_gro: str,
    input_top: str,
    output_dir: str,
    device: str = "cuda",
    gpu_id: int = 0,
    max_iterations: int = 5000,
    optimization_level: str = "medium",
    ff_type: str = "amber",
    ff_name: str = "amber99sb-ildn",
    gb_model: str = "GBn2",
    salt_conc: float = 0.15,
    cutoff: float = 1.1,
    tolerance: float = 100.0
):
    """
    Run multi-step energy minimization.
    
    Args:
        input_gro: Input GRO file path
        input_top: Input TOP file path
        output_dir: Output directory path
        device: Compute device (cuda, cpu, opencl)
        gpu_id: GPU device index
        max_iterations: Maximum minimization iterations per step
        optimization_level: high/medium/low for softcore steps
        ff_type: Force field type for OpenMM (amber or charmm)
        ff_name: Full force field name for reference
        gb_model: GB model for implicit solvent (GBn2 or OBC2)
        salt_conc: Salt concentration in M
        cutoff: Nonbonded cutoff in nm
        tolerance: Minimization tolerance
    """
    input_gro_path = Path(input_gro).resolve()
    input_top_path = Path(input_top).resolve()
    output_path = Path(output_dir).resolve()
    
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Copy input files to output directory
    conf_gro = output_path / 'conf.gro'
    minimize_top = output_path / 'topol.top'
    
    if str(input_gro_path) != str(conf_gro):
        shutil.copy2(str(input_gro_path), str(conf_gro))
    if str(input_top_path) != str(minimize_top):
        shutil.copy2(str(input_top_path), str(minimize_top))
    
    # Change to output directory
    original_cwd = os.getcwd()
    os.chdir(output_path)
    
    try:
        # Load structure
        conf = GromacsGroFile('conf.gro')
        gro_position = conf.getPositions()
        box_vectors = conf.getPeriodicBoxVectors()
        
        # Create topology
        top = GromacsTopFileWithSoftcore(
            'topol.top',
            periodicBoxVectors=box_vectors,
            forcefield_type=ff_type.upper()
        )
        
        # Setup platform
        platform_name, properties = resolve_openmm_platform(device, gpu_id, precision="mixed")
        platform = Platform.getPlatformByName(platform_name)
        
        # Store current positions
        current_positions = gro_position
        
        # ===== Step 1: Gaussian Repulsion =====
        system_gaussian = top.createSystem(
            nonbondedCutoff=cutoff * unit.nanometer,
            nonbondedMethod=ff.CutoffPeriodic,
            nonbonded_type=NONBONDED_GAUSSIAN,
            add_implicit_solvent=False
        )
        
        integrator = LangevinIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picosecond,
            0.002 * unit.picosecond
        )
        integrator.setRandomNumberSeed(0)
        
        simulation = Simulation(top.topology, system_gaussian, integrator, platform, properties)
        simulation.context.setPositions(current_positions)
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        
        simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=50.0)
        
        state_gaussian = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
        output_file = 'step1_gaussian.pdb'
        PDBFile.writeFile(top.topology, state_gaussian.getPositions(), open(output_file, 'w'))
        print(f"  Generated {output_file}")
        
        current_positions = state_gaussian.getPositions()
        
        # Clean up
        del system_gaussian, simulation, state_gaussian
        gc.collect()
        
        # ===== Step 2: Softcore Optimization with Lambda =====
        lambda_values = get_lambda_values(optimization_level)
        
        for step_num, lambda_val in enumerate(lambda_values, 1):
            system_softcore = top.createSystem(
                nonbondedCutoff=cutoff * unit.nanometer,
                nonbondedMethod=ff.CutoffPeriodic,
                nonbonded_type=NONBONDED_SOFTCORE,
                add_implicit_solvent=False,
                soft_lambda=lambda_val
            )
            
            integrator = LangevinIntegrator(
                300 * unit.kelvin,
                1.0 / unit.picosecond,
                0.002 * unit.picosecond
            )
            integrator.setRandomNumberSeed(0)
            
            simulation = Simulation(top.topology, system_softcore, integrator, platform, properties)
            simulation.context.setPositions(current_positions)
            simulation.context.setPeriodicBoxVectors(*box_vectors)
            
            simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=tolerance)
            
            state_softcore = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
            current_positions = state_softcore.getPositions()
            
            output_file = f'step2_softcore_{step_num}.pdb'
            PDBFile.writeFile(top.topology, current_positions, open(output_file, 'w'))
            print(f"  Generated {output_file} (lambda={lambda_val})")
            
        # Clean up
        del system_softcore, simulation, state_softcore
        gc.collect()

        # ===== Step 3: Final Minimization =====
        # For CHARMM: use softcore with lambda=0.99999 (bypasses GB issues)
        # For AMBER: use standard with implicit solvent
        if ff_type.upper() == 'CHARMM':
            # CHARMM: use softcore with very high lambda
            final_lambda = 0.99999
            print(f"  Using softcore final minimize (lambda={final_lambda}) for CHARMM")

            system_final = top.createSystem(
                nonbondedCutoff=cutoff * unit.nanometer,
                nonbondedMethod=ff.CutoffPeriodic,
                nonbonded_type=NONBONDED_SOFTCORE,
                add_implicit_solvent=False,
                soft_lambda=final_lambda
            )
        else:
            # AMBER: use standard with implicit solvent
            gb_constant = IMPLICIT_GBN2 if gb_model.upper() == 'GBN2' else IMPLICIT_OBC2

            system_final = top.createSystem(
                nonbondedCutoff=cutoff * unit.nanometer,
                nonbondedMethod=ff.CutoffPeriodic,
                nonbonded_type=NONBONDED_STANDARD,
                add_implicit_solvent=True,
                gb_model=gb_constant,
                salt_conc=salt_conc
            )

        integrator = LangevinIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picosecond,
            0.002 * unit.picosecond
        )
        integrator.setRandomNumberSeed(0)

        simulation = Simulation(top.topology, system_final, integrator, platform, properties)
        simulation.context.setPositions(current_positions)
        simulation.context.setPeriodicBoxVectors(*box_vectors)

        simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=tolerance)

        state_final = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)

        # Save final output
        PDBFile.writeFile(top.topology, state_final.getPositions(), open('minimize_final.pdb', 'w'))
        print(f"  Generated minimize_final.pdb")

        # Clean up
        del system_final, simulation, state_final
        gc.collect()
        
        # Calculate total optimization steps
        total_steps = 1 + len(lambda_values) + 1  # Step 1 + Step 2 (N steps) + Step 3
        step_info = f"{optimization_level.capitalize()} ({total_steps} steps)"
        
        print(f"  Optimization: {step_info}")
        print(f"\n  Done!")
        
        # Return optimization information
        return {
            'step_info': step_info,
            'total_steps': total_steps,
        }
        
    finally:
        os.chdir(original_cwd)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Minimize Worker Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python minimize_worker.py -i conf.gro -t topol.top -o output_dir
    python minimize_worker.py --input-gro conf.gro --input-top topol.top --output my_output --device cuda
        """
    )
    
    parser.add_argument('-i', '--input-gro', required=True, 
                        help='Input GRO file')
    parser.add_argument('-t', '--input-top', required=True, 
                        help='Input topology file')
    parser.add_argument('-o', '--output', required=True, 
                        help='Output directory')
    parser.add_argument('-d', '--device', default='cuda', choices=['cuda', 'cpu', 'opencl'],
                        help='Compute device (default: cuda)')
    parser.add_argument('-g', '--gpu-id', type=int, default=0,
                        help='GPU device index (default: 0, only used when device is cuda or opencl)')
    parser.add_argument('--iter', type=int, default=5000,
                        help='Maximum iterations per step (default: 5000)')
    parser.add_argument('--tolerance', type=float, default=100.0,
                        help='Minimization tolerance (default: 100.0)')
    parser.add_argument('-l', '--level', default='medium', choices=['high', 'medium', 'low'],
                        help='Optimization level (default: medium)')
    parser.add_argument('--ff-type', default='amber', choices=['amber', 'charmm'],
                        help='Force field type (default: amber)')
    parser.add_argument('--ff-name', default='amber99sb-ildn',
                        help='Force field name for reference (default: amber99sb-ildn)')
    parser.add_argument('--gb-model', default='GBn2', choices=['GBn2', 'OBC2'],
                        help='GB model for implicit solvent (default: GBN2)')
    parser.add_argument('--salt-conc', type=float, default=0.15,
                        help='Salt concentration in M (default: 0.15)')
    parser.add_argument('--cutoff', type=float, default=2.0,
                        help='Nonbonded cutoff in nm (default: 2.0)')
    
    args = parser.parse_args()
    
    run_minimization(
        input_gro=args.input_gro,
        input_top=args.input_top,
        output_dir=args.output,
        device=args.device,
        gpu_id=args.gpu_id,
        max_iterations=args.iter,
        optimization_level=args.level,
        ff_type=args.ff_type,
        ff_name=args.ff_name,
        gb_model=args.gb_model,
        salt_conc=args.salt_conc,
        cutoff=args.cutoff,
        tolerance=args.tolerance
    )
