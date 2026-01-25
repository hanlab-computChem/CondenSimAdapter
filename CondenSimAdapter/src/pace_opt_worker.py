#!/usr/bin/env python3
"""
PACE Optimization Worker Script

This script runs the complete PACE optimization workflow using CUDA platform.
It uses pre-generated topology files to avoid GROMACS issues.
Supports retry logic for each optimization step.
"""

import sys
import os
import gc
from pathlib import Path
import shutil

# Add paths
sys.path.insert(0, '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/CondenSimAdapter/src')
sys.path.insert(0, '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square')

import PACE2openmm
from openmm.app import GromacsGroFile, Simulation, PDBFile
from openmm import Platform, LangevinIntegrator
import openmm.unit as unit


# Optimization level to softcore configs mapping
def get_softcore_configs(level: str):
    """Get softcore configurations based on optimization level.

    Args:
        level: One of 'high', 'medium', 'low'

    Returns:
        List of (alpha_vdw, alpha_coul, step_num) tuples
    """
    configs = {
        'high': [
            (0.2, 0.2, 1),
            (0.15, 0.15, 2),
            (0.1, 0.1, 3),
            (0.075, 0.075, 4),
            (0.05, 0.05, 5),
            (0.025, 0.025, 6),
            (0.01, 0.01, 7),
        ],
        'medium': [
            (0.2, 0.2, 1),
            (0.1, 0.1, 2),
            (0.05, 0.05, 3),
            (0.025, 0.025, 4),
            (0.01, 0.01, 5),
        ],
        'low': [
            (0.15, 0.15, 1),
            (0.05, 0.05, 2),
            (0.025, 0.025, 3),
        ],
    }
    level = level.lower()
    if level not in configs:
        raise ValueError(f"Unknown optimization level: {level}. Valid: {list(configs.keys())}")
    return configs[level]


def run_pace_optimization(
    input_gro: str,
    input_top: str,
    output_dir: str,
    device: str = "cuda",
    max_iterations: int = 5000,
    optimization_level: str = "high"
):
    """
    Run the complete PACE optimization workflow.
    """
    input_gro_path = Path(input_gro).resolve()
    input_top_path = Path(input_top).resolve()
    output_path = Path(output_dir).resolve()
    
    print(f"\n{'='*60}")
    print(f"PACE Force Field Optimization")
    print(f"{'='*60}")
    print(f"\n  Input GRO: {input_gro_path}")
    print(f"  Input TOP: {input_top_path}")
    print(f"  Output: {output_path}")
    print(f"  Device: {device.upper()}")
    print(f"  Optimization level: {optimization_level}")
    print(f"  Max iterations: {max_iterations}")
    
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Copy input files
    print(f"\n[1/4] Preparing files...")
    conf_gro = output_path / 'conf.gro'
    pace_top = output_path / 'PACE.top'
    
    if str(input_gro_path) != str(conf_gro):
        shutil.copy2(str(input_gro_path), str(conf_gro))
    else:
        print(f"  Input GRO is already in output directory, using as-is")
    
    if str(input_top_path) != str(pace_top):
        shutil.copy2(str(input_top_path), str(pace_top))
    else:
        print(f"  Input TOP is already in output directory, using as-is")
    
    # Change to output directory
    original_cwd = os.getcwd()
    os.chdir(output_path)
    
    try:
        # Load structure
        print(f"[2/4] Loading structure...")
        conf = GromacsGroFile('conf.gro')
        gro_position = conf.getPositions()
        box_vectors = conf.getPeriodicBoxVectors()
        print(f"  Loaded {len(gro_position)} atoms")
        
        # Create PACE topology
        gmxlib = '/mnt/hdd1/tianxj_out/gmx2023.2.cuda/share/gromacs/top'
        top = PACE2openmm.PACETopFile(
            'PACE.top',
            periodicBoxVectors=box_vectors,
            defines={},
            epsilon_r=15.0,
            includeDir=gmxlib
        )
        print(f"  PACE topology loaded")
        
        # Setup platform
        platform = Platform.getPlatformByName(device.upper())
        if device.upper() == 'CUDA':
            properties = {'Precision': 'double'}
        else:
            properties = {}
        print(f"  Platform: {platform.getName()}")
        
        print(f"\n[3/4] Running OpenMM optimization...")
        
        # ===== Step 1: Gaussian Optimization =====
        print(f"\n  Step 1: Gaussian repulsion optimization...")
        
        max_retries = 5
        for attempt in range(1, max_retries + 1):
            try:
                system_gaussian = top.create_system(
                    nonbonded_cutoff=1.1 * unit.nanometer,
                    add_nonbonded_force=True,
                    nonbonded_type='gaussian',
                    add_ring_constraints=False  # Ring constraints disabled for Gaussian step
                )
                
                integrator = LangevinIntegrator(
                    325 * unit.kelvin,
                    10.0 / unit.picosecond,
                    0.004 * unit.picosecond
                )
                integrator.setRandomNumberSeed(0)
                
                simulation = Simulation(top.topology, system_gaussian, integrator, platform, properties)
                simulation.context.setPositions(gro_position)
                
                state_before = simulation.context.getState(getEnergy=True)
                print(f"    Energy before: {state_before.getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole):.2f} kJ/mol")
                
                simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=50.0)
                
                state_gaussian = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
                print(f"    Energy after: {state_gaussian.getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole):.2f} kJ/mol")
                
                PDBFile.writeFile(top.topology, state_gaussian.getPositions(), open('opti_gaussian.pdb', 'w'))
                print(f"    Saved: opti_gaussian.pdb")
                
                current_positions = state_gaussian.getPositions()
                break
                
            except Exception as e:
                error_str = str(e)
                if any(x in error_str.lower() for x in ['nan', 'segmentation', 'core dumped', 'memory']):
                    print(f"    ⚠️ Attempt {attempt} failed with fatal error: {error_str}")
                    raise RuntimeError(f"Gaussian optimization failed with fatal error: {error_str}")
                else:
                    print(f"    ⚠️ Attempt {attempt} failed: {error_str}")
                    if attempt < max_retries:
                        print(f"    Retrying...")
                        gc.collect()
                    else:
                        raise RuntimeError(f"Gaussian optimization failed after {max_retries} attempts: {error_str}")
        
        # ===== Step 2: Softcore Optimization =====
        softcore_configs = get_softcore_configs(optimization_level)
        total_softcore_steps = len(softcore_configs)

        for alpha, alpha_coul, step_num in softcore_configs:
            # Use same alpha value for exception force as VDW alpha
            alpha_exc = alpha
            # Enable ring constraints for step 2.1 onwards (all softcore steps)
            add_ring_constraints = True
            print(f"\n  Step 2.{step_num}: Softcore optimization (alpha={alpha})...")
            
            for attempt in range(1, max_retries + 1):
                try:
                    system_softcore = top.create_system(
                        nonbonded_cutoff=1.1 * unit.nanometer,
                        add_nonbonded_force=True,
                        nonbonded_type='softcore',
                        add_ring_constraints=add_ring_constraints
                    )
                    
                    integrator = LangevinIntegrator(
                        325 * unit.kelvin,
                        10.0 / unit.picosecond,
                        0.004 * unit.picosecond
                    )
                    integrator.setRandomNumberSeed(0)
                    
                    simulation = Simulation(top.topology, system_softcore, integrator, platform, properties)
                    simulation.context.setPositions(current_positions)
                    
                    simulation.context.setParameter('soft_alpha', alpha)
                    simulation.context.setParameter('soft_alpha_coul', alpha_coul)
                    simulation.context.setParameter('soft_alpha_exc', alpha_exc)
                    
                    state_before = simulation.context.getState(getEnergy=True)
                    print(f"    Energy before: {state_before.getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole):.2f} kJ/mol")
                    
                    # Use tolerance=200 for step 2.1 onwards
                    simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=200)
                    
                    state_softcore = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
                    print(f"    Energy after: {state_softcore.getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole):.2f} kJ/mol")
                    
                    current_positions = state_softcore.getPositions()
                    
                    output_file = f'opti_softcore_{step_num}.pdb'
                    PDBFile.writeFile(top.topology, current_positions, open(output_file, 'w'))
                    print(f"    Saved: {output_file}")
                    
                    break
                    
                except Exception as e:
                    error_str = str(e)
                    if any(x in error_str.lower() for x in ['nan', 'segmentation', 'core dumped', 'memory']):
                        print(f"    ⚠️ Attempt {attempt} failed with fatal error: {error_str}")
                        raise RuntimeError(f"Softcore step {step_num} failed with fatal error: {error_str}")
                    else:
                        print(f"    ⚠️ Attempt {attempt} failed: {error_str}")
                        if attempt < max_retries:
                            print(f"    Retrying...")
                            gc.collect()
                        else:
                            raise RuntimeError(f"Softcore step {step_num} failed after {max_retries} attempts: {error_str}")
        
        # ===== Step 3: Standard Optimization =====
        print(f"\n  Step 3: Standard force field optimization...")
        
        for attempt in range(1, max_retries + 1):
            try:
                system_standard = top.create_system(
                    nonbonded_cutoff=1.1 * unit.nanometer,
                    add_nonbonded_force=True,
                    nonbonded_type='standard',
                    add_ring_constraints=True  # Enable ring constraints for full potential
                )
                
                integrator = LangevinIntegrator(
                    325 * unit.kelvin,
                    10.0 / unit.picosecond,
                    0.004 * unit.picosecond
                )
                integrator.setRandomNumberSeed(0)
                
                simulation = Simulation(top.topology, system_standard, integrator, platform, properties)
                simulation.context.setPositions(current_positions)
                
                state_before = simulation.context.getState(getEnergy=True)
                print(f"    Energy before: {state_before.getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole):.2f} kJ/mol")
                
                # Step 3: Use tolerance=200
                simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=200)
                
                state_final = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
                print(f"    Energy after: {state_final.getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole):.2f} kJ/mol")
                
                PDBFile.writeFile(top.topology, state_final.getPositions(), open('opti_final.pdb', 'w'))
                print(f"    Saved: opti_final.pdb")
                
                break
                
            except Exception as e:
                error_str = str(e)
                if any(x in error_str.lower() for x in ['nan', 'segmentation', 'core dumped', 'memory']):
                    print(f"    ⚠️ Attempt {attempt} failed with fatal error: {error_str}")
                    raise RuntimeError(f"Standard optimization failed with fatal error: {error_str}")
                else:
                    print(f"    ⚠️ Attempt {attempt} failed: {error_str}")
                    if attempt < max_retries:
                        print(f"    Retrying...")
                        gc.collect()
                    else:
                        raise RuntimeError(f"Standard optimization failed after {max_retries} attempts: {error_str}")
        
        print(f"\n[4/4] Complete!")
        print(f"\n  Output files in: {output_path}")
        print(f"  - opti_gaussian.pdb")
        for alpha, _, step_num in softcore_configs:
            print(f"  - opti_softcore_{step_num}.pdb (alpha={alpha})")
        print(f"  - opti_final.pdb (full potential)")
        
        print(f"\n{'='*60}")
        print(f"PACE optimization completed successfully!")
        print(f"{'='*60}\n")
        
    finally:
        os.chdir(original_cwd)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='PACE Optimization Worker Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python pace_opt_worker.py -i processed.gro -t PACE.top -o output_dir
    python pace_opt_worker.py --input-gro conf.gro --input-top top.top --output my_output --device cuda
        """
    )
    
    parser.add_argument('-i', '--input-gro', required=True, 
                        help='Input GRO file')
    parser.add_argument('-t', '--input-top', required=True, 
                        help='Input PACE topology file')
    parser.add_argument('-o', '--output', required=True, 
                        help='Output directory')
    parser.add_argument('-d', '--device', default='cuda', choices=['cuda', 'cpu', 'opencl'],
                        help='Compute device (default: cuda)')
    parser.add_argument('--iter', type=int, default=5000,
                        help='Maximum iterations (default: 5000)')
    parser.add_argument('-l', '--level', default='high', choices=['high', 'medium', 'low'],
                        help='Optimization level: high (7 steps), medium (5 steps), low (3 steps). Default: high')

    args = parser.parse_args()

    run_pace_optimization(
        input_gro=args.input_gro,
        input_top=args.input_top,
        output_dir=args.output,
        device=args.device,
        max_iterations=args.iter,
        optimization_level=args.level
    )
