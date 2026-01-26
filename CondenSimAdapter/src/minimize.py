#!/usr/bin/env python3
"""
Minimize Module

Handles energy minimization of all-atom structures using AMBER/CHARMM force fields
with implicit solvent (GBn2). Uses gromacs pdb2gmx for topology generation and
top_to_softcore_system.py for OpenMM system creation and minimization.

Workflow:
1. Prepare input (backmap output or user PDB)
2. Generate topology using pdb2gmx (multi-component support via merge_topologies)
3. Generate structure using pdb2gmx
4. Run multi-step OpenMM minimization (Gaussian -> Softcore -> Standard)
5. Output minimized PDB

Force fields: AMBER99SB-ILDN, CHARMM36-jul2021
Solvent models: GBn2 (implicit), OBC2 (implicit)
"""

import os
import sys
import shutil
import subprocess
import glob
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from dataclasses import dataclass, field

import numpy as np
import MDAnalysis as mda

try:
    import openmm.unit as unit
    import openmm as mm
    from openmm.app import GromacsGroFile, GromacsTopFile, Simulation, PDBFile
    from openmm import Platform
except ImportError:
    import simtk.unit as unit
    import simtk.openmm as mm
    from simtk.openmm.app import GromacsGroFile, GromacsTopFile, Simulation, PDBFile
    from simtk.openmm import Platform

import click

from .cg import CGSimulationConfig, CGComponent, ComponentType
from .backmap import BackmapSimulator, SourceType
from .pdb2gmx_utils import (
    generate_all_atom_topology,
    merge_topologies,
    run_pdb2gmx_for_structure,
    load_components_from_yaml,
    load_config_from_yaml,
)
from ..forcefield.registry import (
    REGISTRY,
    list_force_fields,
    get_water_model,
    get_force_field_path,
)


# =============================================================================
# Force Field Constants
# =============================================================================

# Get available force fields from registry
AVAILABLE_FORCE_FIELDS = list_force_fields()

# Force field types (for validation)
AMBER_FORCEFIELD = "amber"
CHARMM_FORCEFIELD = "charmm"

FORCE_FIELD_FAMILIES = [AMBER_FORCEFIELD, CHARMM_FORCEFIELD]

# Registry for force field lookups
FORCE_FIELD_REGISTRY = REGISTRY

# Backward compatibility: FORCE_FIELD_TYPES is now the list of all available force fields
FORCE_FIELD_TYPES = AVAILABLE_FORCE_FIELDS

# Legacy mapping for pdb2gmx (kept for backward compatibility, prefer using registry)
# Supports both numbered names (1-a99SBdisp) and pdb2gmx names (a99SBdisp, charmm36m)
FORCEFIELD_TO_PDB2GMX = {
    # Numbered names
    "1-a99sbdisp": "a99SBdisp",
    "2-amber03wsc": "amber03wsc",
    "3-amber99sbws-stqp": "amber99sbws-STQp",
    "4-amber99sbws-stq": "amber99sbws-stq",
    "5-des-amber": "des-amber",
    "6-des-amber-sf1.0": "des-amber-SF1.0",
    "7-amber99sb-ildn": "amber99sb-ildn",
    "8-amber14sb": "amber14sb_parmbsc1",
    "9-charmm36m": "charmm36-jul2021",
    # pdb2gmx names (for backward compatibility)
    "amber99sb-ildn": "amber99sb-ildn",
    "amber99sbws-STQp": "amber99sbws-STQp",
    "amber99sbws-stq": "amber99sbws-stq",
    "amber03wsc": "amber03wsc",
    "amber14sb": "amber14sb_parmbsc1",
    "a99SBdisp": "a99SBdisp",
    "des-amber": "des-amber",
    "des-amber-SF1.0": "des-amber-SF1.0",
    "charmm36": "charmm36-jul2021",
    "charmm36m": "charmm36-jul2021",
}


# =============================================================================
# GB Model Constants
# =============================================================================

GB_GBN2 = "GBn2"
GB_OBC2 = "OBC2"

GB_MODELS = [GB_GBN2, GB_OBC2]


# =============================================================================
# Optimization Mode Constants
# =============================================================================

OPTIMIZATION_MODE_HIGH = "high"
OPTIMIZATION_MODE_MEDIUM = "medium"
OPTIMIZATION_MODE_LOW = "low"

OPTIMIZATION_MODE_PARAMS = {
    OPTIMIZATION_MODE_HIGH: {
        "steps": 4,
        "lambda_values": [0.65, 0.75, 0.85, 0.95],
    },
    OPTIMIZATION_MODE_MEDIUM: {
        "steps": 3,
        "lambda_values": [0.75, 0.85, 0.95],
    },
    OPTIMIZATION_MODE_LOW: {
        "steps": 2,
        "lambda_values": [0.85, 0.95],
    },
}


# =============================================================================
# Configuration and Result Classes
# =============================================================================


@dataclass
class MinimizeConfig:
    """Minimize simulation configuration.
    
    Attributes:
        forcefield_type: Force field name (e.g., 'amber99sb-ildn', 'charmm36')
        gb_model: GB model for implicit solvent (GBn2 or OBC2)
        box_resize_enabled: Whether to resize the box
        box_resize_dimensions: New box dimensions in nm [x, y, z]
        temperature: Temperature in K
        time_step: Time step in ps
        platform: OpenMM platform (CUDA, OpenCL, CPU)
        precision: Precision mode (single, mixed, double)
        gpu_id: GPU device index
        optimization_mode: Optimization level (high, medium, low)
        softcore_steps: Number of softcore steps
        softcore_lambda_values: Lambda values for softcore optimization
        tolerance: Minimization tolerance in kJ/(mol·nm)
        max_iterations: Maximum minimization iterations
        salt_conc: Salt concentration in M for implicit solvent
        nonbonded_cutoff: Non-bonded cutoff in nm
        disable_disulfide: Disable disulfide detection in pdb2gmx (-ss)
        his_type: Histidine type selection for pdb2gmx (-his), 0=HID, 1=HIE
    """
    
    # Force field selection (use registry for available options, numbered 1-9)
    forcefield_type: str = "1-a99SBdisp"  # Default: a99SBdisp (recommended)
    
    # GB model for implicit solvent
    gb_model: str = GB_GBN2  # GBn2 or OBC2
    
    # Box resize
    box_resize_enabled: bool = False
    box_resize_dimensions: Optional[List[float]] = None
    
    # OpenMM settings
    temperature: float = 300.0  # K
    time_step: float = 0.002  # ps
    platform: str = "CUDA"  # CUDA, OpenCL, CPU
    precision: str = "mixed"  # single, mixed, double
    gpu_id: int = 0  # GPU device index (only used when platform is CUDA or OpenCL)
    
    # Optimization settings
    optimization_mode: str = OPTIMIZATION_MODE_MEDIUM  # high, medium, low
    softcore_steps: int = 5
    softcore_lambda_values: List[float] = field(default_factory=lambda: [0.75, 0.85, 0.95])
    
    # Minimization parameters
    tolerance: float = 100.0  # kJ/(mol·nm)
    max_iterations: int = 5000
    
    # Implicit solvent parameters
    salt_conc: float = 0.15  # M
    
    # Non-bonded cutoff
    nonbonded_cutoff: float = 2.0  # nm
    
    # Explicit solvation parameters
    solvate_enabled: bool = False
    ion_concentration: float = 0.15  # M

    # pdb2gmx behavior
    disable_disulfide: bool = False
    his_type: Optional[int] = None
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        # Validate force field using registry
        is_valid, message = FORCE_FIELD_REGISTRY.validate(self.forcefield_type)
        if not is_valid:
            raise ValueError(message)
        
        # Validate GB model
        if self.gb_model not in GB_MODELS:
            raise ValueError(
                f"Unknown gb_model: {self.gb_model}. "
                f"Valid models: {GB_MODELS}"
            )
        
        # Validate optimization mode
        if self.optimization_mode not in OPTIMIZATION_MODE_PARAMS:
            raise ValueError(
                f"Unknown optimization_mode: {self.optimization_mode}. "
                f"Valid modes: {list(OPTIMIZATION_MODE_PARAMS.keys())}"
            )

        if self.his_type is not None and self.his_type not in (0, 1):
            raise ValueError("his_type must be 0 or 1 if provided")
    
    def set_optimization_mode(self, mode: str):
        """Set optimization mode and update lambda values accordingly.
        
        Args:
            mode: One of 'high', 'medium', 'low'
        """
        if mode not in OPTIMIZATION_MODE_PARAMS:
            raise ValueError(f"Unknown mode: {mode}. Valid: {list(OPTIMIZATION_MODE_PARAMS.keys())}")
        
        self.optimization_mode = mode
        params = OPTIMIZATION_MODE_PARAMS[mode]
        self.softcore_steps = params["steps"]
        self.softcore_lambda_values = params["lambda_values"]


@dataclass
class MinimizeResult:
    """Minimize simulation result."""
    
    success: bool
    output_pdb: str
    input_pdb: str
    errors: List[str] = field(default_factory=list)
    
    # Optimization step information
    step_info: str = ""  # e.g., "Medium (5 steps)"
    total_steps: int = 0
    
    # Optional: intermediate files
    intermediate_files: List[str] = field(default_factory=list)


# =============================================================================
# Minimize Simulator
# =============================================================================


class MinimizeSimulator:
    """Minimize simulator for AMBER/CHARMM all-atom force fields.
    
    This class orchestrates the minimize workflow:
    1. Load components from YAML config
    2. Generate multi-component topology using pdb2gmx + merge_topologies
    3. Generate structure using pdb2gmx from input PDB
    4. Run multi-step OpenMM minimization via worker process
    5. Return minimized structure
    
    Note: This class is an orchestrator only. All OpenMM calls are made
    in a separate worker process to avoid memory corruption issues.
    
    Force fields are managed through the ForceFieldRegistry, which provides
    access to pdb2gmx names and water models for each force field.
    """
    
    def __init__(self, 
                 minimize_config: Optional[MinimizeConfig] = None,
                 components: Optional[List[CGComponent]] = None,
                 system_name: str = "MinimizedSystem"):
        """
        Initialize minimize simulator.
        
        Args:
            minimize_config: MinimizeConfig (optional, uses defaults if not provided)
            components: List of CGComponent for multi-component support
            system_name: System name for [ system ] section in topology
        """
        self.minimize_config = minimize_config or MinimizeConfig()
        self.components = components or []
        self.system_name = system_name
        
        # Use registry to get force field info
        ff_info = FORCE_FIELD_REGISTRY.get_force_field(self.minimize_config.forcefield_type)
        if ff_info:
            self.forcefield_type = ff_info.family.upper()  # 'AMBER' or 'CHARMM' for OpenMM
            self.pdb2gmx_name = ff_info.pdb2gmx_name
            self.water_model = ff_info.water_model
        else:
            # Fallback for backward compatibility
            self.forcefield_type = "AMBER"
            self.pdb2gmx_name = "amber99sb-ildn"
            self.water_model = "tip3p"
    
    @classmethod
    def from_yaml(cls, yaml_path: str, 
                  minimize_config: Optional[MinimizeConfig] = None) -> 'MinimizeSimulator':
        """
        Create MinimizeSimulator from YAML configuration file.
        
        Args:
            yaml_path: Path to YAML configuration file
            minimize_config: Optional MinimizeConfig overrides
        
        Returns:
            MinimizeSimulator instance
        """
        # Load system_name and components from YAML
        system_name, components = load_config_from_yaml(yaml_path)
        
        return cls(minimize_config=minimize_config, components=components, system_name=system_name)
    
    def run(self, input_pdb: str, output_dir: Optional[str] = None) -> MinimizeResult:
        """
        Execute minimize workflow.
        
        Output structure:
        {system_name}_minimize/
        ├── minimize_final.pdb          # Final minimized structure
        ├── topol.top                   # Topology for final structure
        ├── topology/                   # Topology generation intermediate files
        ├── structure/                  # Structure processing intermediate files
        └── minimize/                   # Minimization intermediate files
        
        Args:
            input_pdb: Input PDB file path (backmap output or user PDB)
            output_dir: Output directory (optional, defaults to {system_name}_minimize)
            
        Returns:
            MinimizeResult object
        """
        result = MinimizeResult(
            success=False,
            output_pdb="",
            input_pdb=input_pdb,
            errors=[]
        )
        
        try:
            # Determine output directory
            if output_dir is None:
                # Use self.system_name for output directory
                output_dir = f"{self.system_name}_minimize"
            
            output_path = Path(output_dir).resolve()
            output_path.mkdir(parents=True, exist_ok=True)
            
            # Create subdirectories for intermediate files
            topology_dir = output_path / "topology"
            structure_dir = output_path / "structure"
            minimize_dir = output_path / "minimize"
            
            topology_dir.mkdir(exist_ok=True)
            structure_dir.mkdir(exist_ok=True)
            minimize_dir.mkdir(exist_ok=True)

            # Step counter for progress display
            step_num = 0

            # Copy force field folder to output directory for pdb2gmx
            ff_path = get_force_field_path(self.minimize_config.forcefield_type)
            if ff_path and Path(ff_path).exists():
                ff_folder_name = Path(ff_path).name
                target_ff_path = output_path / ff_folder_name
                if not target_ff_path.exists():
                    step_num += 1
                    click.echo(f"\n[Step {step_num}] Copying force field to output directory...")
                    shutil.copytree(ff_path, target_ff_path)
                    click.echo(f"  Force field: {ff_folder_name}")
                else:
                    step_num += 1
                    click.echo(f"\n[Step {step_num}] Force field already present: {ff_folder_name}")
            else:
                step_num += 1
                click.echo(f"\n[Step {step_num}] Warning: Force field folder not found for '{self.minimize_config.forcefield_type}'")
                click.echo(f"  pdb2gmx may fail if it cannot find the force field")

            # Step: Generate multi-component topology
            step_num += 1
            click.echo(f"\n[Step {step_num}] Generating topology from components...")
            if self.components:
                total_nmol = sum(comp.nmol for comp in self.components)
                his_repeat_count = max(total_nmol * 30, 30)
                component_topologies, water_model = generate_all_atom_topology(
                    self.components,
                    self.pdb2gmx_name,
                    topology_dir,
                    disable_disulfide=self.minimize_config.disable_disulfide,
                    his_type=self.minimize_config.his_type,
                    his_repeat_count=his_repeat_count
                )
                merged_top = merge_topologies(
                    component_topologies, 
                    topology_dir,
                    self.pdb2gmx_name,
                    water_model,
                    self.system_name
                )
                click.echo(f"  Merged topology: {merged_top.name}")
            else:
                raise ValueError("No components configured for topology generation")
            
            # 2. Generate structure from input PDB
            step_num += 1
            click.echo(f"\n[Step {step_num}] Generating structure from input PDB...")
            structure_gro = run_pdb2gmx_for_structure(
                Path(input_pdb),
                structure_dir,
                self.pdb2gmx_name,
                water_model="none",
                disable_disulfide=self.minimize_config.disable_disulfide,
                his_type=self.minimize_config.his_type,
                his_repeat_count=his_repeat_count
            )
            
            # 3. Box resize if enabled
            if self.minimize_config.box_resize_enabled:
                step_num += 1
                click.echo(f"\n[Step {step_num}] Resizing box...")
                resized_pdb = self.resize_box(
                    str(structure_gro), 
                    structure_dir,
                    self.minimize_config.box_resize_dimensions
                )
                structure_gro = resized_pdb
            
            # 4. Run OpenMM minimization in worker process (implicit solvent)
            step_num += 1
            click.echo(f"\n[Step {step_num}] Running implicit solvent minimization...")
            try:
                minimize_result = self.run_openmm_minimization(
                    str(structure_gro),
                    str(merged_top),
                    minimize_dir  # Intermediate files go here
                )
                
                # Copy final results to main directory
                final_pdb = output_path / "minimize_final.pdb"
                final_top = output_path / "topol.top"
                
                shutil.copy2(minimize_result['output_pdb'], final_pdb)
                shutil.copy2(merged_top, final_top)
                
                click.echo(f"  Final structure (implicit): {final_pdb.name}")
                click.echo(f"  Final topology: {final_top.name}")
                
                # Store optimization step information
                if 'step_info' in minimize_result:
                    result.step_info = minimize_result['step_info']
                if 'total_steps' in minimize_result:
                    result.total_steps = minimize_result['total_steps']
                
                # Generate plumed.dat for MDP components
                step_num += 1
                click.echo(f"\n[Step {step_num}] Checking for MDP components and generating plumed.dat...")
                from .plumed_generator import generate_plumed_for_minimize
                
                try:
                    plumed_generated = generate_plumed_for_minimize(
                        components=self.components,
                        topology_dir=topology_dir,
                        output_file=str(output_path / "plumed.dat"),
                        verbose=True
                    )
                    
                    if plumed_generated:
                        click.echo(f"  ✓ Generated plumed.dat for MDP restraints")
                    else:
                        click.echo(f"  ⊗ No MDP components found, skipping plumed.dat")
                except Exception as plumed_error:
                    click.echo(f"  Warning: Failed to generate plumed.dat: {plumed_error}")
                    # Don't fail the entire minimization if plumed generation fails
                    
            except Exception as e:
                click.echo(f"Error in OpenMM minimization: {e}", err=True)
                import traceback
                traceback.print_exc()
                raise
            
            # 5. Solvate if enabled
            if self.minimize_config.solvate_enabled:
                step_num += 1
                click.echo(f"\n[Step {step_num}] Explicit solvation...")
                solvate_dir = output_path / "solvate"
                solvate_dir.mkdir(exist_ok=True)
                
                # Solvate the optimized structure
                solvated_gro, solvated_top = self.solvate_system(
                    str(final_pdb),
                    str(final_top),
                    solvate_dir,
                    self.minimize_config.ion_concentration
                )
                
                # Update final outputs to solvated version
                final_gro_solvated = output_path / "minimize_final_solvated.gro"
                final_top_output = output_path / "topol.top"
                
                shutil.copy2(solvated_gro, final_gro_solvated)
                shutil.copy2(solvated_top, final_top_output)
                
                click.echo(f"  Final structure (solvated): {final_gro_solvated.name}")
                click.echo(f"  Final topology: {final_top_output.name}")
                
                # Build OpenMM system from solvated structure and save PDB
                step_num += 1
                click.echo(f"\n[Step {step_num}] Add chain label for structure for better visualization...")
                final_pdb_solvated = output_path / "minimize_final_solvated.pdb"
                
                # Use OpenMM to read GRO and TOP, build system, and save PDB
                # Force field files are already copied to output directory
                gro_file = GromacsGroFile(solvated_gro)
                top_file = GromacsTopFile(
                    solvated_top,
                    periodicBoxVectors=gro_file.getPeriodicBoxVectors(),
                    includeDir=str(output_path)  # Force field files are in output directory
                )
                # Build system (we don't need the system object for anything,
                # just need to validate topology)
                _ = top_file.createSystem(
                    nonbondedMethod=mm.app.PME,
                    nonbondedCutoff=1*unit.nanometer,
                    constraints=mm.app.HBonds
                )
                
                # Save coordinates to PDB
                with open(final_pdb_solvated, 'w') as f:
                    PDBFile.writeFile(
                        top_file.topology,
                        gro_file.getPositions(asNumpy=True),
                        f
                    )
                
                click.echo(f"  Final structure (PDB): {final_pdb_solvated.name}")
                
                result.output_pdb = str(final_pdb_solvated)
                result.intermediate_files = minimize_result.get('intermediate_files', [])
            else:
                result.output_pdb = str(final_pdb)
                result.intermediate_files = minimize_result.get('intermediate_files', [])
            
            result.success = True
            
        except Exception as e:
            result.errors.append(str(e))
            import traceback
            result.errors.append(traceback.format_exc())
        
        return result
    
    def resize_box(self, pdb_path: str, output_dir: Path,
                   new_dimensions: Optional[List[float]]) -> Path:
        """
        Resize box dimensions.
        
        Args:
            pdb_path: Input PDB path
            output_dir: Output directory
            new_dimensions: New box dimensions in nm [x, y, z]
            
        Returns:
            Path to resized PDB
        """
        if new_dimensions is None:
            raise ValueError("box_resize_dimensions must be provided when box_resize_enabled=True")
        
        # Import helper functions from backmap module
        from .backmap import BackmapSimulator
        
        backmap_sim = BackmapSimulator()
        return backmap_sim.resize_box(pdb_path, output_dir, new_dimensions)
    
    def solvate_system(self, structure_pdb: str, topology_top: str, 
                       solvate_dir: Path, ion_concentration: float) -> Tuple[str, str]:
        """
        Solvate system with explicit water and add ions.
        
        Args:
            structure_pdb: Input structure file (PDB)
            topology_top: Input topology file
            solvate_dir: Directory for solvation
            ion_concentration: Ion concentration in M
            
        Returns:
            Tuple of (solvated_gro, solvated_top)
        """
        # Determine water model based on force field
        ff_info = FORCE_FIELD_REGISTRY.get_force_field(self.minimize_config.forcefield_type)
        
        # Use registry-provided water model and gmx solvate -cs string
        water_cs = ff_info.water_model if ff_info else "tip3p"
        solvate_cs = ff_info.solvate_cs if ff_info else "spc216"
        click.echo(f"  Water model: {water_cs} (gmx -cs: {solvate_cs})")
        click.echo(f"  Ion concentration: {ion_concentration} M")
        
        # Convert PDB to GRO using gmx editconf
        system_gro = solvate_dir / "system.gro"
        system_top = solvate_dir / "topol.top"
        
        # Convert PDB to GRO
        editconf_cmd = [
            'gmx', 'editconf',
            '-f', str(structure_pdb),
            '-o', str(system_gro)
        ]
        
        result = subprocess.run(editconf_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"gmx editconf failed: {result.stderr}")
        
        shutil.copy2(topology_top, system_top)

        # Copy force field folder to solvate directory for #include resolution
        ff_path = get_force_field_path(self.minimize_config.forcefield_type)
        if ff_path:
            ff_folder_name = Path(ff_path).name
            target_ff_path = solvate_dir / ff_folder_name
            if not target_ff_path.exists():
                shutil.copytree(ff_path, target_ff_path)
        
        # Create em.mdp for grompp
        em_mdp_content = self._create_em_mdp()
        em_mdp = solvate_dir / "em_steep.mdp"
        with open(em_mdp, 'w') as f:
            f.write(em_mdp_content)
        
        # Step 1: Solvate with gmx solvate
        click.echo(f"  Adding water box...")
        solvate_cmd = [
            'gmx', 'solvate',
            '-cp', str(system_gro),
            '-cs', solvate_cs,
            '-o', str(solvate_dir / 'solvated.gro'),
            '-p', str(system_top)
        ]
        
        result = subprocess.run(solvate_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"gmx solvate failed: {result.stderr}")
        
        click.echo(f"  ✓ Water box added")
        
        # Step 2: Generate ions.tpr for genion
        click.echo(f"  Preparing for ion addition...")
        grompp_cmd = [
            'gmx', 'grompp',
            '-f', str(em_mdp),
            '-c', str(solvate_dir / 'solvated.gro'),
            '-p', str(system_top),
            '-o', str(solvate_dir / 'ions.tpr'),
            '-maxwarn', '2'
        ]
        
        result = subprocess.run(grompp_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"gmx grompp failed: {result.stderr}")
        
        # Step 3: Add ions with gmx genion
        click.echo(f"  Adding ions...")
        genion_cmd = [
            'gmx', 'genion',
            '-s', str(solvate_dir / 'ions.tpr'),
            '-o', str(solvate_dir / 'system_ions.gro'),
            '-p', str(system_top),
            '-conc', str(ion_concentration),
            '-neutral'
        ]
        
        result = subprocess.run(
            genion_cmd, 
            input='SOL\n',  # Select SOL group for ion replacement
            capture_output=True, 
            text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"gmx genion failed: {result.stderr}")
        
        click.echo(f"  ✓ Ions added")
        
        solvated_gro = str(solvate_dir / 'system_ions.gro')
        solvated_top = str(system_top)
        
        return solvated_gro, solvated_top
    
    def _create_em_mdp(self) -> str:
        """
        Create energy minimization MDP file content for solvation.
        
        Returns:
            MDP file content as string
        """
        em_mdp = """
; Energy minimization for solvated system
integrator          = steep
nsteps              = 50000
emtol               = 2000
emstep              = 0.01

; Output control
nstlog              = 100
nstenergy           = 100

; Neighbor searching
cutoff-scheme       = Verlet
nstlist             = 10
ns_type             = grid
rlist               = 1.2

; Electrostatics
coulombtype         = reaction-field
rcoulomb            = 1.1
epsilon_r           = 15

; van der Waals
vdwtype             = Cut-off
vdw_modifier        = Force-switch
rvdw                = 1.1
rvdw_switch         = 0.9

; Temperature and pressure coupling
Tcoupl              = no
Pcoupl              = no

; Velocity generation
gen_vel             = no

; Constraints
constraints         = none

; Periodic boundary conditions
pbc                 = xyz
"""
        return em_mdp
    
    def run_openmm_minimization(self, structure_gro: str, topology_top: str,
                                 minimize_dir: Path) -> Dict:
        """
        Run multi-step OpenMM minimization using subprocess.
        
        This method runs the minimization in a separate Python process
        to avoid memory corruption issues that can occur when repeatedly
        creating OpenMM systems in the same process.
        
        All intermediate files are stored in minimize_dir.
        
        Args:
            structure_gro: Input structure file (GRO)
            topology_top: Input topology file (TOP)
            minimize_dir: Directory for minimization files (all intermediate files go here)
        
        Returns:
            Dict with 'output_pdb' and 'intermediate_files'
        
        Steps:
        1. Gaussian repulsion
        2. Softcore (multiple steps)
        3. Standard force field (NONBONDED_STANDARD)
        """
        minimize_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy files to minimize directory
        gro_path = Path(structure_gro)
        top_path = Path(topology_top)
        
        conf_gro = minimize_dir / "conf.gro"
        minimize_top = minimize_dir / "topol.top"
        
        shutil.copy2(str(gro_path), str(conf_gro))
        shutil.copy2(str(top_path), str(minimize_top))
        
        # Copy force field folder to minimize directory for #include resolution
        ff_path = get_force_field_path(self.minimize_config.forcefield_type)
        if ff_path:
            ff_folder_name = Path(ff_path).name
            target_ff_path = minimize_dir / ff_folder_name
            if not target_ff_path.exists():
                shutil.copytree(ff_path, target_ff_path)
        
        # Find the minimize_worker.py script
        script_dir = Path(__file__).parent
        minimize_script = script_dir / "minimize_worker.py"
        
        if not minimize_script.exists():
            raise FileNotFoundError(f"Cannot find minimize_worker.py at {minimize_script}")
        
        # Run minimization in subprocess
        click.echo(f"  Running OpenMM minimization in isolated process...")
        
        # Build command arguments
        cmd_args = [
            sys.executable, str(minimize_script),
            '--input-gro', str(conf_gro),
            '--input-top', str(minimize_top),
            '--output', str(minimize_dir),
            '--device', self.minimize_config.platform.lower(),
            '--gpu-id', str(self.minimize_config.gpu_id),
            '--iter', str(self.minimize_config.max_iterations),
            '--tolerance', str(self.minimize_config.tolerance),
            '--level', self.minimize_config.optimization_mode,
            '--ff-type', self.forcefield_type.lower(),  # 'amber' or 'charmm' for OpenMM
            '--ff-name', self.pdb2gmx_name,  # Full pdb2gmx name for reference
            '--gb-model', self.minimize_config.gb_model,
            '--salt-conc', str(self.minimize_config.salt_conc),
            '--cutoff', str(self.minimize_config.nonbonded_cutoff),
        ]
        
        result = subprocess.run(
            cmd_args,
            stdout=None,  # Show output in real-time
            stderr=subprocess.PIPE,
            text=True
        )
        
        if result.returncode != 0:
            error_msg = f"Minimization failed with return code {result.returncode}"
            if result.stderr:
                error_msg += f"\n{result.stderr}"
            raise RuntimeError(error_msg)
        
        # Parse output to find generated files
        output_pdb = minimize_dir / "minimize_final.pdb"
        if not output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {output_pdb}")
        
        # Collect intermediate files
        intermediate_files = []
        pdb_files = sorted(minimize_dir.glob("*.pdb"))
        for pdb in pdb_files:
            if pdb.name != "minimize_final.pdb":
                intermediate_files.append(str(pdb))
        
        return {
            'output_pdb': str(output_pdb),
            'intermediate_files': intermediate_files
        }
