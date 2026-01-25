#!/usr/bin/env python3
"""
PACE Optimization Module

Handles PACE force field optimization of backmapped structures using OpenMM.
Supports multi-component systems and multi-step optimization (Gaussian → Softcore → Standard).
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
# Use openmm (compatible with simtk.openmm used by PACE2openmm)
try:
    import openmm.unit as unit
    import openmm as mm
    from openmm.app import GromacsGroFile, Simulation, PDBFile
    from openmm import Platform
except ImportError:
    # Fallback to simtk.openmm for older versions
    import simtk.unit as unit
    import simtk.openmm as mm
    from simtk.openmm.app import GromacsGroFile, Simulation, PDBFile
    from simtk.openmm import Platform
import click

from .cg import CGSimulationConfig, CGComponent, ComponentType
from .backmap import BackmapSimulator, SourceType

# Import PACE2openmm
# Add src directory (current directory) to path
src_path = Path(__file__).parent
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Import PACE2openmm module
import importlib.util
spec = importlib.util.spec_from_file_location("PACE2openmm", src_path / "PACE2openmm.py")
PACE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PACE)


# =============================================================================
# PACE Force Field Type Constants
# =============================================================================

# PACE force field variants
PACE_FORCEFIELD_NEW = "pace-new"  # Standard PACE force field
PACE_FORCEFIELD_ASM = "pace-asm"  # Alternative PACE variant for condensate structures

# List of valid PACE force field types
PACE_FORCEFIELD_TYPES = [PACE_FORCEFIELD_NEW, PACE_FORCEFIELD_ASM]


# Optimization mode constants
OPTIMIZATION_MODE_HIGH = "high"
OPTIMIZATION_MODE_MEDIUM = "medium"
OPTIMIZATION_MODE_LOW = "low"

# Mapping from optimization mode to softcore parameters
OPTIMIZATION_MODE_PARAMS = {
    OPTIMIZATION_MODE_HIGH: {
        "steps": 7,
        "alpha_values": [0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.01],
    },
    OPTIMIZATION_MODE_MEDIUM: {
        "steps": 5,
        "alpha_values": [0.2, 0.1, 0.05, 0.025, 0.01],
    },
    OPTIMIZATION_MODE_LOW: {
        "steps": 3,
        "alpha_values": [0.15, 0.05, 0.025],
    },
}


@dataclass
class PaceOptConfig:
    """PACE optimization configuration"""
    
    # Force field type selection
    forcefield_type: str = PACE_FORCEFIELD_NEW  # pace-new or pace-asm
    
    # Box resize
    box_resize_enabled: bool = False
    box_resize_dimensions: Optional[List[float]] = None

    # OpenMM settings
    temperature: float = 325.0  # K
    pressure: float = 1.0  # bar
    time_step: float = 0.004  # ps
    platform: str = "CUDA"  # CUDA, OpenCL, CPU
    precision: str = "double"  # single, mixed, double

    # Softcore optimization settings
    optimization_mode: str = OPTIMIZATION_MODE_HIGH  # high, medium, low
    softcore_steps: int = 7
    softcore_alpha_values: List[float] = field(default_factory=lambda: [0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.01])  # VDW softness (LJ)
    softcore_alpha_coul_values: List[float] = field(default_factory=lambda: [0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.01])  # Electrostatic softness
    softcore_alpha_exc_values: List[float] = field(default_factory=lambda: [0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.01])  # Exception force softness (now same as VDW alpha)

    # GROMACS optimization
    gromacs_enabled: bool = True
    gromacs_num_cpus: int = 1
    gromacs_constraints: str = "hbonds"  # hbonds, none

    # Epsilon_r for PACE
    epsilon_r: float = 15.0

    def __post_init__(self):
        """Validate forcefield_type after initialization."""
        if self.forcefield_type not in PACE_FORCEFIELD_TYPES:
            raise ValueError(
                f"Unknown forcefield_type: {self.forcefield_type}. "
                f"Valid types: {PACE_FORCEFIELD_TYPES}"
            )

    def set_optimization_mode(self, mode: str):
        """Set softcore parameters based on optimization mode.

        Args:
            mode: One of 'high', 'medium', 'low'
        """
        mode = mode.lower()
        if mode not in OPTIMIZATION_MODE_PARAMS:
            raise ValueError(f"Unknown optimization mode: {mode}. Valid modes: {list(OPTIMIZATION_MODE_PARAMS.keys())}")

        params = OPTIMIZATION_MODE_PARAMS[mode]
        self.optimization_mode = mode
        self.softcore_steps = params["steps"]
        self.softcore_alpha_values = params["alpha_values"]
        self.softcore_alpha_coul_values = params["alpha_values"]
        self.softcore_alpha_exc_values = params["alpha_values"]


@dataclass
class PaceOptResult:
    """PACE optimization result"""
    success: bool
    output_pdb: str
    input_pdb: str
    errors: List[str] = field(default_factory=list)
    intermediate_files: List[str] = field(default_factory=list)


class PaceOptSimulator:
    """PACE optimization simulator supporting multiple force field variants."""
    
    def __init__(self, config: Optional[CGSimulationConfig] = None, 
                 pace_opt_config: Optional[PaceOptConfig] = None):
        """
        Initialize PACE optimization simulator
        
        Args:
            config: CGSimulationConfig (required for multi-component support)
            pace_opt_config: PaceOptConfig (optional, uses defaults if not provided)
        """
        self.config = config
        self.pace_opt_config = pace_opt_config or PaceOptConfig()
        
        # Backmap simulator for detecting input type
        self.backmap_simulator = BackmapSimulator(config=config)
        
        # Store forcefield type for easy access
        self.forcefield_type = self.pace_opt_config.forcefield_type
        
        # Validate forcefield type
        if self.forcefield_type not in PACE_FORCEFIELD_TYPES:
            raise ValueError(
                f"Unknown forcefield_type: {self.forcefield_type}. "
                f"Valid types: {PACE_FORCEFIELD_TYPES}"
            )
        
        # Dispatch table for topology generation based on forcefield type
        self._topology_generators = {
            PACE_FORCEFIELD_NEW: self._generate_topology_pace_new,
            PACE_FORCEFIELD_ASM: self._generate_topology_pace_asm,
        }
        
        # Dispatch table for structure conversion based on forcefield type
        self._structure_converters = {
            PACE_FORCEFIELD_NEW: self._convert_structure_pace_new,
            PACE_FORCEFIELD_ASM: self._convert_structure_pace_asm,
        }
    
    def _get_topology_generator(self):
        """Get the topology generator function for the current forcefield type."""
        generator = self._topology_generators.get(self.forcefield_type)
        if generator is None:
            raise NotImplementedError(
                f"Topology generator not implemented for forcefield type: {self.forcefield_type}"
            )
        return generator
    
    def _get_structure_converter(self):
        """Get the structure converter function for the current forcefield type."""
        converter = self._structure_converters.get(self.forcefield_type)
        if converter is None:
            raise NotImplementedError(
                f"Structure converter not implemented for forcefield type: {self.forcefield_type}"
            )
        return converter
    
    def run(self, input_path: str, config_path: Optional[str] = None, 
            output_dir: Optional[str] = None) -> PaceOptResult:
        """
        Execute PACE optimization workflow
        
        Args:
            input_path: Input path (backmap output directory or PDB file)
            config_path: Config YAML path (optional)
            output_dir: Output directory (optional)
            
        Returns:
            PaceOptResult object
        """
        result = PaceOptResult(
            success=False,
            output_pdb="",
            input_pdb=input_path,
            errors=[]
        )
        
        try:
            # 1. Prepare input
            click.echo(f"\n[1/6] Preparing input...")
            prepared_input = self.prepare_input(input_path, config_path)
            result.input_pdb = prepared_input['pdb_path']
            
            # Determine output directory
            if output_dir is None:
                if self.config:
                    output_dir = f"{self.config.system_name}_pace_opt"
                else:
                    output_dir = "pace_opt_output"
            
            output_path = Path(output_dir).resolve()
            output_path.mkdir(parents=True, exist_ok=True)
            
            # 2. Box resize if enabled
            if self.pace_opt_config.box_resize_enabled:
                click.echo(f"\n[2/6] Resizing box...")
                resized_pdb = self.resize_box(
                    prepared_input['pdb_path'], 
                    output_path,
                    self.pace_opt_config.box_resize_dimensions
                )
                prepared_input['pdb_path'] = resized_pdb
            
            # 3. Generate topologies for all components
            click.echo(f"\n[3/6] Generating PACE topologies...")
            topology_path = self.generate_topologies(output_path)
            
            # 4. Process structure with pdb2gmx
            click.echo(f"\n[4/6] Processing structure with pdb2gmx...")
            try:
                processed_gro = self.run_pdb2gmx(prepared_input['pdb_path'], output_path)
            except Exception as e:
                click.echo(f"Error in pdb2gmx: {e}", err=True)
                import traceback
                traceback.print_exc()
                raise
            
            # Force garbage collection after pdb2gmx to avoid memory issues
            import gc
            gc.collect()
            
            # 5. Run OpenMM optimization
            click.echo(f"\n[5/6] Running OpenMM optimization...")
            try:
                optimized_pdb = self.run_openmm_optimization(
                    processed_gro, 
                    topology_path, 
                    output_path
                )
            except Exception as e:
                click.echo(f"Error in OpenMM optimization: {e}", err=True)
                import traceback
                traceback.print_exc()
                raise
            result.output_pdb = optimized_pdb
            result.intermediate_files.append(optimized_pdb)
            
            # 6. Run GROMACS optimization (optional)
            if self.pace_opt_config.gromacs_enabled:
                click.echo(f"\n[6/6] Running GROMACS optimization...")
                final_pdb = self.run_gromacs_optimization(
                    optimized_pdb,
                    topology_path,
                    output_path
                )
                result.output_pdb = final_pdb
                result.intermediate_files.append(final_pdb)
            
            result.success = True
            
        except Exception as e:
            result.errors.append(str(e))
            import traceback
            result.errors.append(traceback.format_exc())
        
        return result
    
    def prepare_input(self, input_path: str, config_path: Optional[str] = None) -> Dict:
        """
        Prepare input structure (similar to backmap)
        
        Returns:
            Dict with 'pdb_path' and 'config' keys
        """
        input_path_obj = Path(input_path)
        
        # Load config if needed
        if not self.config and config_path:
            self.config = CGSimulationConfig.from_yaml(config_path)
        
        # Check if it's a backmap output directory (has .aa.pdb files)
        if input_path_obj.is_dir():
            # Look for backmapped structure files
            candidates = [
                "backmapped_aa_final_ter.pdb",
                "backmapped_aa_final.pdb",
                "backmapped_aa_box.pdb",
                "final.aa.pdb",  # Alternative naming
            ]
            
            # Try exact matches first
            for candidate in candidates:
                candidate_path = input_path_obj / candidate
                if candidate_path.exists():
                    return {
                        'pdb_path': str(candidate_path.resolve()),
                        'config': self.config
                    }
            
            # Try pattern match for .aa.pdb files
            aa_pdb_files = list(input_path_obj.glob("*.aa.pdb"))
            if aa_pdb_files:
                # Use the first one found, prefer final.aa.pdb if exists
                final_aa = input_path_obj / "final.aa.pdb"
                if final_aa.exists():
                    return {
                        'pdb_path': str(final_aa.resolve()),
                        'config': self.config
                    }
                return {
                    'pdb_path': str(aa_pdb_files[0].resolve()),
                    'config': self.config
                }
        
        # If it's a file, treat as user provided PDB
        elif input_path_obj.is_file():
            if not config_path:
                raise ValueError("User provided mode requires config.yaml via -f option")
            
            if not self.config:
                self.config = CGSimulationConfig.from_yaml(config_path)
            
            if not input_path_obj.exists():
                raise FileNotFoundError(f"Input PDB not found: {input_path}")
            
            return {
                'pdb_path': str(input_path_obj.resolve()),
                'config': self.config
            }
        
        else:
            raise FileNotFoundError(f"Input path does not exist: {input_path}")
    
    def resize_box(self, pdb_path: str, output_dir: Path, 
                   new_dimensions: Optional[List[float]]) -> str:
        """
        Resize box dimensions (reuse logic from backmap)
        
        Args:
            pdb_path: Input PDB path
            output_dir: Output directory
            new_dimensions: New box dimensions in nm [x, y, z]
            
        Returns:
            Path to resized PDB
        """
        if new_dimensions is None:
            raise ValueError("box_resize_dimensions must be provided when box_resize_enabled=True")
        
        # Import helper functions from old_code/backmap.py
        old_code_path = Path(__file__).parent.parent / "old_code"
        if str(old_code_path) not in sys.path:
            sys.path.insert(0, str(old_code_path))
        import backmap as old_backmap
        
        u = mda.Universe(pdb_path)
        protein_dims = old_backmap.get_protein_dimensions(u)  # Angstroms
        protein_dims_nm = protein_dims / 10.0  # Convert to nm
        
        new_dims = np.array(new_dimensions)
        
        if not old_backmap.check_protein_fits_in_box(protein_dims_nm, new_dims):
            # Search trajectory if available
            # This would require access to CG output directory
            raise ValueError(
                f"Protein dimensions ({protein_dims_nm}) exceed box size ({new_dims}). "
                "Please adjust box dimensions or provide a suitable frame."
            )
        
        # Resize box
        resized_pdb = output_dir / "resized_input.pdb"
        final_universe = mda.Universe(pdb_path)
        final_universe.dimensions = list(new_dims) + [90, 90, 90]
        old_backmap.write_pdb_with_bfactors(final_universe, str(resized_pdb))
        
        # Write CRYST1 with new dimensions
        final_dims_A = new_dims * 10.0  # Convert nm to Angstroms
        cryst1_fixed_pdb = output_dir / "resized_input_cryst1.pdb"
        old_backmap.write_cryst1_only(str(resized_pdb), str(cryst1_fixed_pdb), final_dims_A)
        
        return str(cryst1_fixed_pdb)
    
    def generate_topologies(self, output_dir: Path) -> str:
        """
        Generate PACE topologies for all components and merge them.
        Uses the appropriate topology generator based on forcefield_type.
        
        Returns:
            Path to merged topology file
        """
        if not self.config:
            raise ValueError("Config is required for topology generation")
        
        # Get the appropriate topology generator for this forcefield type
        topology_generator = self._get_topology_generator()
        
        # Generate topologies using the selected method
        component_topologies = topology_generator(output_dir)
        
        # Merge topologies
        topology_dir = output_dir / "topology"
        merged_topology = self._merge_topologies(component_topologies, topology_dir)
        
        return str(merged_topology)
    
    # =============================================================================
    # PACE-NEW Topology Generation (Current Implementation)
    # =============================================================================
    
    def _generate_topology_pace_new(self, output_dir: Path) -> List[Dict]:
        """
        Generate PACE-NEW topologies for all components.
        This is the original implementation using prepare_peptide.py.
        
        Args:
            output_dir: Output directory for topology files
            
        Returns:
            List of component topology dictionaries
        """
        click.echo(f"  Using PACE-NEW force field topology generator")
        
        topology_dir = output_dir / "topology"
        topology_dir.mkdir(parents=True, exist_ok=True)
        
        component_topologies = []
        
        # Generate topology for each component
        for component in self.config.components:
            click.echo(f"  Generating topology for component: {component.name}")
            
            # Get sequence
            sequence = self._get_component_sequence(component)
            
            # Generate topology
            comp_top_pdb, comp_top_top = self._generate_component_topology(
                sequence, 
                component.name,
                topology_dir
            )
            
            # Modify molecule name in topology
            self._modify_topology_molecule_name(comp_top_top, component.name)
            
            component_topologies.append({
                'name': component.name,
                'topology': comp_top_top,
                'nmol': component.nmol
            })
        
        return component_topologies
    
    def _convert_structure_pace_new(self, pdb_path: str, output_dir: Path) -> str:
        """
        Convert input structure for PACE-NEW.
        Currently this is a pass-through (no conversion needed).
        
        Args:
            pdb_path: Input PDB path
            output_dir: Output directory
            
        Returns:
            Path to converted PDB file
        """
        # PACE-NEW doesn't require structure conversion
        # Just copy the file if needed
        return pdb_path
    
    # =============================================================================
    # PACE-ASM Topology Generation (Now using prepare_pace_asm.py)
    # =============================================================================
    
    def _generate_topology_pace_asm(self, output_dir: Path) -> List[Dict]:
        """
        Generate PACE-ASM topologies for all components.
        
        PACE-ASM is an alternative PACE variant for condensate structures.
        This method uses prepare_pace_asm.py which has the same interface as
        prepare_peptide.py (PACE-NEW), so the logic is nearly identical.
        
        Args:
            output_dir: Output directory for topology files
            
        Returns:
            List of component topology dictionaries
        """
        click.echo(f"  Using PACE-ASM force field topology generator")
        
        topology_dir = output_dir / "topology"
        topology_dir.mkdir(parents=True, exist_ok=True)
        
        component_topologies = []
        
        # Generate topology for each component
        for component in self.config.components:
            click.echo(f"  Generating topology for component: {component.name}")
            
            # Get sequence
            sequence = self._get_component_sequence(component)
            
            # Generate topology using PACE-ASM script
            comp_top_pdb, comp_top_top = self._generate_asm_component_topology(
                sequence, 
                component.name,
                topology_dir
            )
            
            # Modify molecule name in topology
            self._modify_topology_molecule_name(comp_top_top, component.name)
            
            component_topologies.append({
                'name': component.name,
                'topology': comp_top_top,
                'nmol': component.nmol
            })
        
        return component_topologies
    
    def _generate_asm_component_topology(self, sequence: str, component_name: str, 
                                          output_dir: Path) -> Tuple[str, str]:
        """
        Generate PACE-ASM topology for a single component using prepare_pace_asm.py.
        
        This has the same interface as _generate_component_topology() for PACE-NEW.
        
        Args:
            sequence: Amino acid sequence
            component_name: Name of the component
            output_dir: Output directory
            
        Returns:
            Tuple of (PDB path, topology path)
        """
        import CondenSimAdapter
        pace_asm_builder_dir = Path(CondenSimAdapter.__path__[0]) / "pace_asm_top_builder"
        prepare_script = pace_asm_builder_dir / "prepare_pace_asm.py"
        
        if not prepare_script.exists():
            raise FileNotFoundError(f"prepare_pace_asm.py not found: {prepare_script}")
        
        original_cwd = os.getcwd()
        comp_topology_dir = output_dir / component_name
        comp_topology_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            os.chdir(pace_asm_builder_dir)
            
            cmd = [sys.executable, str(prepare_script), sequence, "--name", component_name]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(
                    f"Error running prepare_pace_asm.py:\n"
                    f"STDOUT: {result.stdout}\n"
                    f"STDERR: {result.stderr}"
                )
            
            # Move output files (same structure as PACE-NEW)
            out_dir = pace_asm_builder_dir / "out"
            if out_dir.exists():
                source_top = out_dir / "topol.top"
                dest_top = comp_topology_dir / f"{component_name}.top"
                source_pdb = out_dir / "system.pdb"
                dest_pdb = comp_topology_dir / f"{component_name}.pdb"
                
                if source_top.exists():
                    shutil.move(str(source_top), str(dest_top))
                if source_pdb.exists():
                    shutil.move(str(source_pdb), str(dest_pdb))
                
                shutil.rmtree(out_dir)
            
            return str(dest_pdb), str(dest_top)
            
        finally:
            os.chdir(original_cwd)
    
    def _convert_structure_pace_asm(self, pdb_path: str, output_dir: Path) -> str:
        """
        Convert input structure for PACE-ASM.
        
        PACE-ASM may require structure conversion from condensate representation
        to standard PACE representation.
        
        NOTE: This method follows the same pattern as _convert_structure_pace_new().
        The difference is in HOW the conversion is done.
        
        TODO: Implement structure conversion for PACE-ASM.
        
        Args:
            pdb_path: Input PDB path
            output_dir: Output directory
            
        Returns:
            Path to converted PDB file
            
        Raises:
            NotImplementedError: This method is not yet implemented
        """
        click.echo(f"  Converting structure for PACE-ASM")
        click.echo(f"  WARNING: PACE-ASM structure conversion is not yet implemented")
        
        # TODO: Implement PACE-ASM structure conversion
        #
        # This should handle condensate-specific conversions:
        #   - Map condensate atoms to PACE representation
        #   - Handle non-standard residues (RNA, DNA, small molecules)
        #   - Ensure proper protonation states for PACE compatibility
        #
        # return converted_pdb_path
        
        raise NotImplementedError(
            f"PACE-ASM structure conversion is not yet implemented. "
            f"Please implement _convert_structure_pace_asm() in pace_opt.py"
        )
    
    def _get_component_sequence(self, component: CGComponent) -> str:
        """Get sequence for a component"""
        if component.type == ComponentType.IDP:
            if not component.ffasta:
                raise ValueError(f"Component {component.name}: IDP requires ffasta")
            
            sequences = self._read_fasta_sequence(component.ffasta)
            if component.name in sequences:
                return sequences[component.name]
            elif len(sequences) == 1:
                return list(sequences.values())[0]
            else:
                raise ValueError(f"Could not determine sequence for {component.name}")
        
        elif component.type == ComponentType.MDP:
            if not component.fpdb:
                raise ValueError(f"Component {component.name}: MDP requires fpdb")
            
            return self._extract_sequence_from_pdb(component.fpdb)
        
        else:
            raise ValueError(f"Unknown component type: {component.type}")
    
    def _read_fasta_sequence(self, fasta_path: str) -> Dict[str, str]:
        """Read sequences from FASTA file"""
        sequences = {}
        if not os.path.exists(fasta_path):
            return sequences
        
        name = None
        seq_chunks = []
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if name is not None:
                        sequences[name] = ''.join(seq_chunks)
                    name = line[1:].strip()
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
        if name is not None:
            sequences[name] = ''.join(seq_chunks)
        return sequences
    
    def _extract_sequence_from_pdb(self, pdb_path: str) -> str:
        """Extract sequence from PDB file"""
        u = mda.Universe(pdb_path)
        protein = u.select_atoms("protein")
        residues = protein.residues
        
        three_to_one = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
            'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
            'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
            'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        
        sequence = ""
        for residue in residues:
            resname = residue.resname.strip()
            if resname in three_to_one:
                sequence += three_to_one[resname]
            else:
                sequence += 'X'
        
        return sequence
    
    def _generate_component_topology(self, sequence: str, component_name: str, 
                                     output_dir: Path) -> Tuple[str, str]:
        """Generate PACE topology for a single component"""
        import CondenSimAdapter
        pace_builder_dir = Path(CondenSimAdapter.__path__[0]) / "pace_top_builder"
        prepare_script = pace_builder_dir / "prepare_peptide.py"
        
        if not prepare_script.exists():
            raise FileNotFoundError(f"prepare_peptide.py not found: {prepare_script}")
        
        original_cwd = os.getcwd()
        comp_topology_dir = output_dir / component_name
        comp_topology_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            os.chdir(pace_builder_dir)
            
            cmd = [sys.executable, str(prepare_script), sequence, "--name", component_name]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(
                    f"Error running prepare_peptide.py:\n"
                    f"STDOUT: {result.stdout}\n"
                    f"STDERR: {result.stderr}"
                )
            
            # Move output files
            out_dir = pace_builder_dir / "out"
            if out_dir.exists():
                source_top = out_dir / "topol.top"
                dest_top = comp_topology_dir / f"{component_name}.top"
                source_pdb = out_dir / "system.pdb"
                dest_pdb = comp_topology_dir / f"{component_name}.pdb"
                
                if source_top.exists():
                    shutil.move(str(source_top), str(dest_top))
                if source_pdb.exists():
                    shutil.move(str(source_pdb), str(dest_pdb))
                
                shutil.rmtree(out_dir)
            
            return str(dest_pdb), str(dest_top)
            
        finally:
            os.chdir(original_cwd)
    
    def _modify_topology_molecule_name(self, topol_path: str, molecule_name: str):
        """Modify molecule name in topology file"""
        with open(topol_path, 'r') as f:
            lines = f.readlines()
        
        in_moleculetype = False
        for i, line in enumerate(lines):
            if line.strip().startswith('[ moleculetype ]'):
                in_moleculetype = True
                continue
            if in_moleculetype and line.strip() and not line.startswith(';'):
                # Replace molecule name
                parts = line.strip().split()
                if len(parts) >= 1:
                    parts[0] = molecule_name
                    lines[i] = ' '.join(parts) + '\n'
                break
        
        with open(topol_path, 'w') as f:
            f.writelines(lines)
    
    def _merge_topologies(self, component_topologies: List[Dict], output_dir: Path) -> Path:
        """Merge multiple component topologies into one.

        PACE topology structure:
        1. #include "pace-new.ff/forcefield.itp"
        2. [ moleculetype ] section (molecule definition)
        3. ... (bonds, angles, dihedrals, etc.)
        4. ; Include water topology  ← 遇到这行停止，上一行是结尾
        5. [ molecules ] section

        Args:
            component_topologies: List of component topology dicts with 'name', 'topology', 'nmol'
            output_dir: Output directory for merged topology

        Returns:
            Path to merged topology file
        """
        merged_top = output_dir / "PACE.top"

        # Find the water topology marker line index
        water_marker = "; Include water topology"

        # Extract forcefield include and molecule definitions from each component
        forcefield_include = None
        molecule_definitions = {}

        for comp in component_topologies:
            with open(comp['topology'], 'r') as f:
                lines = f.readlines()

            found_forcefield_include = False
            in_moleculetype = False
            molecule_lines = []
            found_water_marker = False

            for i, line in enumerate(lines):
                # Skip empty lines and comment lines for detection
                stripped = line.strip()

                # 1. Find the forcefield include line
                if not found_forcefield_include and stripped.startswith('#include'):
                    # Keep the forcefield include (should be #include "pace-new.ff/forcefield.itp")
                    forcefield_include = line
                    found_forcefield_include = True
                    continue

                # 2. Detect start: [ moleculetype ]
                if stripped.startswith('[ moleculetype ]'):
                    in_moleculetype = True
                    molecule_lines = [line]
                    continue

                # 3. Detect end: ; Include water topology (stop before this line)
                if in_moleculetype and stripped == water_marker:
                    found_water_marker = True
                    break

                # 4. Collect molecule definition lines
                if in_moleculetype:
                    molecule_lines.append(line)

            # If we found moleculetype but not water marker, take all remaining lines
            if in_moleculetype and not found_water_marker:
                # Take everything up to [ molecules ] or [ system ]
                final_lines = []
                for ml in molecule_lines:
                    if ml.strip().startswith('[ molecules ]') or ml.strip().startswith('[ system ]'):
                        break
                    final_lines.append(ml)
                molecule_lines = final_lines

            # Store molecule definition
            if molecule_lines:
                molecule_definitions[comp['name']] = molecule_lines

        # Write merged topology
        with open(merged_top, 'w') as f:
            # 1. Write forcefield include
            if forcefield_include:
                f.write(forcefield_include)
            else:
                # Fallback: write default include
                f.write('#include "pace-new.ff/forcefield.itp"\n')
            f.write('\n')

            # 2. Write all molecule definitions
            for comp in component_topologies:
                if comp['name'] in molecule_definitions:
                    f.writelines(molecule_definitions[comp['name']])
                    f.write('\n')

            # 3. Write [ system ] section
            f.write('\n[ system ]\n')
            f.write('; Name\n')
            f.write('MergedSystem\n\n')

            # 4. Write [ molecules ] section
            f.write('[ molecules ]\n')
            f.write('; name\tnumber\n')
            for comp in component_topologies:
                f.write(f"{comp['name']}\t{comp['nmol']}\n")

        return merged_top
    
    def run_pdb2gmx(self, input_pdb: str, output_dir: Path) -> str:
        """Run pdb2gmx on the input structure
        
        Uses the force field specified in self.forcefield_type:
        - PACE-NEW: ff="pace-new"
        - PACE-ASM: ff="pace-asm"
        """
        # Determine force field based on forcefield_type
        if self.forcefield_type == PACE_FORCEFIELD_NEW:
            ff_name = "pace-new"
        elif self.forcefield_type == PACE_FORCEFIELD_ASM:
            ff_name = "pace-asm"
        else:
            raise ValueError(f"Unknown forcefield_type: {self.forcefield_type}")
        
        print(f"  Using force field: {ff_name}")
        
        # Import gromacs wrapper - need to avoid conflict with old_code/gromacs.py
        # Temporarily remove old_code from path to avoid import conflict
        old_code_path = Path(__file__).parent.parent / "old_code"
        old_code_str = str(old_code_path.resolve())
        removed_from_path = False
        
        # Check all occurrences in sys.path (might be added multiple times)
        path_indices_to_remove = []
        for i, path_item in enumerate(sys.path):
            if str(Path(path_item).resolve()) == old_code_str:
                path_indices_to_remove.append(i)
        
        # Remove from end to beginning to maintain indices
        for i in reversed(path_indices_to_remove):
            sys.path.pop(i)
            removed_from_path = True
        
        try:
            # Try importing GromacsWrapper package
            import gromacs
        except ImportError as e:
            # Restore path if we removed it
            if removed_from_path:
                sys.path.insert(0, old_code_str)
            raise ImportError(
                f"GromacsWrapper is not installed. Please install it: pip install GromacsWrapper\n"
                f"Original error: {e}"
            )
        
        # Restore path if we removed it
        if removed_from_path:
            sys.path.insert(0, old_code_str)
        
        gromacs_dir = output_dir / "gromacs"
        gromacs_dir.mkdir(parents=True, exist_ok=True)
        
        input_pdb_path = Path(input_pdb).resolve()
        local_pdb = gromacs_dir / input_pdb_path.name
        output_gro = gromacs_dir / "processed.gro"
        
        # Remove OT2 atoms and rename OT1 to O (PACE force field doesn't support OT2/OXT)
        # Import MDAnalysis for proper atom manipulation
        try:
            import MDAnalysis as mda
        except ImportError:
            raise ImportError("MDAnalysis is required for PDB preprocessing. Install it: pip install MDAnalysis")
        
        # Load universe
        universe = mda.Universe(str(input_pdb_path))
        
        # Find all OT2 atoms
        ot2_atoms = universe.select_atoms("name OT2")
        if len(ot2_atoms) > 0:
            print(f"  Found {len(ot2_atoms)} OT2 atoms to remove")
        
        # Find all OT1 atoms and rename them to 'O'
        ot1_atoms = universe.select_atoms("name OT1")
        if len(ot1_atoms) > 0:
            print(f"  Found {len(ot1_atoms)} OT1 atoms to rename to 'O'")
            for atom in ot1_atoms:
                atom.name = 'O'
        
        # Select all atoms except OT2
        non_ot2_atoms = universe.select_atoms("not name OT2")
        
        # Write cleaned PDB
        with mda.Writer(str(local_pdb), multiframe=False) as writer:
            writer.write(non_ot2_atoms)
        
        # Clean up MDAnalysis objects to avoid memory issues
        del universe, ot2_atoms, ot1_atoms, non_ot2_atoms
        
        print(f"  Cleaned PDB saved to: {local_pdb}")
        
        # Run pdb2gmx with the appropriate force field
        gromacs.pdb2gmx(
            f=str(local_pdb),
            ff=ff_name,  # Dynamically selected based on forcefield_type
            water="no",
            ignh=True,
            o=str(output_gro),
            p=str(gromacs_dir / "topol.top"),
            i=str(gromacs_dir / "posre.itp")
        )
        
        return str(output_gro)
    
    def run_openmm_optimization(self, structure_gro: str, topology_top: str, 
                               output_dir: Path) -> str:
        """
        Run multi-step OpenMM optimization using subprocess.
        
        This method runs the optimization in a separate Python process
        to avoid memory corruption issues that can occur when repeatedly
        creating OpenMM systems in the same process.
        
        Steps:
        1. Gaussian repulsion
        2. Softcore (multiple steps)
        3. Standard force field
        """
        import subprocess
        
        openmm_dir = output_dir / "openmm"
        openmm_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy files to openmm directory
        gro_path = Path(structure_gro)
        top_path = Path(topology_top)
        
        conf_gro = openmm_dir / "conf.gro"
        pace_top = openmm_dir / "PACE.top"
        
        shutil.copy2(str(gro_path), str(conf_gro))
        shutil.copy2(str(top_path), str(pace_top))
        
        # Find the pace_opt_worker.py script
        script_dir = Path(__file__).parent
        pace_opt_script = script_dir / "pace_opt_worker.py"
        
        if not pace_opt_script.exists():
            raise FileNotFoundError(f"Cannot find pace_opt_worker.py at {pace_opt_script}")
        
        # Run optimization in subprocess
        click.echo(f"  Running OpenMM optimization in isolated process...")
        
        result = subprocess.run(
            [
                sys.executable, str(pace_opt_script),
                '--input-gro', str(conf_gro),
                '--input-top', str(pace_top),
                '--output', str(openmm_dir),
                '--device', self.pace_opt_config.platform.lower(),
                '--iter', str(5000),
                '--level', self.pace_opt_config.optimization_mode
            ],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            click.echo(f"Error running OpenMM optimization:")
            click.echo(result.stdout)
            click.echo(result.stderr)
            raise RuntimeError(f"OpenMM optimization failed with exit code {result.returncode}")
        
        # Print the output
        click.echo(result.stdout)
        
        # Find the output file
        final_pdb = openmm_dir / "opti_final.pdb"
        if not final_pdb.exists():
            raise FileNotFoundError(f"Expected output file {final_pdb} not found")
        
        return str(final_pdb)
            # Only set Precision property for CUDA/OpenCL platforms
        """Find GROMACS topology directory"""
        gmxshare_path = os.environ.get('GMXSHARE')
        if gmxshare_path:
            top_path = Path(gmxshare_path) / 'gromacs' / 'top'
            if top_path.is_dir():
                return str(top_path)
        
        gmxdata_path = os.environ.get('GMXDATA')
        if gmxdata_path:
            top_path = Path(gmxdata_path) / 'top'
            if top_path.is_dir():
                return str(top_path)
        
        # Fallback locations
        fallback_candidates = [
            '/usr/share/gromacs/top/',
            '/usr/local/share/gromacs/top/',
            '/opt/gromacs/share/gromacs/top/',
            os.path.expanduser('~/gromacs/share/gromacs/top/'),
        ]
        
        for candidate in fallback_candidates:
            if os.path.exists(candidate):
                return candidate
        
        raise RuntimeError("Could not find GROMACS topology directory. Please set GMXSHARE or GMXDATA.")
    
    def run_gromacs_optimization(self, input_pdb: str, topology_top: str, 
                                 output_dir: Path) -> str:
        """Run GROMACS vacuum optimization (optional)"""
        # This would implement the vacuum optimization from run_solvent.py
        # For now, return the input PDB
        # TODO: Implement GROMACS optimization
        return input_pdb
