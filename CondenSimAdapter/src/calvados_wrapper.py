#!/usr/bin/env python3
"""
CALVADOS Wrapper

Convert CGSimulationConfig to ms2_calvados format and run simulation.

Topology mapping:
- CUBIC → CALVADOS 'grid'
- SLAB → CALVADOS 'slab'

Usage:
    from CondenSimAdapter.src import CGSimulationConfig
    from CondenSimAdapter.src.calvados_wrapper import run_calvados
    
    config = CGSimulationConfig.from_yaml("config.yaml")
    result = run_calvados(config, output_dir="output/", gpu_id=0)
"""

import os
from pathlib import Path
from typing import Dict, Optional

from .cg import (
    CGSimulationConfig,
    CGComponent,
    ComponentType,
    TopologyType,
    ComputePlatform,
    SimulationResult,
)


class CalvadosWrapper:
    """
    CALVADOS simulation wrapper
    
    Convert CGSimulationConfig to ms2_calvados format and run simulation.
    
    Attributes:
        config: CGSimulationConfig instance
        output_dir: Output directory
        ms2_config: ms2_calvados Config object
        ms2_components: ms2_calvados Components object
    """
    
    def __init__(self, config: CGSimulationConfig):
        """
        Initialize wrapper
        
        Args:
            config: CGSimulationConfig instance
        """
        self.config = config
        self.output_dir: Optional[str] = None
        
        # Get residue file path
        self._residues_path = self._get_residues_path()
    
    def _get_residues_path(self) -> str:
        """Get residue parameter file path
        
        Residues files are loaded from ms2_calvados package's data directory:
        - residues_CALVADOS2.csv: For pure IDP systems
        - residues_CALVADOS3.csv: For systems containing MDP
        """
        from CondenSimAdapter.extern.ms2_calvados.calvados import data as calvados_data
        
        has_mdp = any(c.type == ComponentType.MDP for c in self.config.components)
        residues_file = 'residues_CALVADOS3.csv' if has_mdp else 'residues_CALVADOS2.csv'
        
        # calvados_data is a namespace package, use __path__ instead of __file__
        data_path = calvados_data.__path__[0]
        residues_path = Path(data_path) / residues_file
        
        if not residues_path.exists():
            raise FileNotFoundError(f"Residues file not found: {residues_path}")
        
        return str(residues_path)
    
    def _topol_to_calvados(self) -> str:
        """Convert TopologyType to CALVADOS topology string"""
        if self.config.topol == TopologyType.CUBIC:
            return 'grid'
        elif self.config.topol == TopologyType.SLAB:
            return 'slab'
        elif self.config.topol == TopologyType.DROPLET:
            return 'grid'
        else:
            return 'grid'  # Default

    def _is_droplet_topology(self) -> bool:
        if isinstance(self.config.topol, TopologyType):
            return self.config.topol == TopologyType.DROPLET
        return str(self.config.topol).lower() == 'droplet'

    def _get_droplet_params(self) -> tuple:
        box = self.config.box
        center = (box[0] / 2.0, box[1] / 2.0, box[2] / 2.0)
        radius = self.config.droplet_radius if self.config.droplet_radius is not None else box[0] / 2.0
        return radius, center

    def _get_droplet_force_expr(self) -> str:
        radius, center = self._get_droplet_params()
        k = 1.0
        return (
            f"{k} * step(r - {radius}) * (r - {radius}); "
            f"r = sqrt((x-{center[0]})^2 + (y-{center[1]})^2 + (z-{center[2]})^2)"
        )
    
    def _platform_to_string(self) -> str:
        """Convert ComputePlatform to string"""
        if isinstance(self.config.simulation.platform, ComputePlatform):
            return self.config.simulation.platform.value
        return str(self.config.simulation.platform)
    
    def create_config(self) -> 'ms2_config.Config':
        """Create ms2_calvados Config object
        
        Notes:
            - Only pass parameters that users need to modify, let Config class use default_config.yaml defaults
            - CALVADOS physical constants (eps_lj, cutoff_lj, friction_coeff, etc.) remain unchanged
            - slab_width: Automatically calculated as box[2] / 2 for SLAB topology
        """
        from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Config
        
        sim_params = self.config.simulation
        
        # SLAB topology: automatically calculate slab_width = box_z / 2
        # Other topologies: use CALVADOS default (100)
        if self.config.topol.value == 'slab':
            slab_width = self.config.box[2] / 2
        else:
            slab_width = None  # Use CALVADOS default
        
        # Only pass parameters that users actually configured
        params = {
            'sysname': self.config.system_name,
            'box': self.config.box,
            'temp': self.config.temperature,
            'ionic': self.config.ionic,
            'pH': 7.0,
            'topol': self._topol_to_calvados(),
            'wfreq': sim_params.wfreq,
            'steps': sim_params.steps,
            'platform': self._platform_to_string(),
            'verbose': sim_params.verbose,
        }
        
        # SLAB topology needs to specify slab_width
        if slab_width is not None:
            params['slab_width'] = slab_width
        
        if self._is_droplet_topology():
            params['ext_force'] = True
            params['ext_force_expr'] = self._get_droplet_force_expr()
        
        return Config(**params)
    
    def create_components(self) -> 'ms2_config.Components':
        """Create ms2_calvados Components object"""
        from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Components
        
        first_comp = self.config.components[0] if self.config.components else None
        
        defaults = {
            'molecule_type': 'protein',
            'nmol': first_comp.nmol if first_comp else 1,
            'restraint': first_comp.restraint if first_comp else False,
            'charge_termini': first_comp.charge_termini if first_comp else 'both',
            'fresidues': self._residues_path,
        }
        
        ms2_components = Components(**defaults)
        
        for comp in self.config.components:
            comp_dict = {
                'name': comp.name,
                'nmol': comp.nmol,
                'restraint': comp.restraint,
                'charge_termini': comp.charge_termini,
                'fresidues': self._residues_path,
            }
            
            if comp.type == ComponentType.IDP:
                if comp.ffasta:
                    comp_dict['ffasta'] = comp.ffasta
            
            elif comp.type == ComponentType.MDP:
                if comp.fpdb:
                    comp_dict['fpdb'] = comp.fpdb
                    comp_dict['pdb_folder'] = os.path.dirname(os.path.abspath(comp.fpdb))
                
                if comp.fdomains:
                    comp_dict['fdomains'] = comp.fdomains
                
                if comp.restraint:
                    comp_dict['restraint_type'] = comp.restraint_type
                    comp_dict['use_com'] = comp.use_com
                    comp_dict['k_harmonic'] = comp.k_harmonic
                    comp_dict['colabfold'] = comp.colabfold
            
            ms2_components.add(**comp_dict)
        
        return ms2_components
    
    def write(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
        """
        Write configuration files
        
        Args:
            output_dir: Output directory
            overwrite: Whether to overwrite
            
        Returns:
            Dictionary of generated file paths
        """
        output_dir = os.path.abspath(output_dir)
        
        if os.path.exists(output_dir) and not overwrite:
            raise FileExistsError(f"Output directory exists: {output_dir}")
        
        os.makedirs(output_dir, exist_ok=True)
        self.output_dir = output_dir
        
        # Create and write config
        ms2_config = self.create_config()
        ms2_config.write(output_dir, name='config.yaml')
        
        # Create and write components
        ms2_components = self.create_components()
        ms2_components.write(output_dir, name='components.yaml')
        
        return {
            'config': os.path.join(output_dir, 'config.yaml'),
            'components': os.path.join(output_dir, 'components.yaml'),
            'run_script': os.path.join(output_dir, 'run.py'),
        }
    
    def _generate_config_yaml(self, gpu_id: int = 0, verbose: bool = False, continue_from: str = None) -> str:
        """Generate CALVADOS config.yaml content

        Strategy:
        1. Load CALVADOS default_config.yaml as base configuration
        2. Only override parameters that users actually configured
        3. This avoids hardcoding all physical constants (eps_lj, cutoff_lj, etc.)

        This approach maintains consistent design philosophy with the original CALVADOS Config class.

        Args:
            gpu_id: GPU device ID (user specified GPU)
            verbose: Whether to output detailed logs
            continue_from: PDB file path for continuing simulation (PDB format)
        """
        import yaml
        from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Config
        
        sim_params = self.config.simulation
        
        # SLAB topology: auto calculate slab_width = box_z / 2
        # Other topologies: use CALVADOS default
        if self.config.topol == TopologyType.SLAB:
            slab_width = self.config.box[2] / 2
        else:
            slab_width = None
        
        # Use Config class to load default configuration
        config_obj = self.create_config()
        config_dict = config_obj.config.copy()
        
        # Override user configured parameters (including gpu_id and verbose)
        config_dict.update({
            'sysname': self.config.system_name,
            'box': self.config.box,
            'temp': self.config.temperature,
            'ionic': self.config.ionic,
            'pH': 7.0,
            'topol': self._topol_to_calvados(),
            'wfreq': sim_params.wfreq,
            'steps': sim_params.steps,
            'platform': self._platform_to_string(),
            'verbose': verbose,  # Control CALVADOS detailed output
            'gpu_id': gpu_id,  # User specified GPU ID
        })
        
        # If continue_from is provided, set restart='pdb' and frestart
        if continue_from:
            config_dict.update({
                'restart': 'pdb',
                'frestart': continue_from,
            })
        
        # SLAB topology needs to specify slab_width
        if slab_width is not None:
            config_dict['slab_width'] = slab_width
        
        if self._is_droplet_topology():
            config_dict['topol'] = 'grid'
            config_dict['ext_force'] = True
            config_dict['ext_force_expr'] = self._get_droplet_force_expr()
        
        return yaml.dump(config_dict, default_flow_style=False, sort_keys=False)
    
    def _generate_components_yaml(self) -> str:
        """Generate CALVADOS components.yaml content
        
        Handle fpdb and pdb_folder:
        - CALVADOS expects pdb_folder (directory) and name (filename without extension)
        - Our config uses fpdb (full file path)
        
        Notes:
            Added all default parameters from original CALVADOS default_component.yaml:
            - periodic: false
            - cutoff_restr: 0.9
            - k_go: 15.
            - use_com: true
            - colabfold: 0
        """
        import yaml
        
        first_comp = self.config.components[0] if self.config.components else None
        
        # Calculate pdb_folder (extracted from first MDP component's fpdb)
        pdb_folder = None
        for comp in self.config.components:
            if comp.type.value == 'mdp' and comp.fpdb:
                pdb_folder = os.path.dirname(os.path.abspath(comp.fpdb))
                break
        
        components = {
            'defaults': {
                'molecule_type': 'protein',
                'nmol': first_comp.nmol if first_comp else 1,
                'restraint': first_comp.restraint if first_comp else False,
                'charge_termini': first_comp.charge_termini if first_comp else 'both',
                'fresidues': self._residues_path,
                'alpha': 0,
                'kb': 8033.0,
                'pdb_folder': pdb_folder,
                # Parameters from original CALVADOS default_component.yaml
                'periodic': False,
                'cutoff_restr': 0.9,
                'k_go': 15.0,
                'use_com': True,
                'colabfold': 0,
            },
            'system': {}
        }
        
        for comp in self.config.components:
            # For MDP, verify fpdb file exists
            # CALVADOS expects file path as {pdb_folder}/{name}.pdb
            if comp.type == ComponentType.MDP and comp.fpdb:
                pdb_folder = os.path.dirname(os.path.abspath(comp.fpdb))
                expected_pdb = os.path.join(pdb_folder, f"{comp.name}.pdb")
                actual_pdb = os.path.abspath(comp.fpdb)
                
                # Verify file exists
                if not os.path.exists(actual_pdb):
                    raise FileNotFoundError(
                        f"PDB file not found: {actual_pdb}\n"
                        f"  Component: {comp.name}\n"
                        f"  Expected by CALVADOS: {expected_pdb}\n"
                        f"\n"
                        f"Choose one of the following solutions:\n"
                        f"  1. Rename PDB file: mv h1.pdb H1.pdb\n"
                        f"  2. Modify component name: name: h1  (lowercase)"
                    )
                
                # Warn if actual file path doesn't match CALVADOS expectation
                if actual_pdb != expected_pdb:
                    print(f"\n  Warning: Component '{comp.name}' has fpdb='{comp.fpdb}'")
                    print(f"      CALVADOS expects file named: {comp.name}.pdb")
                    print(f"      This mismatch WILL cause errors!")
                    print(f"\n  Solution:")
                    print(f"    mv {comp.fpdb} {os.path.join(pdb_folder, comp.name)}.pdb")
                    print(f"")
            
            comp_dict = {
                'name': comp.name,
                'molecule_type': 'protein',
                'nmol': comp.nmol,
                'ffasta': comp.ffasta,
                'fdomains': comp.fdomains,
                'fpdb': comp.fpdb,
                'restraint': comp.restraint,
                'restraint_type': comp.restraint_type,
                'use_com': comp.use_com,
                'k_harmonic': comp.k_harmonic,
                'colabfold': comp.colabfold,
                'charge_termini': comp.charge_termini,
            }
            # Remove None values
            comp_dict = {k: v for k, v in comp_dict.items() if v is not None}
            components['system'][comp.name] = comp_dict
        
        return yaml.dump(components, default_flow_style=False, sort_keys=False)
    
    def _write_to_dir(self, output_dir: str, gpu_id: int = 0, verbose: bool = False, continue_from: str = None) -> Dict[str, str]:
        """Write configuration files to specified directory (return file path dictionary)

        Support two fdomains formats:
        1. File path: 'TDP43_domains.yaml' - directly copy to output directory
        2. Inline YAML: 'TDP43:\n  - [3, 76]\n...' - write to temp file

        Args:
            output_dir: Output directory
            gpu_id: GPU device ID (for writing config.yaml)
            verbose: Whether to output detailed logs
            continue_from: PDB file path for continuing simulation (will be copied to output directory)
        """
        import tempfile
        import shutil

        os.makedirs(output_dir, exist_ok=True)
        self.output_dir = output_dir

        # Write config.yaml (pass gpu_id, verbose, and continue_from)
        config_file = os.path.join(output_dir, 'config.yaml')
        with open(config_file, 'w') as f:
            f.write(self._generate_config_yaml(gpu_id=gpu_id, verbose=verbose, continue_from=continue_from))
        
        # Handle components.yaml, support inline fdomains
        components_yaml = self._generate_components_yaml()
        
        # Check if there are inline fdomains to process
        components_yaml = self._process_inline_fdomains(components_yaml, output_dir)
        
        # Write components.yaml
        components_file = os.path.join(output_dir, 'components.yaml')
        with open(components_file, 'w') as f:
            f.write(components_yaml)
        
        return {
            'config': config_file,
            'components': components_file,
        }
    
    def _process_inline_fdomains(self, components_yaml: str, output_dir: str) -> str:
        """Process inline fdomains, write to temp file if YAML content"""
        import yaml

        # Parse YAML
        components = yaml.safe_load(components_yaml)

        for name, props in components.get('system', {}).items():
            fdomains = props.get('fdomains')
            if fdomains and isinstance(fdomains, str):
                # Remove YAML quotes (single or double quotes)
                stripped = fdomains.strip()
                if stripped.startswith('"') and stripped.endswith('"'):
                    stripped = stripped[1:-1]
                elif stripped.startswith("'") and stripped.endswith("'"):
                    stripped = stripped[1:-1]

                # Check if it's inline YAML (not a file path)
                is_inline = False
                if stripped.startswith('{') or stripped.startswith('['):
                    is_inline = True
                elif '\n' in stripped and (':' in stripped or stripped.startswith('-')):
                    # Multi-line content with YAML features
                    is_inline = True
                elif ':' in stripped and not stripped.endswith('.yaml') and not stripped.endswith('.yml'):
                    # Contains colon but doesn't look like file path
                    is_inline = True

                if is_inline:
                    try:
                        # Try to parse as YAML
                        domains_data = yaml.safe_load(stripped)

                        # Ensure parse result is a dictionary
                        if isinstance(domains_data, dict):
                            # Only write current protein's domain data, use protein name as key
                            protein_domains = {name: domains_data.get(name, [])}
                        elif isinstance(domains_data, list):
                            # Direct domain list [[3, 76], ...]
                            protein_domains = {name: domains_data}
                        else:
                            continue

                        # Write to temp file
                        domains_file = os.path.join(output_dir, f'{name}_domains.yaml')
                        with open(domains_file, 'w') as f:
                            yaml.dump(protein_domains, f, default_flow_style=False)

                        # Replace with file path
                        props['fdomains'] = domains_file

                    except yaml.YAMLError:
                        # Not valid YAML, keep as is (might be file path)
                        pass

        return yaml.dump(components, default_flow_style=False, sort_keys=False)
    
    def run(self, output_dir: str = None, gpu_id: int = 0, verbose: bool = False) -> SimulationResult:
        """
        Run CALVADOS simulation
        
        Unified output structure:
        {system_name}_CG/
        ├── final.pdb                   # Final structure
        ├── trajectory.xtc              # Simulation trajectory
        ├── simulation.log              # High-level log
        └── raw/                        # Native output
            ├── config.yaml
            ├── components.yaml
            ├── *.xtc, *.xml, *.pdb, *.chk, *.txt
        
        Args:
            output_dir: Output directory (use config's output_dir by default, use directly if provided)
            gpu_id: GPU device ID
            
        Returns:
            SimulationResult
        """
        from CondenSimAdapter.extern.ms2_calvados.calvados import sim as calvados_sim
        import shutil
        from datetime import datetime
        import time
        
        if output_dir is None:
            output_dir = self.config.output_dir
        
        # Unified _CG suffix addition
        task_name = f"{self.config.system_name}_CG"
        output_dir = os.path.join(output_dir, task_name)
        raw_dir = os.path.join(output_dir, 'raw')
        
        # If directory exists, backup and recreate
        if os.path.exists(output_dir):
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{output_dir}_backup_{timestamp}"
            shutil.move(output_dir, backup_dir)
            print(f"  Backup old results to: {backup_dir}")
        
        os.makedirs(raw_dir, exist_ok=True)

        # Write configuration files to raw directory (pass gpu_id and verbose)
        files = self._write_to_dir(raw_dir, gpu_id=gpu_id, verbose=verbose)
        
        result = SimulationResult()
        result.output_dir = output_dir
        
        start_time = time.time()
        
        try:
            print(f"\n[CALVADOS] Running simulation...")
            print(f"  GPU ID: {gpu_id}")
            print(f"  Task: {task_name}")
            print(f"  Raw output: {raw_dir}")
            print(f"  Topology: {self._topol_to_calvados()}")
            
            # Run simulation (output to raw directory)
            calvados_sim.run(
                path=raw_dir,
                fconfig='config.yaml',
                fcomponents='components.yaml'
            )
            
            # Organize output files
            self._organize_output(raw_dir, output_dir, task_name)
            
            result.success = True
            elapsed = time.time() - start_time
            print(f"  CALVADOS simulation completed ({elapsed:.1f}s)")
            
        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            print(f"  CALVADOS simulation failed: {e}")
            elapsed = time.time() - start_time
        
        # Write high-level simulation log
        self._write_simulation_log(output_dir, task_name, elapsed, result.success)
        
        # Set result file paths
        result.trajectory = os.path.join(output_dir, 'trajectory.xtc')
        result.structure = os.path.join(output_dir, 'final.pdb')
        
        for key in ['trajectory', 'structure']:
            path = getattr(result, key)
            if path and not os.path.exists(path):
                setattr(result, key, None)
        
        return result
    
    def _organize_output(self, raw_dir: str, output_dir: str, task_name: str):
        """
        Organize output files to unified structure
        
        Unified naming rules:
        - trajectory.xtc  <- {task_name}.xtc
        - final.pdb       <- Timestamped pdb or checkpoint.pdb
        """
        import shutil
        
        sysname = self.config.system_name
        
        # 1. Process trajectory file
        src_xtc = os.path.join(raw_dir, f'{sysname}.xtc')
        dst_xtc = os.path.join(output_dir, 'trajectory.xtc')
        if os.path.exists(src_xtc):
            shutil.copy2(src_xtc, dst_xtc)
            print(f"  trajectory.xtc")
        
        # 2. Find and copy final structure (prefer checkpoint.pdb, otherwise look for timestamped PDB)
        src_pdb = os.path.join(raw_dir, 'checkpoint.pdb')
        if not os.path.exists(src_pdb):
            # Look for timestamped PDB
            for f in os.listdir(raw_dir):
                if f.endswith('.pdb') and f != 'top.pdb':
                    src_pdb = os.path.join(raw_dir, f)
                    break
        
        dst_pdb = os.path.join(output_dir, 'final.pdb')
        if os.path.exists(src_pdb):
            shutil.copy2(src_pdb, dst_pdb)
            print(f"  final.pdb")
        
        # 3. Copy important files to raw directory (if not already there)
        important_files = [
            (f'{sysname}.xml', 'system.xml'),
            ('top.pdb', 'top.pdb'),
            ('restart.chk', 'restart.chk'),
            ('checkpoint.pdb', 'checkpoint.pdb'),
        ]
        
        for src_name, dst_name in important_files:
            src = os.path.join(raw_dir, src_name)
            if os.path.exists(src):
                dst = os.path.join(raw_dir, dst_name)
                if src != dst:
                    shutil.copy2(src, dst)
        
        # 4. Rename log files
        for f in os.listdir(raw_dir):
            if f.endswith('.log') or f.endswith('.txt'):
                pass  # Keep original
        
        print(f"  Raw output organized to: {raw_dir}")
    
    def _write_simulation_log(self, output_dir: str, task_name: str, elapsed: float, success: bool):
        """Write high-level simulation log"""
        from datetime import datetime
        
        log_file = os.path.join(output_dir, 'simulation.log')
        
        status = "SUCCESS" if success else "FAILED"
        components_info = []
        for comp in self.config.components:
            comp_info = f"  - {comp.name}: {comp.type.value}, nmol={comp.nmol}"
            if comp.type == ComponentType.IDP:
                comp_info += f", seq={comp.ffasta}"
            elif comp.type == ComponentType.MDP:
                comp_info += f", pdb={comp.fpdb}"
            components_info.append(comp_info)
        
        log_content = f"""# CondenSimAdapter CG Simulation Log
# ============================

Task: {task_name}
Force Field: CALVADOS
Date: {datetime.now().isoformat()}

Status: {status}
Duration: {elapsed:.2f} seconds

System Configuration:
  Box: {self.config.box} nm
  Temperature: {self.config.temperature} K
  Ionic Strength: {self.config.ionic} M
  Topology: {self.config.topol.value}

Components ({len(self.config.components)}):
{chr(10).join(components_info)}

Output Files:
  - final.pdb: Final structure
  - trajectory.xtc: Simulation trajectory
  - raw/: Native simulation output files
"""
        with open(log_file, 'w') as f:
            f.write(log_content)
        
        print(f"  simulation.log")


def run_calvados(config: CGSimulationConfig, output_dir: str = None, gpu_id: int = 0) -> SimulationResult:
    """
    Convenient function to run CALVADOS simulation
    
    Args:
        config: CGSimulationConfig instance
        output_dir: Output directory
        gpu_id: GPU device ID
        
    Returns:
        SimulationResult
    """
    wrapper = CalvadosWrapper(config)
    return wrapper.run(output_dir=output_dir, gpu_id=gpu_id)


__all__ = [
    'CalvadosWrapper',
    'run_calvados',
]

