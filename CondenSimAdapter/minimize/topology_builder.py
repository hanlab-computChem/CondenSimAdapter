#!/usr/bin/env python3
"""
pdb2gmx Utility Functions for AMBER/CHARMM Force Fields

Provides utility functions for:
- Running PCcli to generate initial structure from sequence
- Running pdb2gmx for topology and structure generation
- Merging multi-component topologies (similar to PACE-opt)
"""

import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Optional

from ..core.config import Component as CGComponent
from ..forcefield.registry import REGISTRY, get_force_field, get_force_field_path

# Pre-fill "no" for any disulfide bond prompts from pdb2gmx.
# Large enough to cover unknown number of prompts.
PDB2GMX_SS_NO_INPUT = "n\n" * 200


def _build_pdb2gmx_input(
    disable_disulfide: bool,
    his_type: Optional[int],
    his_repeat_count: int
) -> Optional[str]:
    """
    Build input string for pdb2gmx interactive prompts.

    Args:
        disable_disulfide: Pre-fill "n" for disulfide prompts (-ss)
        his_type: 0 for HID, 1 for HIE; None to skip -his prompts
        his_repeat_count: Number of histidine selections to pre-fill

    Returns:
        Combined input string or None if no interactive input is needed.
    """
    parts = []
    if disable_disulfide:
        parts.append(PDB2GMX_SS_NO_INPUT)
    if his_type is not None:
        if his_type not in (0, 1):
            raise ValueError(f"his_type must be 0 or 1, got: {his_type}")
        parts.append(f"{his_type}\n" * max(his_repeat_count, 1))
    return "".join(parts) if parts else None


# =============================================================================
# PCcli Functions
# =============================================================================

def run_pccli(sequence: str, output_pdb: Path, protein_name: str) -> bool:
    """
    Run PCcli to generate initial structure from amino acid sequence.
    
    Args:
        sequence: Amino acid sequence
        output_pdb: Output PDB file path
        protein_name: Name of the protein
    
    Returns:
        True if successful, False otherwise
    """
    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        'PCcli',
        '-s', sequence,
        '-o', str(output_pdb),
        '-ss', 'l'  # single letter chain
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    PCcli failed: {result.stderr[-200:]}")
            return False
        
        if output_pdb.exists():
            print(f"    PCcli generated: {output_pdb.name}")
            return True
        else:
            print(f"    PCcli output not found: {output_pdb}")
            return False
    
    except Exception as e:
        print(f"    PCcli error: {e}")
        return False


# =============================================================================
# PDB2GMX Topology Functions
# =============================================================================

def run_pdb2gmx_for_topology(
    input_pdb: Path, 
    output_dir: Path, 
    ff_name: str, 
    molecule_name: str,
    water_model: str = "none",  # Use "none" and manually add water/ions in merged topology
    disable_disulfide: bool = False,
    his_type: Optional[int] = None,
    his_repeat_count: int = 30
) -> Path:
    """
    Run pdb2gmx to generate topology from input PDB.
    
    This is used for the FIRST pdb2gmx call in dual pdb2gmx strategy:
    - Generates topology (topol.top)
    - Generates structure (processed.gro) which is discarded
    
    Args:
        input_pdb: Input PDB file
        output_dir: Output directory for topology files
        ff_name: Force field name (amber99sb-ildn, charmm36-jul2021, etc.)
        molecule_name: Name for the molecule in topology
        water_model: Water model for pdb2gmx (tip3p, tip4p2005s, etc.)
        disable_disulfide: Disable disulfide detection in pdb2gmx (-ss)
        his_type: Histidine type selection for pdb2gmx (-his), 0=HID, 1=HIE
        his_repeat_count: Number of histidine selections to pre-fill
    
    Returns:
        Path to generated topology file (topol.top)
    """
    input_pdb = Path(input_pdb)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_gro = output_dir / "structure.gro"
    final_top = output_dir / "topol.top"
    
    # Clean up previous files
    for f in output_dir.glob('*.top'):
        f.unlink()
    for f in output_dir.glob('*.gro'):
        f.unlink()
    for f in output_dir.glob('*.itp'):
        f.unlink()
    
    # Preprocess PDB: remove OXT, HETATM, USER lines
    cleaned_pdb = output_dir / "cleaned_for_top.pdb"
    with open(input_pdb, 'r') as infile, open(cleaned_pdb, 'w') as outfile:
        for line in infile:
            # Remove OXT/O atom names
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                if atom_name in ('OT1', 'OT2', 'OXT'):
                    continue
                # Rename OT1 to O if it exists
                if atom_name == 'OT1':
                    line = line[:12] + ' O' + line[14:]
            # Remove USER lines
            if line.startswith('USER'):
                continue
            outfile.write(line)
    
    # Run pdb2gmx
    # Ensure force field folder exists in the execution directory for pdb2gmx
    # We'll run from output_dir, so copy force field there
    ff_path = get_force_field_path(ff_name)
    if ff_path:
        ff_folder_name = Path(ff_path).name
        target_ff_path = output_dir / ff_folder_name
        if not target_ff_path.exists():
            shutil.copytree(ff_path, target_ff_path)

    cmd = [
        'gmx', 'pdb2gmx',
        '-f', str(cleaned_pdb),
        '-o', str(output_gro),
        '-ff', ff_name,
        '-water', water_model,
        '-ignh', 'yes',
        '-merge', 'all'
    ]
    if disable_disulfide:
        cmd.append('-ss')
    if his_type is not None:
        cmd.append('-his')

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            input=_build_pdb2gmx_input(disable_disulfide, his_type, his_repeat_count),
            cwd=str(output_dir)
        )

        if result.returncode != 0:
            raise RuntimeError(f"pdb2gmx failed: {result.stderr[-300:]}")

        topol_in_cwd = output_dir / "topol.top"
        if topol_in_cwd.exists() and not final_top.exists():
            shutil.move(str(topol_in_cwd), str(final_top))

        if not final_top.exists():
            raise RuntimeError(f"Topology not generated: {final_top}")

        modify_topology_molecule_name(final_top, molecule_name)
        return final_top

    except Exception as e:
        print(f"    Error in pdb2gmx: {e}")
        raise


def modify_topology_molecule_name(topol_path: Path, molecule_name: str):
    """
    Modify molecule name in [ moleculetype ] section of topology file.
    
    Args:
        topol_path: Path to topology file
        molecule_name: New molecule name
    """
    with open(topol_path, 'r') as f:
        lines = f.readlines()
    
    in_moleculetype = False
    for i, line in enumerate(lines):
        if line.strip().startswith('[ moleculetype ]'):
            in_moleculetype = True
            continue
        if in_moleculetype and line.strip() and not line.startswith(';'):
            # Replace molecule name (first word)
            parts = line.strip().split()
            if len(parts) >= 1:
                parts[0] = molecule_name
                lines[i] = ' '.join(parts) + '\n'
            break
    
    with open(topol_path, 'w') as f:
        f.writelines(lines)
    
    pass  # molecule name updated silently


# =============================================================================
# Multi-Component Topology Generation
# =============================================================================

def generate_all_atom_topology(
    components: List[CGComponent], 
    ff_name: str, 
    output_dir: Path,
    disable_disulfide: bool = False,
    his_type: Optional[int] = None,
    his_repeat_count: int = 30
) -> Tuple[List[Dict], str]:
    """
    Generate all-atom topology for multiple components.
    
    For each component:
    - IDP: Get sequence from FASTA → PCcli → pdb2gmx
    - MDP: Use fpdb directly → pdb2gmx (preserves folded structure)
    
    Args:
        components: List of CGComponent
        ff_name: Force field name
        output_dir: Output directory
        disable_disulfide: Disable disulfide detection in pdb2gmx (-ss)
        his_type: Histidine type selection for pdb2gmx (-his), 0=HID, 1=HIE
        his_repeat_count: Number of histidine selections to pre-fill
    
    Returns:
        Tuple of:
        - List of dicts: [{'name': comp.name, 'topology': top_path, 'nmol': comp.nmol}, ...]
        - Water model name for this force field
    """
    from ..core.config import ComponentType, Component
    
    output_dir = Path(output_dir)
    
    component_topologies = []
    
    # Get water model for this force field (used in merged topology)
    ff_info = get_force_field(ff_name)
    water_model = ff_info.water_model if ff_info else "tip3p"
    
    for comp in components:
        # Determine input PDB based on component type
        if comp.comp_type == ComponentType.IDP:
            # IDP: Generate structure from sequence using PCcli
            sequence = _get_component_sequence(comp)
            print(f"  [{comp.name}] IDP  {len(sequence)} residues")
            
            pccli_pdb = output_dir / f"{comp.name}_pccli.pdb"
            if not run_pccli(sequence, pccli_pdb, comp.name):
                raise RuntimeError(f"PCcli failed for component: {comp.name}")
            
            input_pdb = pccli_pdb
            
        elif comp.comp_type == ComponentType.MDP:
            # MDP: Use user-provided fpdb (folded structure)
            if not comp.pdb_path:
                raise ValueError(f"MDP component '{comp.name}' requires fpdb file")
            
            input_pdb = Path(comp.pdb_path)
            if not input_pdb.exists():
                raise FileNotFoundError(f"fpdb not found: {input_pdb}")
            
            print(f"  [{comp.name}] MDP  {input_pdb.name}")
        
        else:
            raise ValueError(f"Unknown component type: {comp.comp_type}")
        
        # Run pdb2gmx for topology (use "none" for water model)
        comp_top_dir = output_dir / comp.name
        comp_top_path = run_pdb2gmx_for_topology(
            input_pdb,
            comp_top_dir,
            ff_name,
            comp.name,
            water_model="none",
            disable_disulfide=disable_disulfide,
            his_type=his_type,
            his_repeat_count=his_repeat_count
        )
        
        component_topologies.append({
            'name': comp.name,
            'topology': str(comp_top_path),
            'nmol': comp.nmol
        })
    
    return component_topologies, water_model


def _get_component_sequence(comp: CGComponent) -> str:
    """
    Get sequence from component.
    
    Args:
        comp: CGComponent instance
    
    Returns:
        Sequence string
    """
    from ..core.config import ComponentType, Component
    
    if comp.sequence:
        return comp.sequence
    
    if comp.comp_type == ComponentType.IDP:
        if comp.fasta_path:
            return _read_fasta(comp.fasta_path, comp.name)
        else:
            raise ValueError(f"Component '{comp.name}' is IDP but no ffasta file")
    elif comp.comp_type == ComponentType.MDP:
        if comp.pdb_path:
            return _seq_from_pdb(comp.pdb_path)
        else:
            raise ValueError(f"Component '{comp.name}' is MDP but no fpdb file")
    
    return ""


def _read_fasta(fasta_path: str, component_name: str = None) -> str:
    """
    Read sequence from FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        component_name: Component name to select sequence
    
    Returns:
        Sequence string
    """
    from Bio import SeqIO
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    if not records:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")
    
    # If component_name specified, try to match
    if component_name:
        for record in records:
            if record.id == component_name or record.name == component_name:
                return str(record.seq)
    
    # Return first sequence
    return str(records[0].seq)


def _seq_from_pdb(pdb_path: str) -> str:
    """
    Extract sequence from PDB file using MDAnalysis.
    
    Args:
        pdb_path: Path to PDB file
    
    Returns:
        Single-letter amino acid sequence string
    """
    import MDAnalysis as mda
    
    # 3-letter to 1-letter amino acid code mapping
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # Common variants
        'HSD': 'H', 'HSE': 'H', 'HSP': 'H',  # Histidine protonation states
        'CYX': 'C',  # Cysteine in disulfide bridge
    }
    
    # Load PDB with MDAnalysis
    universe = mda.Universe(pdb_path)
    protein = universe.select_atoms("protein")
    
    if len(protein) == 0:
        raise ValueError(f"No protein atoms found in PDB: {pdb_path}")
    
    # Get residues and convert each 3-letter code to 1-letter
    residues = protein.residues
    one_letter_seq = ''.join([aa_map.get(r.resname, 'X') for r in residues])
    
    # Warn about unknown residues
    unknown_residues = set([r.resname for r in residues if r.resname not in aa_map])
    if unknown_residues:
        print(f"    Warning: Unknown residues found in PDB (will be marked as 'X'): {unknown_residues}")
    
    return one_letter_seq


# =============================================================================
# Topology Merging (Based on PACE-opt _merge_topologies)
# =============================================================================

def merge_topologies(
    component_topologies: List[Dict], 
    output_dir: Path,
    ff_name: str,
    water_model: str,
    system_name: str = "MinimizedSystem"
) -> Path:
    """
    Merge multiple component topologies into one.
    
    Topology structure:
    1. #include "forcefield.itp"
    2. [ moleculetype ] sections (from [ moleculetype ] to ; Include Position restraint file)
    3. Water and ions includes (manual, since pdb2gmx -water none)
    4. [ system ]
    5. [ molecules ]
    
    Args:
        component_topologies: List of component topology dicts
            [{'name': comp.name, 'topology': top_path, 'nmol': comp.nmol}, ...]
        output_dir: Output directory for merged topology
        ff_name: Force field name (for water/ions include paths)
        water_model: Water model name (for water include path)
        system_name: System name for [ system ] section
    
    Returns:
        Path to merged topology file
    """
    output_dir = Path(output_dir)
    merged_top = output_dir / "topol.top"
    
    # End marker for moleculetype section (before water/posres includes)
    moleculetype_end_marker = "; Include Position restraint file"
    
    forcefield_include = None
    molecule_definitions = {}
    
    # Extract forcefield include and molecule definitions from each component
    for comp in component_topologies:
        with open(comp['topology'], 'r') as f:
            lines = f.readlines()
        
        found_forcefield_include = False
        in_moleculetype = False
        molecule_lines = []
        found_end_marker = False
        
        for i, line in enumerate(lines):
            stripped = line.strip()
            
            # 1. Find forcefield include
            if not found_forcefield_include and stripped.startswith('#include'):
                forcefield_include = line
                found_forcefield_include = True
                continue
            
            # 2. Detect [ moleculetype ]
            if stripped.startswith('[ moleculetype ]'):
                in_moleculetype = True
                molecule_lines = [line]
                continue
            
            # 3. Detect end marker (; Include Position restraint file)
            if in_moleculetype and stripped == moleculetype_end_marker:
                found_end_marker = True
                break
            
            # 4. Collect molecule definition lines
            if in_moleculetype:
                molecule_lines.append(line)
        
        # Handle case without end marker
        if in_moleculetype and not found_end_marker:
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
            # This shouldn't happen for AMBER/CHARMM, but fallback
            raise RuntimeError("No forcefield include found in component topologies")
        f.write('\n')
        
        # 2. Write all molecule definitions
        for comp in component_topologies:
            if comp['name'] in molecule_definitions:
                f.writelines(molecule_definitions[comp['name']])
                f.write('\n')
        
        # 3. Write water and ions includes (manual, since pdb2gmx uses -water none)
        f.write('\n; Include water topology\n')
        f.write(f'#include "./{ff_name}.ff/{water_model}.itp"\n')
        f.write('\n; Include topology for ions\n')
        f.write(f'#include "./{ff_name}.ff/ions.itp"\n')
        
        # 4. Write [ system ] section
        f.write('\n[ system ]\n')
        f.write('; Name\n')
        f.write(f'{system_name}\n\n')
        
        # 5. Write [ molecules ] section
        f.write('[ molecules ]\n')
        f.write('; name\tnumber\n')
        for comp in component_topologies:
            f.write(f"{comp['name']}\t{comp['nmol']}\n")
    
        pass  # merged silently
    return merged_top


# =============================================================================
# Structure Generation (Second pdb2gmx)
# =============================================================================

def run_pdb2gmx_for_structure(
    input_pdb: Path,
    output_dir: Path,
    ff_name: str,
    water_model: str = "none",
    disable_disulfide: bool = False,
    his_type: Optional[int] = None,
    his_repeat_count: int = 30
) -> Path:
    """
    Run pdb2gmx to generate structure from backmap output PDB.
    
    If input_pdb is a directory (backmap output), automatically detects:
    - final.aa.pdb (preferred)
    - final_aa.pdb
    - *.pdb files in the directory
    
    Note: water_model is set to "none" - water and ions are manually added
    in the merged topology file instead.
    
    Args:
        input_pdb: Input PDB file OR backmap output directory
        output_dir: Output directory
        ff_name: Force field name
        water_model: Water model (set to "none", water/ions added manually)
        disable_disulfide: Disable disulfide detection in pdb2gmx (-ss)
        his_type: Histidine type selection for pdb2gmx (-his), 0=HID, 1=HIE
        his_repeat_count: Number of histidine selections to pre-fill
    
    Returns:
        Path to generated structure file (processed.gro)
    """
    input_path = Path(input_pdb).resolve()
    
    # If input is a directory (backmap output), find the PDB file
    if input_path.is_dir():
        print(f"    Input is directory: {input_path}")
        
        # Try common backmap output names
        candidates = [
            input_path / "final.aa.pdb",
            input_path / "final_aa.pdb",
            input_path / "final.pdb",
        ]
        
        # Check candidates
        for candidate in candidates:
            if candidate.exists():
                print(f"    Found PDB file: {candidate.name}")
                input_path = candidate
                break
        else:
            # Fall back to any .pdb file
            pdb_files = list(input_path.glob("*.pdb"))
            if pdb_files:
                # Sort and take the first one (usually final.aa.pdb comes first alphabetically)
                pdb_files.sort()
                input_path = pdb_files[0]
                print(f"    Found PDB file: {input_path.name}")
            else:
                raise FileNotFoundError(
                    f"No PDB file found in directory: {input_pdb}\n"
                    f"Expected files: final.aa.pdb, final_aa.pdb, or any *.pdb"
                )
    elif not input_path.is_file():
        raise FileNotFoundError(f"Input path does not exist: {input_pdb}")
    
    # input_path is now definitely a file
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Clean up previous files
    for f in output_dir.glob('*.gro'):
        f.unlink()
    
    # Run pdb2gmx directly on the PDB file
    # Ensure force field folder exists in the execution directory for pdb2gmx
    pdb2gmx_cwd = output_dir  # Run from structure/ directory
    ff_path = get_force_field_path(ff_name)
    if ff_path:
        ff_folder_name = Path(ff_path).name
        target_ff_path = pdb2gmx_cwd / ff_folder_name
        if not target_ff_path.exists():
            shutil.copytree(ff_path, target_ff_path)

    output_gro = output_dir / "processed.gro"
    cmd = [
        'gmx', 'pdb2gmx',
        '-f', str(input_path),
        '-o', str(output_gro),
        '-ff', ff_name,
        '-water', water_model,
        '-ignh', 'yes',
        '-merge', 'all'
    ]
    if disable_disulfide:
        cmd.append('-ss')
    if his_type is not None:
        cmd.append('-his')
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            input=_build_pdb2gmx_input(disable_disulfide, his_type, his_repeat_count),
            cwd=str(pdb2gmx_cwd)
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"pdb2gmx failed: {result.stderr[-300:]}")
        
        if not output_gro.exists():
            raise RuntimeError(f"Structure file not found: {output_gro}")
        
        return output_gro

    except Exception as e:
        print(f"    Error in pdb2gmx: {e}")
        raise


# =============================================================================
# Component Loading from YAML
# =============================================================================

def load_config_from_yaml(yaml_path: str) -> Tuple[str, List[CGComponent]]:
    """
    Load system_name and components from YAML configuration file.
    
    Args:
        yaml_path: Path to YAML file
    
    Returns:
        Tuple of (system_name, List of CGComponent)
    """
    import yaml
    import os
    
    yaml_path = Path(yaml_path)
    config_dir = yaml_path.parent
    
    with open(yaml_path, 'r') as f:
        d = yaml.safe_load(f)
    
    # Get system_name from config, fallback to first component name or default
    system_name = d.get('system_name', '')
    
    components = []
    
    # Get components section
    components_data = d.get('components', [])
    if not components_data:
        # Try nested structure
        for key in d:
            if isinstance(d[key], dict) and 'components' in d[key]:
                components_data = d[key]['components']
                break
    
    for comp_data in components_data:
        # Handle relative paths
        if comp_data.get('ffasta'):
            ffasta = comp_data['ffasta']
            if not os.path.isabs(ffasta):
                ffasta = str(config_dir / ffasta)
            comp_data['ffasta'] = ffasta
        
        if comp_data.get('fpdb'):
            fpdb = comp_data['fpdb']
            if not os.path.isabs(fpdb):
                fpdb = str(config_dir / fpdb)
            comp_data['fpdb'] = fpdb
        
        comp = CGComponent.from_dict(comp_data)
        components.append(comp)
    
    # If system_name not specified, use first component name or default
    if not system_name and components:
        system_name = components[0].name
    
    if not system_name:
        system_name = "MinimizedSystem"
    
    return system_name, components


def load_components_from_yaml(yaml_path: str) -> List[CGComponent]:
    """
    Load components from YAML configuration file.
    
    Args:
        yaml_path: Path to YAML file
    
    Returns:
        List of CGComponent
    """
    _, components = load_config_from_yaml(yaml_path)
    return components


# =============================================================================
# Utility Functions
# =============================================================================

def count_chains_in_pdb(pdb_path: Path) -> int:
    """
    Count number of chains in PDB file.
    
    Args:
        pdb_path: Path to PDB file
    
    Returns:
        Number of unique chains
    """
    chains = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21:22].strip()
                if chain_id:
                    chains.add(chain_id)
    return len(chains)


def verify_files_valid(top_path: str, gro_path: str) -> Tuple[bool, str]:
    """
    Verify that both topology and structure files are valid.
    
    Args:
        top_path: Path to topology file
        gro_path: Path to structure file
    
    Returns:
        Tuple of (success: bool, message: str)
    """
    from openmm.app import GromacsGroFile
    
    try:
        # Check files exist
        if not Path(top_path).exists():
            return False, f"Topology file not found: {top_path}"
        if not Path(gro_path).exists():
            return False, f"Structure file not found: {gro_path}"
        
        # Try to load GRO file
        gro = GromacsGroFile(gro_path)
        n_atoms_gro = len(gro.getPositions())
        
        # Count atoms in TOP file (basic check)
        with open(top_path, 'r') as f:
            content = f.read()
            has_atoms_section = '[ atoms ]' in content
            has_molecules = '[ molecules ]' in content
            if not has_atoms_section:
                return False, "TOP file missing [ atoms ] section"
            if not has_molecules:
                return False, "TOP file missing [ molecules ] section"
        
        return True, f"TOP valid, GRO has {n_atoms_gro} atoms"
    
    except Exception as e:
        return False, str(e)
