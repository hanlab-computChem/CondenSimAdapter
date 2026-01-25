#!/usr/bin/env python3
"""
PLUMED Generator for Minimize Post-Processing

Automatically generates plumed.dat for MDP components with contact map restraints.
Uses user-provided reference structures (fpdb) and global atom indices from topology generation.
"""

import os
import yaml
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance

from .cg import CGComponent, ComponentType


# =============================================================================
# Main Function
# =============================================================================

def generate_plumed_for_minimize(
    components: List[CGComponent],
    topology_dir: Path,
    output_file: str,
    cutoff: float = 4.5,  # Angstrom
    kappa: float = 10000.0,
    verbose: bool = True
) -> bool:
    """
    为 minimize 结果生成 plumed.dat
    
    Args:
        components: CGComponent 列表
        topology_dir: Topology 生成目录（包含每个 component 的 structure.gro）
        output_file: 输出 plumed.dat 路径
        cutoff: Contact map 距离阈值（Angstrom）
        kappa: Restraint 力常数
        verbose: 是否输出详细信息
        
    Returns:
        是否成功生成（True if MDP components exist）
    """
    # 1. Check for MDP components
    mdp_components = [comp for comp in components if comp.type == ComponentType.MDP]
    
    if not mdp_components:
        if verbose:
            print("  No MDP components found, skipping plumed.dat generation")
        return False
    
    if verbose:
        print(f"  Found {len(mdp_components)} MDP component(s)")
    
    # 2. Build global index mapping
    try:
        global_index_map = _build_global_index_map(topology_dir, components)
    except Exception as e:
        if verbose:
            print(f"  Error building global index map: {e}")
        return False
    
    # 3. Generate contact maps for each MDP component
    all_contactmaps = []  # List of (comp, copy_id, domain_id, domain_range, pairs)
    
    for comp in mdp_components:
        # Skip if no fdomains
        if not comp.fdomains:
            if verbose:
                print(f"  Warning: MDP component '{comp.name}' has no fdomains, skipping")
            continue
        
        # Skip if no fpdb (shouldn't happen for MDP)
        if not comp.fpdb:
            if verbose:
                print(f"  Warning: MDP component '{comp.name}' has no fpdb, skipping")
            continue
        
        # Parse fdomains
        try:
            domains = _parse_fdomains(comp.fdomains)
        except Exception as e:
            if verbose:
                print(f"  Error parsing fdomains for '{comp.name}': {e}")
            continue
        
        if not domains:
            if verbose:
                print(f"  Warning: No domains found for '{comp.name}', skipping")
            continue
        
        if verbose:
            print(f"  Component '{comp.name}': {len(domains)} domain(s), {comp.nmol} copy(ies)")
        
        # Use structure.gro from topology generation (contains correct H atoms and ordering)
        ref_structure_path = topology_dir / comp.name / "structure.gro"
        if not ref_structure_path.exists():
            if verbose:
                print(f"  Warning: structure.gro not found for '{comp.name}', skipping")
            continue
        
        # Generate contact maps for each domain
        domain_contact_pairs = []  # List of contact pairs for each domain (local indices)
        
        for domain_id, (domain_start, domain_end) in enumerate(domains):
            try:
                pairs = _generate_contactmap_for_domain(
                    str(ref_structure_path), domain_start, domain_end, cutoff
                )
                domain_contact_pairs.append(pairs)
                
                if verbose:
                    print(f"    Domain {domain_id}: residues {domain_start}-{domain_end}, {len(pairs)} contact pairs")
            except Exception as e:
                if verbose:
                    print(f"    Error generating contact map for domain {domain_id}: {e}")
                domain_contact_pairs.append([])
        
        # Convert to global indices for each copy
        for copy_id in range(comp.nmol):
            start_atom_idx, end_atom_idx = global_index_map[(comp.name, copy_id)]
            
            for domain_id, pairs in enumerate(domain_contact_pairs):
                if len(pairs) == 0:
                    continue
                
                # Convert local indices to global indices
                # pairs are (atom1_local, atom2_local, distance) with 1-based indexing
                global_pairs = []
                for atom1_local, atom2_local, dist in pairs:
                    # atom1_local and atom2_local are 1-based (for PLUMED)
                    # Convert to global: add (start_atom_idx - 1) to get 0-based, then +1 for PLUMED
                    atom1_global = atom1_local + (start_atom_idx - 1)
                    atom2_global = atom2_local + (start_atom_idx - 1)
                    global_pairs.append((atom1_global, atom2_global, dist))
                
                all_contactmaps.append({
                    'component': comp.name,
                    'copy_id': copy_id,
                    'domain_id': domain_id,
                    'domain_range': (domain_start, domain_end),
                    'pairs': global_pairs,
                    'ref_structure': str(ref_structure_path)
                })
    
    if not all_contactmaps:
        if verbose:
            print("  No contact maps generated, skipping plumed.dat")
        return False
    
    # 4. Write plumed.dat
    try:
        _write_plumed_dat(all_contactmaps, output_file, kappa, verbose)
    except Exception as e:
        if verbose:
            print(f"  Error writing plumed.dat: {e}")
        return False
    
    if verbose:
        print(f"  Successfully generated plumed.dat with {len(all_contactmaps)} CONTACTMAP(s)")
    
    return True


# =============================================================================
# Helper Functions
# =============================================================================

def _build_global_index_map(
    topology_dir: Path,
    components: List[CGComponent]
) -> Dict[Tuple[str, int], Tuple[int, int]]:
    """
    构建 (component_name, copy_id) -> (start_idx, end_idx) 映射
    
    从每个 component 的 structure.gro 读取原子数，按顺序计算偏移量。
    
    Args:
        topology_dir: topology 生成目录
        components: CGComponent 列表
        
    返回格式：
    {
        ('TDP43', 0): (1, 3138),
        ('TDP43', 1): (3139, 6276),
        ('FUS', 0): (6277, 9000),
        ...
    }
    """
    global_index_map = {}
    current_atom_idx = 1  # PLUMED uses 1-based indexing
    
    for comp in components:
        # Read structure.gro to get atom count
        structure_gro = topology_dir / comp.name / "structure.gro"
        
        if not structure_gro.exists():
            raise FileNotFoundError(
                f"Structure file not found for component '{comp.name}': {structure_gro}"
            )
        
        # Use MDAnalysis to read atom count
        u = mda.Universe(str(structure_gro))
        num_atoms_per_monomer = len(u.atoms)
        
        # Create mapping for each copy
        for copy_id in range(comp.nmol):
            start_idx = current_atom_idx
            end_idx = current_atom_idx + num_atoms_per_monomer - 1
            
            global_index_map[(comp.name, copy_id)] = (start_idx, end_idx)
            
            current_atom_idx = end_idx + 1
    
    return global_index_map


def _generate_contactmap_for_domain(
    ref_structure: str,
    domain_start: int,
    domain_end: int,
    cutoff: float = 4.5
) -> List[Tuple[int, int, float]]:
    """
    为单个 domain 生成 contact pairs（局部索引）
    
    基于 restraint_generator.py 的逻辑实现。
    
    Args:
        ref_structure: 参考结构文件（PDB 或 GRO）
                      应使用 pdb2gmx 生成的 structure.gro，包含正确的氢原子和原子顺序
        domain_start: 残基起始（1-based）
        domain_end: 残基终止（1-based）
        cutoff: 距离阈值（Angstrom）
        
    Returns:
        [(atom1_local, atom2_local, distance), ...]
        原子索引是相对于 ref_structure 的局部索引（1-based for PLUMED）
    """
    # Load reference structure (MDAnalysis auto-detects PDB/GRO)
    ref_universe = mda.Universe(ref_structure)
    
    # Select heavy atoms in domain
    heavy_atoms = ref_universe.select_atoms(
        f'not name H* and resid {domain_start}-{domain_end}'
    )
    
    if len(heavy_atoms) == 0:
        return []
    
    # Calculate pairwise distances with cutoff
    pairs, distances = capped_distance(
        heavy_atoms.positions,
        heavy_atoms.positions,
        max_cutoff=cutoff,
        return_distances=True
    )
    
    # Filter: different residues and sequence distance > 1
    atom1 = heavy_atoms[pairs[:, 0]]
    atom2 = heavy_atoms[pairs[:, 1]]
    mask = (atom1.resids != atom2.resids) & (np.abs(atom1.resids - atom2.resids) > 1)
    
    selected_pairs = pairs[mask]
    selected_distances = distances[mask]
    
    # Get global indices (0-based in MDAnalysis)
    global_indices1 = heavy_atoms[selected_pairs[:, 0]].indices
    global_indices2 = heavy_atoms[selected_pairs[:, 1]].indices
    
    # Remove duplicates
    unique_pairs = set()
    keep_indices = []
    for i, (idx1, idx2) in enumerate(zip(global_indices1, global_indices2)):
        pair = tuple(sorted((idx1, idx2)))
        if pair not in unique_pairs:
            unique_pairs.add(pair)
            keep_indices.append(i)
    
    # Build result list with 1-based indices (for PLUMED)
    result = []
    for i in keep_indices:
        atom1_idx = int(global_indices1[i]) + 1  # Convert to 1-based
        atom2_idx = int(global_indices2[i]) + 1  # Convert to 1-based
        distance = float(selected_distances[i])
        result.append((atom1_idx, atom2_idx, distance))
    
    return result


def _parse_fdomains(fdomains: str, config_dir: str = None) -> List[Tuple[int, int]]:
    """
    解析 fdomains 配置（支持文件路径和内联 YAML）
    
    复用 CGSimulator._parse_fdomains 的逻辑
    
    Args:
        fdomains: fdomains 配置（文件路径或内联 YAML）
        config_dir: 配置文件目录（用于解析相对路径）
        
    Returns:
        [(start, end), ...] 列表，1-based
    """
    def _is_inline_yaml(text: str) -> bool:
        """检查是否是内联 YAML（而不是文件路径）"""
        if not text:
            return False
        stripped = text.strip()
        # 以 {、[、字母开头，或者包含换行符且有 YAML 特征
        if stripped.startswith('{') or stripped.startswith('['):
            return True
        if '\n' in stripped and (':' in stripped or stripped.startswith('-')):
            return True
        return False
    
    # 检查是否是内联YAML
    if _is_inline_yaml(fdomains):
        # 直接解析字符串
        data = yaml.safe_load(fdomains)
    else:
        # 解析文件
        fdomains_abs = fdomains
        if not os.path.isabs(fdomains):
            if config_dir:
                fdomains_abs = os.path.join(config_dir, fdomains)
        
        if not os.path.exists(fdomains_abs):
            raise FileNotFoundError(f"Domains file not found: {fdomains_abs}")
        
        with open(fdomains_abs, 'r') as f:
            data = yaml.safe_load(f)
    
    # 解析域定义
    domains = []
    if isinstance(data, dict):
        for protein_name, domain_list in data.items():
            if isinstance(domain_list, list):
                for domain in domain_list:
                    if isinstance(domain, (list, tuple)) and len(domain) == 2:
                        domains.append((domain[0], domain[1]))
    
    return domains


def _write_plumed_dat(
    contactmaps: List[Dict],
    output_file: str,
    kappa: float,
    verbose: bool
):
    """
    写入 plumed.dat 文件（严格按照 contactmap_restrain.py 格式）
    
    Args:
        contactmaps: Contact map 列表
        output_file: 输出文件路径
        kappa: Restraint 力常数
        verbose: 是否输出详细信息
    """
    with open(output_file, 'w') as f:
        # Write CONTACTMAP definitions
        q_labels = []
        for q_label, cm in enumerate(contactmaps):
            pairs = cm['pairs']
            
            if len(pairs) == 0:
                continue
            
            weight = 1.0 / len(pairs)
            
            # 开头：CONTACTMAP  ...（注意两个空格）
            f.write("CONTACTMAP  ...\n")
            
            # Write contact pairs - 每个独占一行
            for i, (atom1, atom2, distance) in enumerate(pairs, 1):
                ref_value = distance / 10.0  # Angstrom to nm
                # 注意：REF 值后面有空格和右括号: REF={value} }
                line = (f"ATOMS{i}={atom1},{atom2} "
                       f"SWITCH{i}={{Q R_0=0.01 BETA=20 LAMBDA=1.5 REF={ref_value} }} "
                       f"WEIGHT{i}={weight}\n")
                f.write(line)
            
            # 结尾：LABEL、SUM、... CONTACTMAP、PRINT
            f.write(f"LABEL=Q{q_label}\n")
            f.write("SUM\n")
            f.write("... CONTACTMAP\n")
            f.write(f"PRINT ARG=Q{q_label} FILE=COLVAR{q_label}\n")
            
            q_labels.append(q_label)
        
        # Write restraints（使用大写 Q）
        for q_label in q_labels:
            f.write(f"RESTRAINT ARG=Q{q_label} AT=1.0 KAPPA={int(kappa)} SLOPE=0.\n")
    
    if verbose:
        print(f"  Written plumed.dat to: {output_file}")
