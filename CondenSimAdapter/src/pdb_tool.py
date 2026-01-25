#!/usr/bin/env python3
"""
PDB Utilities for Multiscale Simulations

Provides:
- Chain ID generation (supports 100+ chains)
- Coordinate extraction from PDB files

Chain ID Scheme:
- 0-61: Single character (A-Z, a-z, 0-9)
- 62+: Double character combinations (AA, AB, ..., 00, 01, ...)
"""

import os
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Sequence
from dataclasses import dataclass


# =============================================================================
# Amino Acid Conversion Utilities
# =============================================================================

# Standard amino acid 1-letter to 3-letter code conversion
AA_1_TO_3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
    # Selenocysteine and pyrrolysine (rare)
    'U': 'SEC', 'O': 'PYL',
}

def _one_to_three_letter(one_letter: str) -> str:
    """Convert 1-letter amino acid code to 3-letter code.
    
    Args:
        one_letter: Single letter amino acid code
        
    Returns:
        Three letter amino acid code (default 'UNK' for unknown)
    """
    return AA_1_TO_3.get(one_letter.upper(), 'UNK')


# =============================================================================
# Chain Label Generator
# =============================================================================

class ChainLabel:
    """
    Generate unique chain IDs supporting 100+ chains.
    
    Chain ID Scheme:
    - 0-61: Single character (A-Z, a-z, 0-9)
    - 62+: Double character combinations (AA, AB, ..., zz, 00, 01, ...)
    
    Examples:
        0 -> 'A'
        25 -> 'Z'
        26 -> 'a'
        51 -> 'z'
        52 -> '0'
        61 -> '9'
        62 -> 'AA'
        100 -> 'AQ'
    """
    
    # 62 unique single characters
    SINGLE_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    
    # 62*62 = 3844 possible double combinations
    DOUBLE_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    
    @classmethod
    def get(cls, index: int) -> str:
        """
        Get chain ID for given index.
        
        Args:
            index: Chain index (0-based)
            
        Returns:
            Chain ID string
        """
        if index < 0:
            raise ValueError(f"Chain index must be non-negative, got {index}")
        
        if index < 62:
            return cls.SINGLE_CHARS[index]
        
        # Double character: first 62*62 chains use combinations
        idx = index - 62
        first = idx // 62
        second = idx % 62
        return cls.DOUBLE_CHARS[first] + cls.DOUBLE_CHARS[second]
    
    @classmethod
    def get_range(cls, start: int, count: int) -> List[str]:
        """Get list of chain IDs for a range"""
        return [cls.get(i) for i in range(start, start + count)]
    
    @classmethod
    def validate(cls, chain_id: str) -> bool:
        """Check if a chain ID is valid"""
        if len(chain_id) == 1:
            return chain_id in cls.SINGLE_CHARS
        elif len(chain_id) == 2:
            return (chain_id[0] in cls.DOUBLE_CHARS and 
                    chain_id[1] in cls.DOUBLE_CHARS)
        return False


# =============================================================================
# Coordinate Data Structure
# =============================================================================

@dataclass
class Residue:
    """A single residue with coordinates"""
    chain_id: str
    resid: int  # 1-based residue number
    resname: str  # e.g., 'ALA', 'GLY'
    atom_name: str  # e.g., 'CA'
    coords: Tuple[float, float, float]  # (x, y, z)


# =============================================================================
# Coordinate Extraction
# =============================================================================

def extract_coordinates_from_pdb(
    pdb_path: str,
    chains: Optional[List[str]] = None,
) -> Dict[str, List[Tuple[float, float, float]]]:
    """
    Extract coordinates from a PDB file.

    Args:
        pdb_path: Path to PDB file
        chains: Optional list of chain IDs to extract (None = all)

    Returns:
        Dict mapping chain_id -> list of (x, y, z) coordinates
    """
    coordinates = {}

    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            # Parse ATOM record
            chain_id = line[21:22].strip()

            # Filter by chains if specified
            if chains is not None and chain_id not in chains:
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if chain_id not in coordinates:
                coordinates[chain_id] = []
            coordinates[chain_id].append((x, y, z))

    return coordinates


# =============================================================================
# Type hints for delayed imports
# =============================================================================

try:
    from typing import TYPE_CHECKING
except ImportError:
    TYPE_CHECKING = False

if TYPE_CHECKING:
    from .cg import CGSimulationConfig, CGComponent


# =============================================================================
# Chain Mapping Utilities
# =============================================================================

def build_chain_mapping(
    components: List,
) -> Dict[str, Tuple[str, int]]:
    """
    Build a mapping from chain index to (component_name, molecule_index).
    
    Returns:
        Dict mapping chain_index -> (component_name, mol_index)
    """
    mapping = {}
    chain_index = 0
    
    for comp in components:
        for mol_idx in range(comp.nmol):
            chain_id = ChainLabel.get(chain_index)
            mapping[chain_id] = (comp.name, mol_idx)
            chain_index += 1
    
    return mapping


def get_chain_info(
    chain_id: str,
    components: List,
) -> Optional[Tuple[str, int, int]]:
    """
    Get information about a chain.
    
    Args:
        chain_id: Chain ID (e.g., 'A', 'AA')
        components: List of CGComponent from config
        
    Returns:
        Tuple of (component_name, nres_per_mol, nmol) or None if not found
    """
    # Find chain index
    chain_index = None
    for idx in range(5000):  # Reasonable upper bound
        if ChainLabel.get(idx) == chain_id:
            chain_index = idx
            break
    
    if chain_index is None:
        return None
    
    # Find which component this chain belongs to
    cumulative = 0
    for comp in components:
        total_in_comp = comp.nmol
        if chain_index < cumulative + total_in_comp:
            mol_idx = chain_index - cumulative
            # Get nres based on type
            nres = comp.nres if hasattr(comp, 'nres') else 0
            return (comp.name, nres, mol_idx)
        cumulative += total_in_comp
    
    return None


# =============================================================================
# Validation Utilities
# =============================================================================

def validate_chain_count(components: List) -> int:
    """Calculate total number of chains needed"""
    return sum(comp.nmol for comp in components)


def validate_pdb_structure(
    pdb_path: str,
    expected_chains: int,
    expected_atoms_per_chain: int,
) -> Tuple[bool, str]:
    """
    Validate that a PDB file matches expected structure.
    
    Args:
        pdb_path: Path to PDB file
        expected_chains: Expected number of chains
        expected_atoms_per_chain: Expected atoms per chain
        
    Returns:
        Tuple of (is_valid, message)
    """
    chains = {}
    atoms_per_chain = {}
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            
            chain_id = line[21:22].strip()
            
            if chain_id not in chains:
                chains[chain_id] = 0
            chains[chain_id] += 1
    
    # Check chain count
    actual_chains = len(chains)
    if actual_chains != expected_chains:
        return False, f"Chain count mismatch: expected {expected_chains}, got {actual_chains}"
    
    # Check atoms per chain
    for chain_id, count in chains.items():
        if count != expected_atoms_per_chain:
            return False, f"Chain {chain_id} has {count} atoms, expected {expected_atoms_per_chain}"
    
    return True, f"Valid: {actual_chains} chains x {expected_atoms_per_chain} atoms"


# =============================================================================
# Module Exports
# =============================================================================

__all__ = [
    # Chain labeling
    'ChainLabel',
    
    # Data structures
    'Residue',
    
    # Coordinate handling
    'extract_coordinates_from_pdb',
    
    # Chain mapping
    'build_chain_mapping',
    'get_chain_info',
    
    # Validation
    'validate_chain_count',
    'validate_pdb_structure',
]


# Type hints for delayed imports (for cg.py integration)
if False:
    from .cg import CGSimulationConfig, CGComponent
