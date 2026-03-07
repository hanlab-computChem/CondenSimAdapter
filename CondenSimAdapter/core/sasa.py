"""
SASA (Solvent Accessible Surface Area) calculation using Shrake-Rupley algorithm.

Internal implementation to replace external mdsim dependency.
"""

import numpy as np
from typing import List, Optional, Dict
from pathlib import Path


# VdW radii from Amber force field (Å)
_VDW_RADII: Dict[str, float] = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.50, 'S': 1.80,
    'P': 1.85, 'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98,
    'FE': 0.64, 'CA': 0.99, 'MG': 0.66, 'ZN': 0.74, 'NA': 1.02,
    'K': 1.38, 'CS': 1.67, 'CU': 0.65, 'MN': 0.65,
}

# Element guess from atom name
def _guess_element(atom_name: str) -> str:
    """Guess element from atom name.
    
    PDB atom naming conventions:
    - CA, CB, CD, etc. are carbon atoms (alpha carbon, beta carbon, etc.)
    - CL, BR, FE, MG are elements (chlorine, bromine, iron, magnesium)
    - Numbers at start indicate alternative conformations (e.g., 1HB, 2HB)
    """
    atom_name = atom_name.strip().upper()
    
    # Remove leading numbers (alternative conformation indicators)
    while atom_name and atom_name[0].isdigit():
        atom_name = atom_name[1:]
    
    if not atom_name:
        return 'C'  # default
    
    # Two-letter elements (must check before single letter)
    two_letter_elements = ('CL', 'BR', 'FE', 'MG', 'ZN', 'NA', 'CS', 'CU', 'MN')
    if len(atom_name) >= 2 and atom_name[:2] in two_letter_elements:
        return atom_name[:2]
    
    # Single letter elements
    if atom_name[0] in _VDW_RADII:
        return atom_name[0]
    
    return 'C'  # default to carbon


def _generate_sphere_points(n_points: int = 1920) -> np.ndarray:
    """
    Generate uniformly distributed points on a unit sphere using golden spiral.
    
    Args:
        n_points: Number of points to generate
        
    Returns:
        Array of shape (n_points, 3) with unit vectors
    """
    points = np.zeros((n_points, 3))
    phi = np.pi * (3.0 - np.sqrt(5.0))  # golden angle
    
    for i in range(n_points):
        y = 1 - (i / float(n_points - 1)) * 2
        radius = np.sqrt(1 - y * y)
        theta = phi * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points[i] = [x, y, z]
    
    return points


def calc_sasa_shrake_rupley(
    coordinates: np.ndarray,
    atom_radii: np.ndarray,
    probe_radius: float = 1.4,
    n_sphere_points: int = 1920,
) -> np.ndarray:
    """
    Calculate SASA for each atom using Shrake-Rupley algorithm.
    
    Args:
        coordinates: (N, 3) array of atom coordinates in Å
        atom_radii: (N,) array of VdW radii in Å
        probe_radius: Water probe radius in Å (default 1.4)
        n_sphere_points: Number of points per sphere (default 1920, matching COCOMO)
        
    Returns:
        (N,) array of SASA values in Å²
    """
    n_atoms = len(coordinates)
    if n_atoms == 0:
        return np.array([])
    
    # Generate sphere points once
    sphere_points = _generate_sphere_points(n_sphere_points)
    
    # Expanded radii (VdW + probe)
    expanded_radii = atom_radii + probe_radius
    
    # Calculate SASA for each atom
    sasa = np.zeros(n_atoms)
    
    for i in range(n_atoms):
        r_i = expanded_radii[i]
        center_i = coordinates[i]
        
        # Generate test points around atom i
        test_points = center_i + sphere_points * r_i
        
        # Check which points are exposed
        exposed = np.ones(n_sphere_points, dtype=bool)
        
        for j in range(n_atoms):
            if i == j:
                continue
            
            # Check if any test point is inside atom j's expanded sphere
            r_j = expanded_radii[j]
            center_j = coordinates[j]
            
            # Quick distance check first
            dist_sq = np.sum((center_i - center_j) ** 2)
            if dist_sq > (r_i + r_j) ** 2:
                continue  # Too far to overlap
            
            # Check each point
            diff = test_points - center_j
            dists_sq = np.sum(diff ** 2, axis=1)
            exposed[dists_sq < r_j * r_j] = False
        
        # Calculate SASA
        n_exposed = np.sum(exposed)
        sasa[i] = 4.0 * np.pi * r_i * r_i * n_exposed / n_sphere_points
    
    return sasa


def calc_sasa_from_pdb(
    pdb_path: str,
    n_sphere_points: int = 1920,
    probe_radius: float = 1.4,
) -> Optional[np.ndarray]:
    """
    Calculate per-residue SASA from a PDB file.
    
    Args:
        pdb_path: Path to PDB file (must be full-atom, not CA-only)
        n_sphere_points: Number of points per sphere (default 1920)
        probe_radius: Water probe radius in Å (default 1.4)
        
    Returns:
        Array of SASA values in nm² per residue, or None if calculation fails
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        return None
    
    # Parse PDB
    atoms = []
    residues = []  # List of (chain_id, res_id, res_name)
    current_residue = None
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip()
            res_id = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            # Track residue boundaries
            res_key = (chain_id, res_id, res_name)
            if res_key != current_residue:
                current_residue = res_key
                residues.append(res_key)
            
            element = _guess_element(atom_name)
            vdw = _VDW_RADII.get(element, 1.7)
            
            atoms.append({
                'coords': np.array([x, y, z]),
                'vdw': vdw,
                'res_idx': len(residues) - 1,
            })
    
    if not atoms:
        return None
    
    # Check if this is CA-only structure
    atom_names = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_names.append(line[12:16].strip())
    
    ca_only = all(name == 'CA' for name in atom_names)
    if ca_only:
        return None
    
    # Prepare arrays
    coords = np.array([a['coords'] for a in atoms])
    radii = np.array([a['vdw'] for a in atoms])
    res_indices = np.array([a['res_idx'] for a in atoms])
    
    # Calculate SASA
    atom_sasa = calc_sasa_shrake_rupley(coords, radii, probe_radius, n_sphere_points)
    
    # Sum per residue
    n_res = len(residues)
    residue_sasa = np.zeros(n_res)
    for i, s in enumerate(atom_sasa):
        residue_sasa[res_indices[i]] += s
    
    # Convert Å² to nm²
    return residue_sasa * 0.01  # 1 Å² = 0.01 nm²
