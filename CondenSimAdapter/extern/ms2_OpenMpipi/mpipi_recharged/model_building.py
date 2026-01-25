"""
Model building utilities for OpenMpipi.

Provides functions for:
- Calculating target box vectors
- Building combined models from multiple chains
- Equilibrating slab configurations
"""

import numpy as np
import openmm.app as app
import openmm as mm
import openmm.unit as unit
import scipy.constants as c
import os
import random
import warnings

from .system_building import get_mpipi_system
from .constants import PLATFORM, PROPERTIES

# Import grid building from CALVADOS (the canonical implementation)
try:
    from CondenSimAdapter.extern.ms2_calvados.calvados.build import build_xyzgrid as calvados_build_xyzgrid
except ImportError:
    raise ImportError("CALVADOS is required for Mpipi-Recharged. Please install CondenSimAdapter.extern.ms2_calvados.")


# ====== CALVADOS-style placement functions ======

def draw_starting_vec(box):
    """
    Draw random position within the simulation box.
    
    Args:
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.
    
    Returns:
        numpy.ndarray: Random 3D position within the box.
    """
    return np.random.random(3) * np.array(box)


def check_walls(x, box):
    """
    Check if molecule position is within the simulation box.
    
    Args:
        x (numpy.ndarray): Positions of all atoms in the molecule.
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.
    
    Returns:
        bool: True if molecule is outside the box (clashes with wall), False otherwise.
    """
    box = np.array(box)
    if np.min(x) < 0:
        return True  # clash with left box wall
    d = box - x
    if np.min(d) < 0:
        return True  # clash with right box wall
    return False  # molecule in box


def check_clash(x, pos_others, box, cutoff=0.75):
    """
    Check for clashes with other particles using distance array.
    
    Args:
        x (numpy.ndarray): Positions of atoms in the new molecule.
        pos_others (list): List of position arrays for other molecules already placed.
                           Each element can have different number of atoms.
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.
        cutoff (float): Minimum distance (in nm) to avoid clash. Default 0.75 nm.
    
    Returns:
        bool: True if clash detected, False otherwise.
    """
    boxfull = np.append(box, [90, 90, 90])  # Add angles for distance calculation
    
    if len(pos_others) == 0:
        return False  # no other particles
    
    # Check each placed molecule separately (they may have different sizes)
    overall_min = float('inf')
    for placed_pos in pos_others:
        # placed_pos is a numpy array of shape (M, 3) where M can vary
        placed_pos_arr = np.array(placed_pos)
        if len(placed_pos_arr) == 0:
            continue
        
        # Calculate minimum distances between x and this placed molecule
        dx = x[:, np.newaxis, :] - placed_pos_arr[np.newaxis, :, :]  # (N, M, 3)
        box_arr = np.array(box[:3])
        dx = dx - box_arr[np.newaxis, np.newaxis, :] * np.round(dx / box_arr[np.newaxis, np.newaxis, :])
        dists = np.sqrt(np.sum(dx**2, axis=2))  # (N, M)
        min_dist_to_this = np.min(dists)
        
        if min_dist_to_this < overall_min:
            overall_min = min_dist_to_this
    
    if overall_min < cutoff:
        return True  # clash with other particles
    else:
        return False  # no clash


def _distance_array(x, y, box):
    """
    Calculate minimum distance from each point in x to any point in y.
    
    Args:
        x (numpy.ndarray): Coordinates of first set (N, 3).
        y (numpy.ndarray): Coordinates of second set (M, 3) or (K, M, 3) for multiple molecules.
        box (numpy.ndarray): Box dimensions [Lx, Ly, Lz, alpha, beta, gamma].
    
    Returns:
        numpy.ndarray: Minimum distance array (N,) if y is 2D, or (K, N) if y is 3D.
    """
    box = np.array(box)
    Lx, Ly, Lz = box[:3]
    box_arr = np.array([Lx, Ly, Lz])
    
    # Handle y shape: if 3D (multiple molecules), flatten to 2D for distance calc
    if y.ndim == 3:
        # y shape: (K, M, 3) - K molecules, M atoms each
        K, M, _ = y.shape
        y_flat = y.reshape(-1, 3)  # (K*M, 3)
        
        # Calculate distances for each molecule separately
        result = np.zeros((K, x.shape[0]))
        for k in range(K):
            y_k = y[k]  # (M, 3)
            dx = x[:, np.newaxis, :] - y_k[np.newaxis, :, :]  # (N, M, 3)
            dx = dx - box_arr[np.newaxis, np.newaxis, :] * np.round(dx / box_arr[np.newaxis, np.newaxis, :])
            dist = np.sqrt(np.sum(dx**2, axis=2))  # (N, M)
            result[k] = np.min(dist, axis=1)  # (N,) - min distance for each atom in x
        return result  # (K, N)
    else:
        # y shape: (M, 3) - single molecule
        dx = x[:, np.newaxis, :] - y[np.newaxis, :, :]  # (N, M, 3)
        dx = dx - box_arr[np.newaxis, np.newaxis, :] * np.round(dx / box_arr[np.newaxis, np.newaxis, :])
        dist = np.sqrt(np.sum(dx**2, axis=2))  # (N, M)
        return np.min(dist, axis=1)  # (N,)


def random_placement(box, pos_others, xinit, ntries=10000, cutoff=0.75, verbose=False):
    """
    Randomly place a molecule within the box, avoiding clashes and walls.
    
    Args:
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.
        pos_others (list): List of position arrays for other molecules already placed.
        xinit (numpy.ndarray): Initial coordinates of the molecule to place (centered at origin).
        ntries (int): Maximum number of placement attempts. Default 10000.
        cutoff (float): Minimum distance to avoid clash. Default 0.75 nm.
        verbose (bool): Print verbose output.
    
    Returns:
        numpy.ndarray: Final positions of the placed molecule.
    
    Raises:
        ValueError: If unable to place molecule within ntries attempts.
    """
    for ntry in range(ntries):
        x0 = draw_starting_vec(box)  # random point in box
        xs = x0 + xinit  # shift molecule to that position
        
        # Check if within box walls
        if check_walls(xs, box):
            if verbose and ntry % 1000 == 0:
                print(f'  Try {ntry}: outside box walls', flush=True)
            continue
        
        # Check for clashes with other molecules
        if check_clash(xs, pos_others, box, cutoff=cutoff):
            if verbose and ntry % 1000 == 0:
                print(f'  Try {ntry}: clash detected', flush=True)
            continue
        
        # Success!
        if verbose and ntry > 0:
            print(f'  Placed after {ntry + 1} attempts', flush=True)
        return xs
    
    raise ValueError(f'Failed to place molecule after {ntries} attempts. '
                     f'Consider increasing box size or reducing molecule count.')


def build_xyzgrid(N, box):
    """
    Create a 3D grid for placing molecules.

    Args:
        N (int): Number of molecules to place.
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.

    Returns:
        numpy.ndarray: Grid positions (N, 3) in nm.
    """
    box = np.array(box)
    r = box / np.sum(box)
    a = np.cbrt(N / np.product(r))
    n = a * r
    nxyz = np.floor(n)
    
    while np.product(nxyz) < N:
        ndeviation = n / nxyz
        devmax = np.argmax(ndeviation)
        nxyz[devmax] += 1
    while np.product(nxyz) > N:
        nmax = np.argmax(nxyz)
        nxyz[nmax] -= 1
        if np.product(nxyz) < N:
            nxyz[nmax] += 1
            break

    xyz = []
    x, y, z = 0., 0., 0.
    ctx, cty, ctz = 0, 0, 0
    dx = box[0] / nxyz[0]
    dy = box[1] / nxyz[1]
    dz = box[2] / nxyz[2]
    zplane = 1
    xyplane = 1

    for _ in np.arange(N):
        if zplane > 0:
            xshift, yshift = 0, 0
        else:
            xshift, yshift = dx/2, dy/2

        if xyplane < 0:
            zshift = dz/2
        else:
            zshift = 0

        xyz.append([x + xshift, y + yshift, z + zshift])

        ctx += 1
        x += dx

        if ctx == nxyz[0]:
            ctx = 0
            x = 0
            cty += 1
            y += dy
            if cty == nxyz[1]:
                ctx = 0
                cty = 0
                x = 0.
                y = 0.
                ctz += 1
                z += dz
                zplane = -zplane

        if (ctx % 2 == cty % 2):
            xyplane = 1
        else:
            xyplane = -1

    xyz = np.asarray(xyz)
    return xyz

# Use CALVADOS's build_xyzgrid (canonical implementation)
build_xyzgrid = calvados_build_xyzgrid

def calculate_target_box_vectors(chain_info, 
                                 box_size=None,
                                 target_density=0.1*unit.gram/unit.centimeter**3):
    """
    Calculate the initial target box vectors for the system.
    
    If box_size is provided, use it directly.
    Otherwise, calculate from total mass and target density.

    Args:
        chain_info (dict): 
            A dictionary with keys = chain objects and values = number of copies to place.
            Each chain object must have:
              - chain_mass (in daltons)
              - min_rg_coords (initial coordinates)
        box_size (list, optional): 
            Box size [Lx, Ly, Lz] in nm. If provided, use this directly.
        target_density (Quantity, optional): 
            Target density as an OpenMM `Quantity` (default 0.1 g/cm^3).

    Returns:
        Quantity: 
            A 3x3 array (in nm) of box vectors suitable for periodic boundary conditions.
    """
    if box_size is not None:
        # Use provided box size directly
        Lx, Ly, Lz = box_size
        return np.array([[Lx, 0, 0], 
                         [0, Ly, 0], 
                         [0, 0, Lz]]) * unit.nanometer
    
    # Calculate from total mass and target density
    total_mass_g = sum([chain.chain_mass.value_in_unit(unit.dalton)/c.Avogadro*unit.gram * n_copies 
                        for chain, n_copies in chain_info.items()], 0*unit.gram)
    
    target_volume = total_mass_g / target_density
    short_side_length = (target_volume / 6) ** (1/3)  # default long_side_scale_factor=6
    short_side_length = short_side_length.in_units_of(unit.nanometer) / unit.nanometer

    return short_side_length * np.array([[6, 0, 0], 
                                         [0, 1, 0], 
                                         [0, 0, 1]]) * unit.nanometer


# ====== Model building ======

def _copy_topology_with_new_chain_id(topology, new_chain_id):
    """
    Create a copy of a topology with a new chain ID.
    
    Args:
        topology (Topology): Source topology to copy.
        new_chain_id (str): New chain ID to use.
    
    Returns:
        Topology: A new Topology object with the updated chain ID.
    """
    new_topo = app.Topology()
    
    for old_chain in topology.chains():
        new_chain = new_topo.addChain(id=new_chain_id)
        
        for old_residue in old_chain.residues():
            new_residue = new_topo.addResidue(old_residue.name, new_chain, old_residue.id, old_residue.insertionCode)
            
            for old_atom in old_residue.atoms():
                new_topo.addAtom(old_atom.name, old_atom.element, new_residue, old_atom.id, old_atom.formalCharge)
    
    return new_topo


def _get_unique_chain_id(index, total_chains):
    """
    Generate a unique chain ID for the given index.
    
    Generates IDs like 'A', 'B', 'C', ..., 'Z', 'A1', 'B1', ... for up to 676 chains.
    
    Args:
        index (int): The 0-based index of the chain.
        total_chains (int): Total number of chains (for logging only).
    
    Returns:
        str: A unique chain ID.
    """
    # Standard PDB chain IDs: A-Z (26), then A1-Z1, A2-Z2, etc.
    if index < 26:
        return chr(ord('A') + index)
    else:
        quotient = index // 26
        remainder = index % 26
        return chr(ord('A') + remainder) + str(quotient)


def _build_model_with_gmx_insert_molecules(chain_info, target_box_vectors, topol='cubic', radius=0.35, verbose=False):
    """
    Build a Modeller object using GROMACS insert-molecules tool.
    
    This function uses gmx insert-molecules to place multiple copies of molecules
    in the simulation box, which is more robust than manual placement.
    
    Args:
        chain_info (dict): 
            A dictionary with keys = chain objects and values = number of copies to place.
            Each chain object must have:
              - topology (OpenMM Topology)
              - min_rg_coords (initial coords)
        target_box_vectors (Quantity): 
            A 3x3 array of box vectors (in nm).
        topol (str): 
            Topology type: 'slab', 'grid', or 'cubic' (default: 'cubic').
        radius (float): 
            Minimum distance (nm) between atoms of different molecules. Default 0.35 nm.
        verbose (bool): 
            Print verbose output during placement.

    Returns:
        Modeller: 
            A Modeller object containing the combined system of all chains, 
            set with appropriately extended periodic box vectors.
    """
    import tempfile
    import shutil
    
    try:
        import gromacs
    except ImportError:
        raise ImportError("GromacsWrapper is required for gmx insert-molecules. Please install it using: pip install GromacsWrapper")
    
    print('Initializing the model using gmx insert-molecules...', flush=True)
    
    # Get box dimensions
    box_vectors = target_box_vectors.value_in_unit(unit.nanometer)
    Lx, Ly, Lz = box_vectors[0, 0], box_vectors[1, 1], box_vectors[2, 2]
    box = np.array([Lx, Ly, Lz])
    
    # Adjust box for slab topology
    # Use full z dimension for SLAB mode to ensure all molecules can be inserted
    if topol == 'slab':
        slab_width = Lz  # Use full z instead of z/2
        insert_box = np.array([Lx, Ly, slab_width])
        print(f'  SLAB mode: using z-range [0, {slab_width:.2f}] nm (full z)', flush=True)
    else:
        insert_box = box
    
    # Create temporary directory for GROMACS files
    temp_dir = tempfile.mkdtemp(prefix='gmx_insert_')
    
    try:
        # Process each unique chain type
        all_chains = []
        for chain_idx, (chain_obj, n_copies) in enumerate(chain_info.items()):
            # Write monomer PDB file
            chain_id = getattr(chain_obj, 'chain_id', f'chain_{chain_idx}')
            monomer_pdb = os.path.join(temp_dir, f'monomer_{chain_idx}.pdb')
            monomer_topology = _copy_topology_with_new_chain_id(chain_obj.topology, 'A')
            monomer_positions = chain_obj.min_rg_coords * unit.nanometer
            
            # Write PDB file
            with open(monomer_pdb, 'w') as f:
                app.PDBFile.writeFile(monomer_topology, monomer_positions, f)
            
            # Verify PDB file was written correctly
            pdb_check = app.PDBFile(monomer_pdb)
            n_atoms_in_pdb = len(list(pdb_check.topology.atoms()))
            n_atoms_expected = len(chain_obj.min_rg_coords)
            if n_atoms_in_pdb != n_atoms_expected:
                raise ValueError(f"PDB file has {n_atoms_in_pdb} atoms, expected {n_atoms_expected}")
            
            if verbose:
                print(f'  Prepared monomer {chain_id}: {n_copies} copies ({n_atoms_in_pdb} atoms)', flush=True)
            
            # Use gmx insert-molecules to place molecules (output PDB format)
            output_pdb = os.path.join(temp_dir, f'output_{chain_idx}.pdb')
            
            # For first chain type, create empty box
            if len(all_chains) == 0:
                # Use subprocess directly (gromacs wrapper has issues with box parameter)
                import subprocess
                cmd = [
                    'gmx', 'insert-molecules',
                    '-ci', monomer_pdb,
                    '-nmol', str(n_copies),
                    '-box', f'{insert_box[0]:.5f}', f'{insert_box[1]:.5f}', f'{insert_box[2]:.5f}',
                    '-radius', f'{radius:.5f}',
                    '-o', output_pdb
                ]
                if verbose:
                    print(f'    Running: {" ".join(cmd)}', flush=True)
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)
                if result.returncode != 0:
                    raise RuntimeError(f"gmx insert-molecules failed:\n{result.stderr}\n{result.stdout}")
                # Verify all molecules were inserted
                if os.path.exists(output_pdb):
                    pdb_check = app.PDBFile(output_pdb)
                    n_atoms_inserted = len(list(pdb_check.topology.atoms()))
                    n_atoms_per_mol = len(chain_obj.min_rg_coords)
                    n_mols_inserted = n_atoms_inserted // n_atoms_per_mol
                    if verbose:
                        print(f'    Inserted {n_mols_inserted} molecules ({n_atoms_inserted} atoms, expected {n_copies * n_atoms_per_mol})', flush=True)
                    if n_mols_inserted != n_copies:
                        raise RuntimeError(f"gmx insert-molecules only inserted {n_mols_inserted} molecules instead of {n_copies}. Box may be too small (current: {insert_box}) or radius too large ({radius} nm).")
            else:
                # For subsequent chain types, insert into existing structure
                previous_pdb = all_chains[-1][2]  # Use the output_pdb from previous chain
                # Use subprocess directly
                import subprocess
                cmd = [
                    'gmx', 'insert-molecules',
                    '-ci', monomer_pdb,
                    '-nmol', str(n_copies),
                    '-f', previous_pdb,
                    '-radius', f'{radius:.5f}',
                    '-o', output_pdb
                ]
                if verbose:
                    print(f'    Running: {" ".join(cmd)}', flush=True)
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)
                if result.returncode != 0:
                    raise RuntimeError(f"gmx insert-molecules failed:\n{result.stderr}\n{result.stdout}")
                # Verify all molecules were inserted
                if os.path.exists(output_pdb):
                    pdb_check = app.PDBFile(output_pdb)
                    n_atoms_inserted = len(list(pdb_check.topology.atoms()))
                    n_atoms_per_mol = len(chain_obj.min_rg_coords)
                    n_mols_inserted = (n_atoms_inserted - len(list(app.PDBFile(previous_pdb).topology.atoms()))) // n_atoms_per_mol
                    if verbose:
                        print(f'    Inserted {n_mols_inserted} molecules ({n_atoms_inserted} total atoms)', flush=True)
                    if n_mols_inserted != n_copies:
                        raise RuntimeError(f"gmx insert-molecules only inserted {n_mols_inserted} molecules instead of {n_copies}. Box may be too small or radius too large.")
            
            all_chains.append((chain_obj, n_copies, output_pdb))
        
        # Read final PDB file using OpenMM PDBFile
        final_pdb = all_chains[-1][2]
        pdb_file = app.PDBFile(final_pdb)
        final_topology = pdb_file.topology
        final_positions = pdb_file.positions
        
        # Apply slab offset if needed
        if topol == 'slab':
            offset_z = Lz / 2.0 - slab_width / 2.0
            positions_array = np.array([[p.x, p.y, p.z] for p in final_positions])
            positions_array[:, 2] += offset_z
            final_positions = [unit.Quantity(np.array([p[0], p[1], p[2]]), unit.nanometer) for p in positions_array]
        
        # Rebuild topology with unique chain IDs
        # Use original chain topologies and assign unique chain IDs to each molecule copy
        new_topology = app.Topology()
        chain_idx = 0
        total_chains = sum(n for _, n, _ in all_chains)
        
        # Count atoms per chain type to split positions
        atoms_per_molecule = {}
        for chain_obj, n_copies, _ in all_chains:
            n_atoms_per_mol = len(chain_obj.min_rg_coords)
            atoms_per_molecule[chain_obj] = n_atoms_per_mol
        
        # Process each chain type
        atom_idx = 0
        for chain_obj, n_copies, _ in all_chains:
            n_atoms_per_mol = atoms_per_molecule[chain_obj]
            
            for copy_idx in range(n_copies):
                chain_id = _get_unique_chain_id(chain_idx, total_chains)
                new_chain = new_topology.addChain(id=chain_id)
                
                # Copy topology from original chain object
                for old_residue in chain_obj.topology.residues():
                    new_residue = new_topology.addResidue(old_residue.name, new_chain, old_residue.id)
                    for old_atom in old_residue.atoms():
                        new_topology.addAtom(old_atom.name, old_atom.element, new_residue, old_atom.id)
                
                atom_idx += n_atoms_per_mol
                chain_idx += 1
        
        # Extract positions for each molecule from PDB file
        # Convert positions to numpy array format
        all_positions_array = np.array([[p.x, p.y, p.z] for p in final_positions])
        
        # Verify total atoms match
        total_atoms_expected = sum(n * atoms_per_molecule[chain_obj] for chain_obj, n, _ in all_chains)
        if len(all_positions_array) != total_atoms_expected:
            raise ValueError(f"Position array size mismatch: expected {total_atoms_expected} atoms, got {len(all_positions_array)}")
        
        new_positions = []
        atom_idx = 0
        
        for chain_obj, n_copies, _ in all_chains:
            n_atoms_per_mol = atoms_per_molecule[chain_obj]
            for copy_idx in range(n_copies):
                if atom_idx + n_atoms_per_mol > len(all_positions_array):
                    raise ValueError(f"Not enough positions: need {atom_idx + n_atoms_per_mol}, have {len(all_positions_array)}")
                mol_positions_array = all_positions_array[atom_idx:atom_idx + n_atoms_per_mol]
                # Convert to OpenMM Quantity format
                for pos in mol_positions_array:
                    new_positions.append(unit.Quantity(pos, unit.nanometer))
                atom_idx += n_atoms_per_mol
        
        # Verify final positions count matches topology
        total_atoms_in_topology = sum(1 for _ in new_topology.atoms())
        if len(new_positions) != total_atoms_in_topology:
            raise ValueError(f"Topology/position mismatch: topology has {total_atoms_in_topology} atoms, positions have {len(new_positions)}")
        
        # Create Modeller with rebuilt topology and positions
        model = app.Modeller(new_topology, new_positions)
        
        # Set box vectors
        model.topology.setPeriodicBoxVectors(target_box_vectors)
        
        print(f'  ✓ Successfully placed {sum(n for _, n, _ in all_chains)} molecules using gmx insert-molecules', flush=True)
        
        return model
        
    finally:
        # Cleanup temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)


def _parse_gro_file(gro_file):
    """Parse GROMACS GRO file and return positions and topology."""
    positions = []
    atoms = []
    residues = []
    chains = []
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    # Skip title line
    n_atoms = int(lines[1].strip())
    
    # Parse atoms
    for i in range(2, 2 + n_atoms):
        line = lines[i]
        resnum = int(line[0:5].strip())
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnum = int(line[15:20].strip())
        
        x = float(line[20:28]) * 0.1  # Convert from Angstrom to nm
        y = float(line[28:36]) * 0.1
        z = float(line[36:44]) * 0.1
        
        positions.append(unit.Quantity(np.array([x, y, z]), unit.nanometer))
        atoms.append(atomname)
        residues.append((resnum, resname))
    
    # Create OpenMM topology
    topology = app.Topology()
    chain = topology.addChain('A')
    
    current_resnum = None
    current_residue = None
    for i, (resnum, resname) in enumerate(residues):
        if resnum != current_resnum:
            current_resnum = resnum
            current_residue = topology.addResidue(resname, chain)
        element = app.Element.getBySymbol(atoms[i][0]) if atoms[i][0].isalpha() else app.Element.getBySymbol('C')
        topology.addAtom(atoms[i], element, current_residue)
    
    return positions, topology


def build_model(chain_info, target_box_vectors, topol='cubic', use_grid_placement=True, clash_cutoff=0.75, verbose=False, use_gmx_insert=False, gmx_radius=0.35):
    """
    Build a Modeller object with the given chain information and box vectors.
    
    This function places chains in the box using either:
    1. GROMACS insert-molecules (if use_gmx_insert=True) - recommended for robustness
    2. CALVADOS-style grid placement (if use_grid_placement=True)
    3. Simple center placement (if use_grid_placement=False)
    
    Args:
        chain_info (dict): 
            A dictionary with keys = chain objects and values = number of copies to place.
            Each chain object must have:
              - topology (OpenMM Topology)
              - min_rg_coords (initial coords)
        target_box_vectors (Quantity): 
            A 3x3 array of box vectors (in nm).
        topol (str): 
            Topology type: 'slab', 'grid', or 'cubic' (default: 'cubic').
        use_grid_placement (bool): 
            If True and use_gmx_insert=False, use CALVADOS-style grid placement.
            If False and use_gmx_insert=False, center each chain at its original position.
        clash_cutoff (float): 
            Minimum distance (nm) between atoms of different molecules. Default 0.75 nm.
            Only used if use_gmx_insert=False.
        verbose (bool): 
            Print verbose output during placement.
        use_gmx_insert (bool):
            If True, use GROMACS insert-molecules tool instead of manual placement.
            Default False for backward compatibility.
        gmx_radius (float):
            Minimum distance (nm) for gmx insert-molecules. Default 0.35 nm.
            Only used if use_gmx_insert=True.

    Returns:
        Modeller: 
            A Modeller object containing the combined system of all chains, 
            set with appropriately extended periodic box vectors.
    """
    # Use GROMACS insert-molecules if requested
    if use_gmx_insert:
        return _build_model_with_gmx_insert_molecules(
            chain_info, target_box_vectors, topol=topol, radius=gmx_radius, verbose=verbose
        )
    print('Initializing the model...', flush=True)
    
    # Get box dimensions
    box_vectors = target_box_vectors.value_in_unit(unit.nanometer)
    Lx, Ly, Lz = box_vectors[0, 0], box_vectors[1, 1], box_vectors[2, 2]
    box = np.array([Lx, Ly, Lz])
    
    # Build chain list (expand copies)
    chain_list = []
    for chain, n_copies in chain_info.items():
        chain_list.extend([(chain, copy_idx) for copy_idx in range(n_copies)])
    
    total_chains = len(chain_list)
    print(f'  Total chains to place: {total_chains}', flush=True)
    print(f'  Topology: {topol}', flush=True)
    
    if total_chains == 0:
        raise ValueError("No chains to place in chain_info")
    
    if use_grid_placement:
        # ====== CALVADOS-style grid placement ======
        print(f'  Using CALVADOS-style grid placement', flush=True)
        
        # Calculate grid dimensions based on topology (CALVADOS approach)
        if topol == 'slab':
            # SLAB mode: use reduced z-range (slab_width = box[2] / 2)
            slab_width = Lz / 2.0
            grid_box = np.array([Lx, Ly, slab_width])
            print(f'  SLAB mode: slab_width = {slab_width:.2f} nm', flush=True)
        else:
            # GRID or CUBIC mode: use full box
            grid_box = box
        
        # Pre-calculate grid positions
        print(f'  Calculating grid positions for {total_chains} molecules...', flush=True)
        grid_positions = build_xyzgrid(total_chains, grid_box)
        
        # Apply offset for SLAB mode (CALVADOS logic: center the slab in z)
        if topol == 'slab':
            offset_z = Lz / 2.0 - slab_width / 2.0
            grid_positions[:, 2] += offset_z
            print(f'  SLAB offset: z += {offset_z:.2f} nm', flush=True)
        
        print(f'  Grid range X: [{grid_positions[:,0].min():.2f}, {grid_positions[:,0].max():.2f}]', flush=True)
        print(f'  Grid range Y: [{grid_positions[:,1].min():.2f}, {grid_positions[:,1].max():.2f}]', flush=True)
        print(f'  Grid range Z: [{grid_positions[:,2].min():.2f}, {grid_positions[:,2].max():.2f}]', flush=True)
        
        # Place chains using grid positions with minimum separation distance
        # Minimum distance between any two atoms of different molecules: 4 Å = 0.4 nm
        min_separation = 0.4  # nm
        
        # First, calculate molecule size from first chain (all chains should be similar)
        first_chain, first_copy_idx = chain_list[0]
        xinit = first_chain.min_rg_coords
        xinit_centered = xinit - xinit.mean(axis=0)
        mol_size = np.abs(xinit_centered).max(axis=0)  # Max extent in each dimension
        
        # Adjust grid positions to ensure molecules fit in box
        # Shift grid positions inward if they're too close to boundaries
        grid_positions_adjusted = grid_positions.copy()
        for i in range(len(grid_positions_adjusted)):
            # Ensure grid position + molecule size fits in box
            grid_positions_adjusted[i] = np.clip(
                grid_positions_adjusted[i],
                mol_size,
                box - mol_size
            )
        
        # Place first chain
        x0 = grid_positions_adjusted[0]
        first_positions = x0 + xinit_centered
        
        # Final check - if still outside, center it (shouldn't happen after adjustment)
        if np.min(first_positions) < 0 or np.any(box - first_positions < 0):
            print(f'  WARNING: First chain still outside box, using safe center position...', flush=True)
            first_positions = box / 2.0 + xinit_centered
        
        # Generate unique chain ID for first chain
        first_chain_id = _get_unique_chain_id(0, total_chains)
        first_topology = _copy_topology_with_new_chain_id(first_chain.topology, first_chain_id)
        
        model = app.Modeller(first_topology, first_positions * unit.nanometer)
        
        # Store placed molecule positions for clash checking
        placed_positions = [first_positions]
        
        if verbose:
            print(f'  Placing chain 1/{total_chains} ({first_chain_id})...', flush=True)
        
        # Place remaining chains with separation check
        for chain_idx in range(1, total_chains):
            chain_obj, copy_idx = chain_list[chain_idx]
            xinit = chain_obj.min_rg_coords
            xinit_centered = xinit - xinit.mean(axis=0)
            
            # Generate unique chain ID for this copy
            chain_id = _get_unique_chain_id(chain_idx, total_chains)
            
            if verbose:
                print(f'  Placing chain {chain_idx + 1}/{total_chains} ({chain_id})...', flush=True)
            
            # Start with adjusted grid position
            x0 = grid_positions_adjusted[chain_idx].copy()
            final_positions = x0 + xinit_centered
            
            # Check for clashes with already placed molecules
            clash_detected = check_clash(final_positions, placed_positions, box, cutoff=min_separation)
            
            # If clash detected, try to adjust position iteratively
            max_adjustment_iterations = 10
            for adjust_iter in range(max_adjustment_iterations):
                if not clash_detected:
                    break
                    
                if verbose and adjust_iter == 0:
                    print(f'    Clash detected, adjusting position to ensure {min_separation:.2f} nm separation...', flush=True)
                
                # Find the closest placed molecule and calculate actual minimum distance
                min_dist = float('inf')
                closest_pos = None
                closest_placed_pos = None
                for placed_pos in placed_positions:
                    # Calculate all pairwise distances
                    dx = final_positions[:, np.newaxis, :] - placed_pos[np.newaxis, :, :]
                    dx = dx - box * np.round(dx / box)
                    dists = np.sqrt(np.sum(dx**2, axis=2))
                    min_dist_to_this = np.min(dists)
                    if min_dist_to_this < min_dist:
                        min_dist = min_dist_to_this
                        closest_pos = placed_pos.mean(axis=0)
                        closest_placed_pos = placed_pos
                
                if closest_pos is not None and min_dist < min_separation:
                    # Calculate direction from closest molecule to new molecule
                    new_center = final_positions.mean(axis=0)
                    direction = new_center - closest_pos
                    # Apply minimum image convention
                    direction = direction - box * np.round(direction / box)
                    
                    # If molecules are too close or direction is zero, use random direction
                    if np.linalg.norm(direction) < 1e-6:
                        # Use a direction based on chain index to avoid all moving same way
                        angles = [0, 2*np.pi/3, 4*np.pi/3]
                        angle = angles[chain_idx % 3]
                        direction = np.array([np.cos(angle), np.sin(angle), 0.5])
                        direction = direction / np.linalg.norm(direction)
                    else:
                        direction = direction / np.linalg.norm(direction)
                    
                    # Move away by required separation distance plus molecule radius
                    # Use actual minimum distance to calculate required movement
                    required_separation = min_separation + 2 * mol_size.max()
                    movement = required_separation - min_dist + 0.1  # Add small buffer
                    new_center = closest_pos + direction * movement
                    
                    # Wrap to box
                    new_center = new_center - box * np.floor(new_center / box)
                    # Ensure new center allows molecule to fit
                    new_center = np.clip(new_center, mol_size.max(), box - mol_size.max())
                    
                    # Update positions
                    final_positions = new_center + xinit_centered
                    
                    # Wrap positions to box
                    final_positions = final_positions - box * np.floor(final_positions / box)
                    
                    # Re-check clash
                    clash_detected = check_clash(final_positions, placed_positions, box, cutoff=min_separation)
                else:
                    break
            
            # If still clashing after adjustments, use a fallback position
            if clash_detected:
                if verbose:
                    print(f'    WARNING: Still clashing after adjustments, using fallback position...', flush=True)
                # Use a position far from all placed molecules
                # Calculate a position that maximizes distance to all placed centers
                placed_centers = [pos.mean(axis=0) for pos in placed_positions]
                # Try positions on a grid around the box center
                grid_offsets = [
                    [0, 0, 0],
                    [min_separation + mol_size.max(), 0, 0],
                    [0, min_separation + mol_size.max(), 0],
                    [0, 0, min_separation + mol_size.max()],
                ]
                best_pos = None
                best_min_dist = 0
                for offset in grid_offsets:
                    candidate_center = box / 2.0 + np.array(offset)
                    candidate_center = candidate_center - box * np.floor(candidate_center / box)
                    candidate_pos = candidate_center + xinit_centered
                    candidate_pos = candidate_pos - box * np.floor(candidate_pos / box)
                    
                    # Check minimum distance to all placed molecules
                    min_dist_to_all = float('inf')
                    for placed_pos in placed_positions:
                        dx = candidate_pos[:, np.newaxis, :] - placed_pos[np.newaxis, :, :]
                        dx = dx - box * np.round(dx / box)
                        dists = np.sqrt(np.sum(dx**2, axis=2))
                        min_dist_to_all = min(min_dist_to_all, np.min(dists))
                    
                    if min_dist_to_all > best_min_dist:
                        best_min_dist = min_dist_to_all
                        best_pos = candidate_pos
                
                if best_pos is not None:
                    final_positions = best_pos
                else:
                    # Last resort: use offset from center
                    safe_offset = (chain_idx % 3) * (min_separation + mol_size.max())
                    final_positions = box / 2.0 + xinit_centered + np.array([safe_offset, safe_offset, safe_offset])
                    final_positions = final_positions - box * np.floor(final_positions / box)
            
            # Final check - ensure in box
            if np.min(final_positions) < 0 or np.any(box - final_positions < 0):
                print(f'  WARNING: Chain {chain_idx + 1} still outside box, using safe position...', flush=True)
                safe_offset = (chain_idx % 3) * (min_separation + mol_size.max())
                final_positions = box / 2.0 + xinit_centered + np.array([safe_offset, safe_offset, safe_offset])
                final_positions = final_positions - box * np.floor(final_positions / box)
            
            # Verify separation (final check)
            final_clash = check_clash(final_positions, placed_positions, box, cutoff=min_separation)
            if final_clash and verbose:
                # Calculate actual minimum distance
                for placed_pos in placed_positions:
                    dx = final_positions[:, np.newaxis, :] - placed_pos[np.newaxis, :, :]
                    dx = dx - box * np.round(dx / box)
                    dists = np.sqrt(np.sum(dx**2, axis=2))
                    min_dist = np.min(dists)
                    print(f'    WARNING: Minimum distance to placed molecule: {min_dist:.3f} nm (target: {min_separation:.3f} nm)', flush=True)
            
            # Create topology with unique chain ID
            chain_topology = _copy_topology_with_new_chain_id(chain_obj.topology, chain_id)
            model.add(chain_topology, final_positions * unit.nanometer)
            
            # Add to placed positions for next iteration
            placed_positions.append(final_positions)
    
    else:
        # ====== Simple center placement (original behavior) ======
        print('  Using center placement (all chains at original positions)', flush=True)
        
        first_chain, first_copy_idx = chain_list[0]
        
        # Generate unique chain ID for first chain
        first_chain_id = _get_unique_chain_id(0, total_chains)
        first_topology = _copy_topology_with_new_chain_id(first_chain.topology, first_chain_id)
        
        model = app.Modeller(first_topology, first_chain.min_rg_coords * unit.nanometer)
        
        for chain_idx in range(1, total_chains):
            chain_obj, copy_idx = chain_list[chain_idx]
            
            # Generate unique chain ID for this copy
            chain_id = _get_unique_chain_id(chain_idx, total_chains)
            chain_topology = _copy_topology_with_new_chain_id(chain_obj.topology, chain_id)
            
            model.add(chain_topology, chain_obj.min_rg_coords * unit.nanometer)

    # Set periodic box vectors
    model.topology.setPeriodicBoxVectors(target_box_vectors)
    
    print('The model is built.', flush=True)
    print('The topology of the model:', model.topology)
    
    return model


# ====== System building ======

def build_mpipi_recharged_system_from_chains(chain_info,
                                              box_size=None,
                                              topol='cubic',
                                              T=280*unit.kelvin,
                                              csx=150,
                                              CM_remover=True,
                                              periodic=True,
                                              use_grid_placement=True,
                                              clash_cutoff=0.75,
                                              verbose=False,
                                              use_gmx_insert=False,
                                              gmx_radius=0.35):
    """
    Build an OpenMM System using Mpipi-Recharged forcefield from chain objects.
    
    This function:
    1. Calls get_compact_model() on each chain to relax and calculate Rg
    2. Builds a combined model using CALVADOS-style grid placement
    3. Adds forces using get_mpipi_system()
    
    Note: To use external coordinates (e.g., from CALVADOS pre-equilibration),
    set them after building the system by calling:
        simulation.context.setPositions(calvados_positions)
    
    Args:
        chain_info (dict):
            A dictionary with keys = chain objects (MDP, IDP, RNA, etc.) and values = number of copies.
            Each chain object must have:
              - topology (OpenMM Topology)
              - initial_coords (initial coordinates)
        box_size (list, optional):
            Box size [Lx, Ly, Lz] in nm. If provided, use this directly.
        topol (str, optional):
            Topology type: 'slab', 'grid', or 'cubic' (default: 'cubic').
            SLAB mode uses reduced z-range for initial placement.
        T (Quantity, optional):
            Temperature in Kelvin (default 280 K).
        csx (float, optional):
            Ionic strength in mM (default 150).
        CM_remover (bool, optional):
            Whether to add center-of-mass motion remover (default True).
        periodic (bool, optional):
            Whether to use periodic boundary conditions (default True).
        use_grid_placement (bool, optional):
            Whether to use CALVADOS-style grid placement (default True, recommended for large systems).
            Ignored if use_gmx_insert=True.
        clash_cutoff (float, optional):
            Minimum distance (nm) between atoms of different molecules (default 0.75 nm).
            Only used if use_gmx_insert=False.
        verbose (bool, optional):
            Print verbose output during placement.
        use_gmx_insert (bool, optional):
            If True, use GROMACS insert-molecules tool instead of manual placement (default False).
            This is more robust for large systems and avoids overlap issues.
        gmx_radius (float, optional):
            Minimum distance (nm) for gmx insert-molecules (default 0.35 nm).
            Only used if use_gmx_insert=True.
    
    Returns:
        tuple: (system, model) where:
            - system: An OpenMM System object ready for simulation
            - model: A Modeller object containing the combined topology and positions
    """
    
    # Calculate target box vectors first
    print('Calculating target box vectors...', flush=True)
    target_box_vectors = calculate_target_box_vectors(chain_info, box_size=box_size)
    print(f'  Box size: {target_box_vectors[0,0]} x {target_box_vectors[1,1]} x {target_box_vectors[2,2]}')
    print(f'  Topology: {topol}', flush=True)
    
    # Always use grid placement (standard workflow)
    # Note: If you want to use external coordinates, do it after this function
    # by calling simulation.context.setPositions() with external coordinates
    
    # ===== Standard workflow: grid placement =====
    print('Relaxing monomers...', flush=True)
    
    # Step 1: Relax each chain and calculate Rg
    for chain in chain_info.keys():
        print(f'  Relaxing {chain.chain_id}...', flush=True)
        chain.get_compact_model(simulation_time=10*unit.nanosecond)
    
    # Step 2: Build combined model
    print('Building combined model...', flush=True)
    
    total_chains = sum(chain_info.values())
    if use_gmx_insert:
        print(f'  Using gmx insert-molecules for {total_chains} chains (radius={gmx_radius} nm)...', flush=True)
        model = build_model(
            chain_info, target_box_vectors, 
            topol=topol,
            use_gmx_insert=True,
            gmx_radius=gmx_radius,
            verbose=verbose
        )
    elif use_grid_placement and total_chains > 1:
        print(f'  Using CALVADOS-style grid placement for {total_chains} chains...', flush=True)
        model = build_model(
            chain_info, target_box_vectors, 
            topol=topol,
            use_grid_placement=True, 
            clash_cutoff=clash_cutoff,
            verbose=verbose
        )
    else:
        print(f'  Placing {total_chains} chain(s) at center...', flush=True)
        model = build_model(
            chain_info, target_box_vectors, 
            topol=topol,
            use_grid_placement=False
        )
    
    # Step 4: Build system with forces
    print('Building Mpipi-Recharged system...', flush=True)
    
    # Build globular_indices_dict: map each topology chain ID to its globular indices
    # This is critical for multi-copy systems where the topology has many chains
    # (e.g., 'A', 'B', 'C', 'A1', 'B1', ...) but chain_info only has unique biomolecule objects
    globular_indices_dict = {}
    
    # Get all chains from the combined topology
    topology_chains = list(model.topology.chains())
    
    # Build chain_list: [(chain_obj, copy_idx), ...] for all copies
    chain_list = []
    for chain_obj, n_copies in chain_info.items():
        for copy_idx in range(n_copies):
            chain_list.append((chain_obj, copy_idx))
    
    if verbose:
        print(f"  [Debug] {len(topology_chains)} topology chains, {len(chain_list)} chain_list entries")
    
    # Safety check: topology chains should match chain_list length
    if len(topology_chains) != len(chain_list):
        raise ValueError(
            f"Mismatch between topology chains ({len(topology_chains)}) and chain_list ({len(chain_list)}). "
            f"This indicates a bug in build_model()."
        )
    
    # Match each topology chain to its corresponding biomolecule and globular indices
    for i, topo_chain in enumerate(topology_chains):
        chain_obj, copy_idx = chain_list[i]
        topo_chain_id = topo_chain.id
        
        # Get the biomolecule's globular_indices (local to that chain)
        globular_indices_dict[topo_chain_id] = chain_obj.globular_indices
        
        if verbose:
            print(f"  [Debug] Chain '{topo_chain_id}' (copy {copy_idx}) -> globular_indices = {chain_obj.globular_indices}")
    
    if verbose:
        print(f"  [Debug] globular_indices_dict has {len(globular_indices_dict)} entries")
        for k, v in list(globular_indices_dict.items())[:3]:
            print(f"    '{k}': {v}")
    
    system = get_mpipi_system(
        np.array(model.positions),
        model.topology,
        globular_indices_dict,
        T.value_in_unit(unit.kelvin),
        csx,
        CM_remover=CM_remover,
        periodic=periodic
    )
    
    # Override box vectors set by get_mpipi_system (which calculates from position ranges)
    # Use the target box vectors we calculated earlier instead
    if periodic:
        box_vecs = [
            mm.Vec3(x=target_box_vectors[0,0].value_in_unit(unit.nanometer), y=0.0, z=0.0),
            mm.Vec3(x=0.0, y=target_box_vectors[1,1].value_in_unit(unit.nanometer), z=0.0),
            mm.Vec3(x=0.0, y=0.0, z=target_box_vectors[2,2].value_in_unit(unit.nanometer))
        ] * unit.nanometer
        system.setDefaultPeriodicBoxVectors(*box_vecs)
    
    print('System built successfully.', flush=True)
    return system, model


# ====== Equilibration functions ======

def equilibrate_slab(model, 
                     target_box_vectors, 
                     chain_info, 
                     T=280*unit.kelvin, 
                     csx=150, 
                     pulling_time=20*unit.nanosecond, 
                     equi_time=400*unit.nanosecond):
    """
    Perform a two-stage equilibration:
    
    1) Pull all chains toward the center of the box to remove large voids (a "slab" pulling).
    2) Remove the pulling force, set the final box vectors, and equilibrate further.

    Args:
        model (Modeller): 
            A Modeller object with the combined system and initial positions.
        target_box_vectors (Quantity): 
            A 3x3 array of box vectors (in nm).
        chain_info (dict): 
            A dictionary with keys = chain objects and values = number of copies. 
            Used to retrieve `globular_indices`.
        T (Quantity, optional): 
            Temperature to run the simulation, default 280 K.
        csx (float, optional): 
            Ionic strength in mM, used to compute Debye length internally.
        pulling_time (Quantity, optional): 
            Duration of the pulling simulation (default 20 ns).
        equi_time (Quantity, optional): 
            Duration of the equilibration simulation after pulling (default 400 ns).
    """
    print('Setting up integrator...', flush=True)
    integrator = mm.LangevinMiddleIntegrator(T, 0.01/unit.picosecond, 10*unit.femtosecond)

    # Prepare the globular indices dictionary
    globular_indices_dict = {
        chain.chain_id: chain.globular_indices for chain in chain_info.keys()
    }

    print('Building mpipi system...', flush=True)
    # Calls your second-code function get_mpipi_system
    system = get_mpipi_system(
        np.array(model.positions), 
        model.topology, 
        globular_indices_dict, 
        T.value_in_unit(unit.kelvin), 
        csx, 
        CM_remover=False, 
        periodic=True
    )

    # Create a gentle pulling force to drive everything toward the box center
    pulling_force = mm.CustomExternalForce(
        'k*periodicdistance(x, y, z, x0, y0, z0)^2'
    )
    midpoint_x = 0.5*target_box_vectors[0][0]
    midpoint_y = 0.5*target_box_vectors[1][1]
    midpoint_z = 0.5*target_box_vectors[2][2]
    
    pulling_force.addGlobalParameter('x0', midpoint_x)
    pulling_force.addGlobalParameter('y0', midpoint_y)
    pulling_force.addGlobalParameter('z0', midpoint_z)
    pulling_force.addPerParticleParameter('k')

    # Add a small spring constant for each atom
    for atom in model.topology.atoms():
        pulling_force.addParticle(atom.index, [0.001])
    system.addForce(pulling_force)
    print('Pulling force is set.', flush=True)

    # Setting up the Simulation
    print('Initializing simulation...', flush=True)
    simulation = app.Simulation(
        model.topology, 
        system, 
        integrator, 
        platform=PLATFORM, 
        platformProperties=PROPERTIES
    )
    simulation.context.setPositions(model.positions)
    simulation.context.setPeriodicBoxVectors(*model.topology.getPeriodicBoxVectors())

    # Minimize
    simulation.minimizeEnergy()
    print('Energy minimized.', flush=True)

    # Pulling stage
    simulation.reporters.append(
        app.StateDataReporter(
            './output_pulling.dat', 10000, 
            step=True, potentialEnergy=True, 
            temperature=True, density=True, elapsedTime=True
        )
    )
    simulation.reporters.append(app.XTCReporter('./traj_pulling.xtc', 10000))
    print('Beginning pulling...', flush=True)
    simulation.step(int(pulling_time/(10*unit.femtosecond)))
    print('Pulling complete.', flush=True)

    # Remove pulling force
    system.removeForce(system.getNumForces()-1)
    print('Pulling force removed.', flush=True)

    # Update box vectors
    simulation.context.setPeriodicBoxVectors(*target_box_vectors)
    print('Box vectors updated.', flush=True)
    
    # Re-initialize the context so that box changes take effect
    state = simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    simulation.context.reinitialize()
    simulation.context.setState(state)
    print('Context reinitialized.', flush=True)

    # Minimize again
    model.topology.setPeriodicBoxVectors(target_box_vectors)
    simulation.minimizeEnergy()

    # Clear old reporters
    simulation.reporters = []

    # Equilibration stage
    simulation.reporters.append(
        app.StateDataReporter(
            './output_equi.dat', 10000, 
            step=True, potentialEnergy=True, 
            temperature=True, elapsedTime=True
        )
    )
    simulation.reporters.append(app.XTCReporter('./traj_equi.xtc', 50000))
    print('Beginning equilibration...', flush=True)
    simulation.step(int(equi_time/(10*unit.femtosecond)))
    
    # Final state
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    app.PDBFile.writeFile(model.topology, positions, open('./equi_model.pdb', 'w'))
    simulation.saveState('equi_state.xml')
    print('Equilibration complete.', flush=True)


def build_and_equilibrate_model(chain_info, 
                                box_size=None,
                                target_density=0.1*unit.gram/unit.centimeter**3, 
                                T=280*unit.kelvin, 
                                csx=150,
                                pulling_time=20*unit.nanosecond, 
                                equi_time=1000*unit.nanosecond,
                                use_random_placement=True,
                                clash_cutoff=0.75,
                                verbose=False):
    """
    Wrapper function to:
      1) Relax each chain monomer in isolation.
      2) Calculate and build the combined model with the desired box size (CALVADOS-style).
      3) Perform slab-pulling followed by a full equilibration.

    Args:
        chain_info (dict):
            A dictionary with keys = chain objects and values = number of copies. 
            Each chain object must have:
              - get_compact_model(simulation_time): a method to relax the chain individually
              - chain_mass, min_rg_coords, max_rg, min_rg
        box_size (list, optional):
            Box size [Lx, Ly, Lz] in nm. If provided, use this directly.
        target_density (Quantity, optional):
            Target density in g/cm^3.
        T (Quantity, optional):
            Simulation temperature.
        csx (float, optional):
            Ionic strength in mM (for Debye length calculations).
        pulling_time (Quantity, optional):
            Duration of the slab-pulling simulation.
        equi_time (Quantity, optional):
            Duration of the final equilibration simulation.
        use_random_placement (bool, optional):
            Whether to use CALVADOS-style random placement (default True).
        clash_cutoff (float, optional):
            Minimum distance (nm) between atoms of different molecules (default 0.75 nm).
        verbose (bool, optional):
            Print verbose output during placement.
    """
    
    print('Relaxing monomers...')
    for chain in chain_info.keys():
        chain.get_compact_model(simulation_time=10*unit.nanosecond)
    
    # Build combined model
    target_box_vectors = calculate_target_box_vectors(chain_info, 
                                                      box_size=box_size,
                                                      target_density=target_density)
    
    model = build_model(
        chain_info, target_box_vectors,
        use_random_placement=use_random_placement,
        clash_cutoff=clash_cutoff,
        verbose=verbose
    )

    # Equilibrate
    equilibrate_slab(model, target_box_vectors, chain_info, 
                     T=T, csx=csx, 
                     pulling_time=pulling_time, 
                     equi_time=equi_time)
