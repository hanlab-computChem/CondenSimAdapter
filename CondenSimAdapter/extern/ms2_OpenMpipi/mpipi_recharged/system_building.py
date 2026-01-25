"""
System building module for OpenMpipi.

Provides the core get_mpipi_system function for building OpenMM systems
using the Mpipi-Recharged forcefield.
"""

import numpy as np

import openmm as mm
import openmm.app as app
import openmm.unit as unit

from scipy.spatial import KDTree
import os
import copy

from .constants import DATA_DIR


def get_harmonic_bonds(positions, topology, globular_indices_dict, IDR_k=8031., verbose=True):
    """
    Generates a list of harmonic bonds for a given topology. The function handles nucleosome chains, DNA chains,
    and other biomolecules and assigns bonds between atoms based on their positions and chain type.

    Args:
        positions (ndarray): An array of atomic positions.
        topology (Topology): An OpenMM topology object containing the chains and atoms of the system.
        globular_indices_dict (dict): A dictionary mapping chain IDs to lists of globular domain indices.
        dyad_positions (list): A list of dyad positions for nucleosomes.
        constraints (str, optional): Type of nucleosome-DNA constraint ('none', 'dyad', 'inner', 'breathing', 'all'). Default is 'inner'.
        IDR_k (float, optional): The force constant for harmonic bonds in IDRs. Default is 8031.
        verbose (bool): Whether to print debugging information. Default is True.
    
    Returns:
        list: A list of harmonic bonds, where each bond is represented as a tuple (atom1, atom2, distance, force constant).
    """
    
    bonds = []

    all_atoms = list(topology.atoms())  # Fetch all atoms in topology once
    atom_tree = KDTree(positions)

    chains = list(topology.chains())  # Fetch chains only once

    if verbose:
        print(f"\n  [ENM Debug] globular_indices_dict = {globular_indices_dict}")

    # Iterate over all chains in topology
    for i, chain in enumerate(chains):
        chain_atoms = list(chain.atoms())
        chain_id = chain.id

        # All IDR regions use the same spacing (0.381 nm for proteins)
        IDR_d = 0.381
        
        globular_indices_list = globular_indices_dict.get(chain_id, [])  # Use default empty list if not found

        if verbose:
            print(f"\n  [ENM Debug] Chain '{chain_id}': {len(chain_atoms)} atoms")
            print(f"  [ENM Debug] globular_indices_list = {globular_indices_list}")

        # Flatten nested lists if needed
        all_globular_indices = []
        for domain in globular_indices_list:
            if isinstance(domain, (list, tuple)):
                # Domain is [start, end] inclusive, 0-based
                all_globular_indices.extend(list(range(domain[0], domain[1] + 1)))
            else:
                # Domain is a single index
                all_globular_indices.append(domain)

        if verbose:
            print(f"  [ENM Debug] all_globular_indices (flattened) = {all_globular_indices[:20]}..." if len(all_globular_indices) > 20 else f"  [ENM Debug] all_globular_indices = {all_globular_indices}")

        # Identify IDR indices (residues NOT in globular regions)
        IDR_indices = [i for i in range(len(chain_atoms)) if i not in all_globular_indices]

        if verbose:
            print(f"  [ENM Debug] IDR indices count: {len(IDR_indices)}")
            print(f"  [ENM Debug] Globular indices count: {len(all_globular_indices)}")

        # Add bonds for IDR regions (harmonic bonds between consecutive residues)
        idr_bond_count = 0
        for idx in range(len(chain_atoms) - 1):
            if idx in IDR_indices or idx + 1 in IDR_indices:
                bonds.append((chain_atoms[idx], chain_atoms[idx + 1], IDR_d, IDR_k))
                idr_bond_count += 1

        if verbose:
            print(f"  [ENM Debug] IDR harmonic bonds added: {idr_bond_count}")

        # Add ENM bonds for globular regions
        total_enm_bonds = 0
        for domain_idx, globular_indices in enumerate(globular_indices_list):
            if isinstance(globular_indices, (list, tuple)):
                # Domain is [start, end] inclusive
                domain_start = globular_indices[0]
                domain_end = globular_indices[1]
                domain_atoms = [chain_atoms[j] for j in range(domain_start, domain_end + 1)]
            else:
                # Single index (shouldn't happen normally)
                domain_atoms = [chain_atoms[globular_indices]]

            if verbose:
                # Print which global indices are being used
                global_indices = [a.index for a in domain_atoms]
                print(f"  [ENM Debug] Domain {domain_idx + 1}: local [{domain_start}, {domain_end}] -> global {global_indices[:10]}..." if len(global_indices) > 10 else f"  [ENM Debug] Domain {domain_idx + 1}: local [{domain_start}, {domain_end}] -> global {global_indices}")

            if len(domain_atoms) >= 2:
                ENM_bonds = get_ENM_bonds(atom_tree, domain_atoms, all_atoms, verbose=verbose)
                bonds.extend(ENM_bonds)
                total_enm_bonds += len(ENM_bonds)
                if verbose:
                    print(f"  [ENM Debug]   ENM bonds added: {len(ENM_bonds)}")
            else:
                if verbose:
                    print(f"  [ENM Debug]   Skipping domain (less than 2 atoms)")

        if verbose:
            print(f"  [ENM Debug] Total ENM bonds for chain '{chain_id}': {total_enm_bonds}")

    if verbose:
        total_bonds = len(bonds)
        print(f"\n  [ENM Debug] Total bonds in system: {total_bonds}")

    return bonds


def get_ENM_bonds(atom_tree, ENM_atoms, all_atoms, cutoff=0.75, k=8031., verbose=False):
    """
    Generates a list of Elastic Network Model (ENM) bonds for a given set of atoms based on a distance cutoff.

    Args:
        atom_tree (KDTree): KDTree of atomic positions for efficient neighbor lookup.
        ENM_atoms (list): List of atoms to consider for ENM bonds.
        all_atoms (list): List of all atoms in the system.
        cutoff (float, optional): Distance cutoff for considering a bond (in nm). Default is 0.75 nm.
        k (float, optional): Force constant for the ENM bonds. Default is 8031.
        verbose (bool): Whether to print debugging information. Default is False.
    
    Returns:
        list: A list of ENM bonds, where each bond is represented as a tuple (atom1, atom2, distance, force constant).
    """
    
    bonds = []
    num_atoms = len(ENM_atoms)
    
    # Precompute positions of ENM atoms from the KDTree
    ENM_positions = {atom.index: atom_tree.data[atom.index] for atom in ENM_atoms}
    
    # Iterate through each atom in ENM_atoms
    for i in range(num_atoms):
        atom1 = ENM_atoms[i]
        atom1_pos = ENM_positions[atom1.index]
        
        # Query nearby atoms within the cutoff distance
        nearby_indices = atom_tree.query_ball_point(atom1_pos, cutoff)
        
        # Check each nearby atom and form a bond
        for j in nearby_indices:
            if j != atom1.index and j in ENM_positions.keys():  # Exclude self
                atom2 = all_atoms[j]
                atom2_pos =  ENM_positions[atom2.index]
             
                r = np.linalg.norm(atom1_pos - atom2_pos)
                bonds.append((atom1, atom2, r, k))
    
    return bonds


def calculate_debye_length(T, csx):
    """
    Calculate the Debye length based on the temperature and ionic strength.

    This function calculates the Debye length based on the temperature and ionic strength of the system.
    The Debye length is used to screen electrostatic interactions in the system.

    Args:
        T (float): The temperature of the system in Kelvin.
        csx (float): The ionic strength of the system in mM.

    Returns:
        float: The Debye length in nanometers.
    """
    # Electrolyte solution ionic strength in m^-3
    cs = (csx / 1000) * 6.022e26
    # Relative dielectric constant of the medium
    er = 5321 / T + 233.76 - 0.9297 * T + 0.001417 * T**2 - 0.0000008292 * T**3
    # Bjerrum length in meters
    bjerrum = (1.671e-5) / (er * T)
    # Debye length in nanometers
    debye_length = np.sqrt(1 / (8 * np.pi * bjerrum * cs)) * 1e9

    return debye_length


def get_mpipi_system(positions, topology, globular_indices_dict, T, csx, CM_remover=True, periodic=True):
    """
    Build an OpenMM System using the Mpipi-Recharged forcefield.

    This function creates a complete OpenMM System with:
    - Harmonic bond forces for IDR regions
    - Elastic Network Model (ENM) for globular domains
    - Custom short-range potential (WCA-like)
    - Yukawa electrostatic screening

    This function supports two topology naming conventions:
    1. OpenMpipi format: atom name = 'p' + single-letter amino acid (e.g., 'pA', 'pN')
    2. Standard PDB format: atom name = 'CA', residue name = three-letter code (e.g., 'ALA')

    Args:
        positions (ndarray): An array of atomic positions (shape: [n_atoms, 3]).
        topology (Topology): An OpenMM topology object.
        globular_indices_dict (dict): Dictionary mapping chain IDs to lists of globular domain indices.
        T (float): Temperature in Kelvin.
        csx (float): Ionic strength in mM.
        CM_remover (bool): Whether to add center-of-mass motion remover. Default is True.
        periodic (bool): Whether to use periodic boundary conditions. Default is True.

    Returns:
        System: A complete OpenMM System ready for simulation.
    """
    
    NB_PARAMETERS = np.loadtxt(os.path.join(DATA_DIR, 'recharged_params.txt'))
    
    # Mpipi-Recharged residue mapping: name -> [mass, index]
    mpipi_mapping = {
        'pM': [131.20, 0], 'pG': [57.05, 1], 'pK': [128.20, 2], 'pT': [101.10, 3],
        'pR': [156.20, 4], 'pA': [71.08, 5], 'pD': [115.10, 6], 'pE': [129.10, 7],
        'pY': [163.20, 8], 'pV': [99.07, 9], 'pL': [113.20, 10], 'pQ': [128.10, 11],
        'pW': [186.20, 12], 'pF': [147.20, 13], 'pS': [87.08, 14], 'pH': [137.10, 15],
        'pN': [114.10, 16], 'pP': [97.12, 17], 'pC': [103.10, 18], 'pI': [113.20, 19],
        'rU': [244.20, 20]
    }
    
    # Standard PDB residue name mapping: three-letter -> single-letter
    pdb_aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    # Check if topology uses OpenMpipi format or standard PDB format
    first_atom = next(iter(topology.atoms()), None)
    if first_atom is not None:
        # Detect format based on atom name
        if first_atom.name == 'CA':
            # Standard PDB format - need to use residue name to determine type
            use_pdb_format = True
        else:
            # OpenMpipi format
            use_pdb_format = False
    else:
        use_pdb_format = False
    
    # Build residue type list for atom index lookup
    atom_residue_types = []
    for atom in topology.atoms():
        if use_pdb_format:
            # Use residue name to determine residue type
            resname = atom.residue.name
            if resname in pdb_aa_map:
                single_letter = pdb_aa_map[resname]
                residue_type = 'p' + single_letter
            else:
                # Unknown residue, default to 'pG' (glycine)
                residue_type = 'pG'
        else:
            # Use atom name directly (OpenMpipi format)
            residue_type = atom.name
        atom_residue_types.append(residue_type)
    
    debye_length = calculate_debye_length(T, csx)
    system = mm.System()
    
    # Add particles using detected residue types
    for residue_type in atom_residue_types:
        if residue_type in mpipi_mapping:
            system.addParticle(mpipi_mapping[residue_type][0])
        else:
            # Fallback: use default mass for unknown residue type
            print(f"  Warning: Unknown residue type '{residue_type}', using default mass")
            system.addParticle(71.08)  # Default to Alanine mass

        
    harm_bonds = get_harmonic_bonds(positions, topology, globular_indices_dict) 
    
    bond_flag = True
    if len(list(topology.bonds())) > 0:
        bond_flag = False
    
    harm_potential = mm.HarmonicBondForce()
    for bond in harm_bonds:
        a1, a2, d, k = bond  
        
        harm_potential.addBond(a1.index, a2.index, d, k)
        if bond_flag == True:
            topology.addBond(a1, a2)
    

    system.addForce(harm_potential)
    
    wf_string = '''
    glob_factor * step(rc-r)*epsilon*alpha*((sigma/r)^(2*mu)-1)*((rc/r)^(2*mu)-1)^2;
    alpha = 2*(3^(2*mu))*((3)/(2*((3^(2*mu))-1)))^3;
    rc = 3*sigma;

    glob_factor = select(globular1*globular2, 0.7, select(globular1+globular2, sqrt(0.7), 1));
    
    epsilon = wf_table(index1, index2, 0);
    sigma   = wf_table(index1, index2, 1);
    mu =  floor(wf_table(index1, index2, 2));
    '''
    yukawa_string = '''
    (A_table(index1, index2)/r) * exp(-kappa*r);
    '''
    
    wf_potential = mm.CustomNonbondedForce(wf_string)
    yukawa_potential = mm.CustomNonbondedForce(yukawa_string)
    wf_potential.addPerParticleParameter('index')
    wf_potential.addPerParticleParameter('globular')
    
    yukawa_potential.addPerParticleParameter('index')
    yukawa_potential.addGlobalParameter('kappa', 1/debye_length)
    
    if periodic == True:
        wf_potential.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        yukawa_potential.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        
    else:
        wf_potential.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        yukawa_potential.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    wf_potential.setCutoffDistance(2.5*unit.nanometer)                        
    
    yukawa_potential.setCutoffDistance(3.5*unit.nanometer)
    yukawa_potential.setForceGroup(1) # to use different cutoff, have to be in different ForceGroup
    
    for chain in topology.chains():
        atoms = list(chain.atoms())
        if 'nuc' in chain.id: 
            
            for i, atom in enumerate(atoms):
                globular = 1 if atom.element.symbol == 'Pt' else 0
                residue_type = atom_residue_types[atom.index]
                if residue_type in mpipi_mapping:
                    index = mpipi_mapping[residue_type][1]
                else:
                    index = 5  # Default to Alanine index
             
                wf_potential.addParticle([index, globular])
                yukawa_potential.addParticle([index])
        
        else:
            for i, atom in enumerate(atoms):
                globular = 1 if atom.element.symbol == 'Pt' else 0
                residue_type = atom_residue_types[atom.index]
                if residue_type in mpipi_mapping:
                    index = mpipi_mapping[residue_type][1]
                else:
                    index = 5  # Default to Alanine index
                
                wf_potential.addParticle([index, globular]) 
                yukawa_potential.addParticle([index])
    
                
    wf_potential.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    yukawa_potential.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    
    wf_table = mm.Discrete3DFunction(21, 21, 3, NB_PARAMETERS[:-21*21])
    wf_potential.addTabulatedFunction('wf_table', wf_table)
    
    yukawa_table = mm.Discrete2DFunction(21, 21, NB_PARAMETERS[-21*21:])
    yukawa_potential.addTabulatedFunction('A_table', yukawa_table)
    
    system.addForce(wf_potential)
    system.addForce(yukawa_potential)
    
    if CM_remover == True: 
        system.addForce(mm.CMMotionRemover(1000))
    
    ranges = np.ptp(positions, axis=0)
    max_range_value = ranges[np.argmax(ranges)]
    box_length = max_range_value + 50
    box_vecs = [mm.Vec3(x=box_length, y=0.0, z=0.0), mm.Vec3(x=0.0, y=box_length, z=0.0), mm.Vec3(x=0.0, y=0.0, z=box_length)]*unit.nanometer
    system.setDefaultPeriodicBoxVectors(*box_vecs)
    
    return system

