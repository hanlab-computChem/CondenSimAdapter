"""
Coordinate building utilities for OpenMpipi.

Provides functions for:
- Generating spiral coordinates for IDPs
- Parsing PDB files for MDP structures
"""

import numpy as np
import openmm.app as app
import openmm as mm
import openmm.unit as unit
import mdtraj as md
import os
import random


def generate_spiral_coords(N, spacing=0.381):
    """
    Generates a spiral configuration of N coordinates in 3D space, often used 
    for initializing coarse-grained biomolecular models.

    Args:
        N (int): The number of points to generate.
        spacing (float, optional): The spacing between adjacent points. Default is 0.381 nm.

    Returns:
        ndarray: An N x 3 array of 3D coordinates representing the spiral, centered around the origin.
    """
    theta = np.sqrt(np.arange(N) / float(N)) * 2 * np.pi
    r = spacing * np.sqrt(np.arange(N))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = np.linspace(-N * spacing / 2, N * spacing / 2, N)
    
    points = np.column_stack((x, y, z))
    
    return points - np.mean(points, axis=0)


def parse_pdb(pdb_file, globular_indices):
    """
    Parses an atomistic PDB file to extract coarse-grained (CG) coordinates based on globular regions and returns
    the sequence of the protein. Residues are mapped to Ca.

    Args:
        pdb_file (str): Path to the PDB file to parse.
        globular_indices (list): List of lists containing indices of globular regions.

    Returns:
        tuple:
            - CG_coords (ndarray): N x 3 array of CG coordinates, centered around the origin.
            - pdb_sequence (str): The sequence of the protein extracted from the PDB file.
    """
    
    pdb = app.PDBFile(pdb_file)
    topology = pdb.topology
    aa_coords = pdb.positions

    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    residues = [residue for residue in topology.residues() if residue.name in aa_map.keys()]
    pdb_sequence = ''.join([aa_map[residue.name] for residue in residues])

    all_globular_indices = [index for domain in globular_indices for index in domain]
    
    CG_coords = []
    for index, residue in enumerate(residues):
        for atom in residue.atoms():
            if atom.name == 'CA':  
                CG_coords.append(aa_coords[atom.index]/unit.nanometer)
                break

    CG_coords = np.array(CG_coords)
    
    return CG_coords - np.mean(CG_coords, axis=0), pdb_sequence

