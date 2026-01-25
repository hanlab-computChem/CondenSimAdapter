import numpy as np
import pandas as pd
from CondenSimAdapter.extern.ms2_openabc.utils import parse_pdb, atomistic_pdb_to_ca_pdb
from CondenSimAdapter.extern.ms2_openabc.lib import _amino_acids, _kcal_to_kj
import sys
import os
import yaml

_hps_amino_acid_mass_dict = dict(ALA=71.08, ARG=156.20, ASN=114.10, ASP=115.10, CYS=103.10, 
                                 GLN=128.10, GLU=129.10, GLY=57.05, HIS=137.10, ILE=113.20, 
                                 LEU=113.20, LYS=128.20, MET=131.20, PHE=147.20, PRO=97.12, 
                                 SER=87.08, THR=101.10, TRP=186.20, TYR=163.20, VAL=99.07)

_hps_amino_acid_charge_dict = dict(ALA=0.0, ARG=1.0, ASN=0.0, ASP=-1.0, CYS=0.0, 
                                   GLN=0.0, GLU=-1.0, GLY=0.0, HIS=0.5, ILE=0.0,
                                   LEU=0.0, LYS=1.0, MET=0.0, PHE=0.0, PRO=0.0,
                                   SER=0.0, THR=0.0, TRP=0.0, TYR=0.0, VAL=0.0)


class HPSParser(object):
    """
    HPS protein parser with MDP (Multi-Domain Protein) support.
    """
    def __init__(self, ca_pdb, fdomains=None, default_parse=True):
        """
        Initialize a protein with HPS model.
        
        Parameters
        ----------
        ca_pdb : str
            Path for the CA pdb file. 
        
        fdomains : str, optional
            Path to domain definition YAML file. 
            Format:
                protein_name:
                  - [start, end]  # first folded domain (1-based, inclusive)
                  - [start, end]  # second folded domain
        
        default_parse : bool
            Whether to parse with default settings. 
        
        """
        self.pdb = ca_pdb
        self.fdomains = fdomains
        self.atoms = parse_pdb(ca_pdb)
        # check if all the atoms are protein CA atoms
        assert ((self.atoms['resname'].isin(_amino_acids)).all() and self.atoms['name'].eq('CA').all())
        
        # Initialize MDM-related attributes
        self.domains = None
        self.domain_distances = None
        self.enm_pairs = None
        
        if default_parse:
            print('Parse molecule with default settings.')
            self.parse_mol()
            if fdomains:
                self.parse_domains()
    
    @classmethod
    def from_atomistic_pdb(cls, atomistic_pdb, ca_pdb, fdomains=None, write_TER=False, default_parse=True):
        """
        Initialize an HPS model protein from atomistic pdb. 
        
        Parameters
        ----------
        atomistic_pdb : str
            Path for the atomistic pdb file. 
        
        ca_pdb : str
            Output path for the CA pdb file. 
        
        fdomains : str, optional
            Path to domain definition YAML file. 
        
        write_TER : bool
            Whether to write TER between two chains. 
        
        default_parse : bool
            Whether to parse with default settings.
        
        Returns
        ------- 
        result : class instance
            A class instance. 
        
        """
        atomistic_pdb_to_ca_pdb(atomistic_pdb, ca_pdb, write_TER)
        result = cls(ca_pdb, fdomains, default_parse)
        return result
    
    def parse_mol(self, exclude12=True, mass_dict=_hps_amino_acid_mass_dict, 
                  charge_dict=_hps_amino_acid_charge_dict):
        """
        Parse molecule. 
        
        Parameters
        ----------
        exclude12 : bool
            Whether to exclude nonbonded interactions between 1-2 atoms. 
        
        mass_dict : dict
            Mass dictionary. 
        
        charge_dict : dict
            Charge dictionary. 
        
        """
        bonds = []
        n_atoms = len(self.atoms.index)
        for atom1 in range(n_atoms):
            chain1 = self.atoms.loc[atom1, 'chainID']
            if atom1 < n_atoms - 1:
                atom2 = atom1 + 1
                chain2 = self.atoms.loc[atom2, 'chainID']
                if chain1 == chain2:
                    bonds.append([atom1, atom2])
        bonds = np.array(bonds)
        self.protein_bonds = pd.DataFrame(bonds, columns=['a1', 'a2'])
        self.protein_bonds.loc[:, 'r0'] = 0.38
        self.protein_bonds.loc[:, 'k_bond'] = 2000*_kcal_to_kj
        if exclude12:
            self.exclusions = self.protein_bonds[['a1', 'a2']].copy()
        else:
            self.exclusions = pd.DataFrame(columns=['a1', 'a2'])
        # set mass and charge
        for i, row in self.atoms.iterrows():
            self.atoms.loc[i, 'mass'] = mass_dict[row['resname']]
            self.atoms.loc[i, 'charge'] = charge_dict[row['resname']]
    
    def parse_domains(self, cutoff_distance=0.9, min_seq_sep=3):
        """
        Parse domain definitions and build reference structures for elastic network.
        
        Parameters
        ----------
        cutoff_distance : float
            Cutoff distance for elastic network (in nm). Default 0.9 nm.
        
        min_seq_sep : int
            Minimum sequence separation for elastic network pairs. Default 3.
        
        """
        if not self.fdomains:
            print('Warning: No domain file specified, skipping domain parsing.')
            return
        
        # 1. Parse domain YAML file
        print(f'Parsing domain file: {self.fdomains}')
        with open(self.fdomains, 'r') as f:
            domain_data = yaml.safe_load(f)
        
        # Get the protein name (use filename without extension as default)
        protein_name = os.path.splitext(os.path.basename(self.fdomains))[0]
        if protein_name in domain_data:
            domains = domain_data[protein_name]
        else:
            # Try to find the first key in the YAML
            protein_name = list(domain_data.keys())[0]
            domains = domain_data[protein_name]
        
        # Convert to 0-based indices
        self.domains = []
        for domain in domains:
            if isinstance(domain[0], list):
                # Subdomain format
                for subdom in domain:
                    start = subdom[0] - 1
                    end = subdom[1]
                    self.domains.append(list(range(start, end)))
            else:
                start = domain[0] - 1
                end = domain[1]
                self.domains.append(list(range(start, end)))
        
        print(f'Found {len(self.domains)} folded domains:')
        for i, domain in enumerate(self.domains):
            print(f'  Domain {i+1}: residues {domain[0]+1} to {domain[-1]+1} ({len(domain)} residues)')
        
        # 2. Extract reference coordinates from PDB (already in self.atoms)
        self.ref_positions = self._get_reference_positions()
        
        # 3. Compute distance matrix within domains
        self.domain_distances = self._compute_domain_distances()
        
        # 4. Build elastic network pairs
        self.enm_pairs = self._build_elastic_network(cutoff_distance, min_seq_sep)
        
        print(f'Elastic network: {len(self.enm_pairs)} pairs added.')
    
    def _get_reference_positions(self):
        """
        Get reference positions from PDB coordinates.
        
        Returns
        -------
        positions : np.ndarray
            Array of shape (n_atoms, 3) with coordinates in nm.
        
        """
        positions = np.zeros((len(self.atoms), 3))
        for i, row in self.atoms.iterrows():
            positions[i, 0] = row['x'] / 10.0  # Convert from Angstrom to nm
            positions[i, 1] = row['y'] / 10.0
            positions[i, 2] = row['z'] / 10.0
        return positions
    
    def _compute_domain_distances(self):
        """
        Compute pairwise distances within each domain.
        
        Returns
        -------
        distances : dict
            Dictionary mapping domain index to distance matrix.
        
        """
        distances = {}
        for i, domain in enumerate(self.domains):
            n_res = len(domain)
            dist_matrix = np.zeros((n_res, n_res))
            for j in range(n_res):
                for k in range(j + 1, n_res):
                    pos_j = self.ref_positions[domain[j]]
                    pos_k = self.ref_positions[domain[k]]
                    dist = np.linalg.norm(pos_j - pos_k)
                    dist_matrix[j, k] = dist
                    dist_matrix[k, j] = dist
            distances[i] = dist_matrix
        return distances
    
    def _build_elastic_network(self, cutoff_distance, min_seq_sep):
        """
        Build elastic network pairs based on reference structure.
        
        Parameters
        ----------
        cutoff_distance : float
            Cutoff distance for elastic network (in nm).
        
        min_seq_sep : int
            Minimum sequence separation for elastic network pairs.
        
        Returns
        -------
        pairs : list
            List of tuples (i, j, r0) for elastic network.
        
        """
        pairs = []
        for domain_idx, domain in enumerate(self.domains):
            n_res = len(domain)
            dist_matrix = self.domain_distances[domain_idx]
            
            for j in range(n_res):
                for k in range(j + min_seq_sep, n_res):
                    r0 = dist_matrix[j, k]
                    if r0 < cutoff_distance:
                        # Convert back to global indices
                        global_i = domain[j]
                        global_j = domain[k]
                        pairs.append((global_i, global_j, r0))
        
        return pairs
    
    def has_elastic_network(self):
        """
        Check if this parser has elastic network constraints.
        
        Returns
        -------
        bool : True if elastic network pairs exist.
        
        """
        return self.enm_pairs is not None and len(self.enm_pairs) > 0
