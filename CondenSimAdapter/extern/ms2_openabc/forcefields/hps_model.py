import numpy as np
import pandas as pd
try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit
from CondenSimAdapter.extern.ms2_openabc.forcefields.cg_model import CGModel
from CondenSimAdapter.extern.ms2_openabc.forcefields import functional_terms
from CondenSimAdapter.extern.ms2_openabc.lib import _amino_acids, _kcal_to_kj
import sys
import os

__location__ = os.path.dirname(os.path.abspath(__file__))

class HPSModel(CGModel):
    """
    The class for HPS model that represents a mixture of HPS model proteins. 
    
    This class inherits CGModel class. 
    """
    def __init__(self):
        """
        Initialize. 
        
        """
        self.atoms = None
        self.bonded_attr_names = ['protein_bonds', 'exclusions']
        self.molecules = []  # Track molecules with their parsers
        
    def add_protein_bonds(self, force_group=1):
        """
        Add protein bonds. 
        
        Parameters
        ----------
        force_group : int
            Force group. 
        
        """
        print('Add protein bonds.')
        force = functional_terms.harmonic_bond_term(self.protein_bonds, self.use_pbc, force_group)
        self.system.addForce(force)
    
    def add_contacts(self, hydropathy_scale='Urry', epsilon=0.2*_kcal_to_kj, mu=1, delta=0.08, force_group=2):
        """
        Add nonbonded contacts. 
        
        The raw hydropathy scale is scaled and shifted by: mu*lambda - delta
        
        Parameters
        ----------
        hydropathy_scale : str
            Hydropathy scale, can be KR or Urry. 
        
        epsilon : float or int
            Contact strength. 
        
        mu : float or int
            Hydropathy scale factor. 
        
        delta : float or int
            Hydropathy shift factor. 
        
        force_group : int
            Force group. 
            
        """
        print('Add nonbonded contacts.')
        resname_list = self.atoms['resname'].tolist()
        atom_types = [_amino_acids.index(x) for x in resname_list]
        if hydropathy_scale == 'KR':
            print('Use KR hydropathy scale.')
            df_contact_parameters = pd.read_csv(f'{__location__}/parameters/HPS_KR_parameters.csv')
        elif hydropathy_scale == 'Urry':
            print('Use Urry hydropathy scale.')
            df_contact_parameters = pd.read_csv(f'{__location__}/parameters/HPS_Urry_parameters.csv')
        else:
            sys.exit(f'Error: hydropathy scale {hydropathy_scale} cannot be recognized!')
        sigma_ah_map, lambda_ah_map = np.zeros((20, 20)), np.zeros((20, 20))
        for i, row in df_contact_parameters.iterrows():
            atom_type1 = _amino_acids.index(row['atom_type1'])
            atom_type2 = _amino_acids.index(row['atom_type2'])
            sigma_ah_map[atom_type1, atom_type2] = row['sigma']
            sigma_ah_map[atom_type2, atom_type1] = row['sigma']
            lambda_ah_map[atom_type1, atom_type2] = row['lambda']
            lambda_ah_map[atom_type2, atom_type1] = row['lambda']
        print(f'Scale factor mu = {mu} and shift delta = {delta}.')
        lambda_ah_map = mu*lambda_ah_map - delta
        force = functional_terms.ashbaugh_hatch_term(atom_types, self.exclusions, self.use_pbc, epsilon, 
                                                    sigma_ah_map, lambda_ah_map, force_group)
        self.system.addForce(force)
    
    def add_dh_elec(self, ldby=1*unit.nanometer, dielectric_water=80.0, cutoff=3.5*unit.nanometer, force_group=3):
        """
        Add Debye-Huckel electrostatic interactions. 
        
        Parameters
        ----------
        ldby : Quantity
            Debye length. 
        
        dielectric_water : float or int
            Dielectric constant of water. 
        
        cutoff : Quantity
            Cutoff distance. 
        
        force_group : int
            Force group. 
        
        """
        print('Add Debye-Huckel electrostatic interactions.')
        print(f'Set Debye length as {ldby.value_in_unit(unit.nanometer)} nm.')
        print(f'Set water dielectric as {dielectric_water}.')
        charges = self.atoms['charge'].tolist()
        force = functional_terms.dh_elec_term(charges, self.exclusions, self.use_pbc, ldby, dielectric_water, 
                                              cutoff, force_group)
        self.system.addForce(force)

    def add_elastic_network(self, force_constant=700.0*unit.kilojoule_per_mole/unit.nanometer**2, force_group=4):
        """
        Add elastic network restraints for folded domains (MDP support).
        
        This method collects elastic network pairs from all molecules that have 
        domain definitions and adds them to the OpenMM system.
        
        Parameters
        ----------
        force_constant : Quantity
            Force constant for the harmonic potential. Default 700 kJ/mol/nm^2 (matching CALVADOS).
        
        force_group : int
            Force group index. Default 4.
        
        """
        # Collect all elastic network pairs from molecules
        all_pairs = []
        for mol in self.molecules:
            if hasattr(mol, 'enm_pairs') and mol.enm_pairs:
                # Adjust indices based on molecule offset
                # The enm_pairs are stored with local indices relative to each molecule
                # We need to add the offset when appending to the system
                all_pairs.extend(mol.enm_pairs)
        
        if not all_pairs:
            print('No elastic network pairs found. Skipping elastic network.')
            return
        
        print(f'Adding elastic network with {len(all_pairs)} pairs.')
        force = functional_terms.elastic_network_term(all_pairs, self.use_pbc, force_constant, force_group)
        self.system.addForce(force)

    def add_elastic_network_go(self, force_constant=1.0, force_group=4):
        """
        Add Go-like elastic network restraints for folded domains (MDP support).
        
        This method collects elastic network pairs from all molecules that have 
        domain definitions and adds them to the OpenMM system using Go-like potential.
        
        Parameters
        ----------
        force_constant : float
            Scaling factor for the Go potential.
        
        force_group : int
            Force group index. Default 4.
        
        """
        # Collect all elastic network pairs from molecules
        all_pairs = []
        for mol in self.molecules:
            if hasattr(mol, 'enm_pairs') and mol.enm_pairs:
                all_pairs.extend(mol.enm_pairs)
        
        if not all_pairs:
            print('No elastic network pairs found. Skipping Go-like elastic network.')
            return
        
        print(f'Adding Go-like elastic network with {len(all_pairs)} pairs.')
        force = functional_terms.elastic_network_go_term(all_pairs, self.use_pbc, force_constant, force_group)
        self.system.addForce(force)

    def append_mol(self, new_mol, verbose=False):
        """
        Append a new molecule and track its parser for elastic network support.
        
        Parameters
        ----------
        new_mol : HPSParser
            The molecule parser object.
        
        verbose : bool
            Whether to print verbose output.
        
        """
        # Call parent append_mol
        super().append_mol(new_mol, verbose)
        
        # Track the molecule for elastic network
        self.molecules.append(new_mol)

    def add_all_default_forces(self):
        """
        Add all the forces with default settings. 
        """
        print('Add all the forces with default settings.')
        self.add_protein_bonds(force_group=1)
        self.add_contacts(force_group=2)
        self.add_dh_elec(force_group=3)
        

