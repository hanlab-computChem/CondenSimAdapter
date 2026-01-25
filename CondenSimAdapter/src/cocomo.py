from __future__ import print_function
import numpy as np
from openmm.unit import *
from openmm import *
from openmm.app import *
from parmed import load_file, unit as u
from parmed.openmm.reporters import NetCDFReporter
import parmed as pmd
from parmed.openmm import topsystem

import argparse
import pickle
import os
import yaml
from pathlib import Path
from typing import Dict, List, Optional

# Import from same package
from .cg import CGSimulationConfig, CGSimulator


class COCOMO:
    """COCOMO force field system builder using CGSimulator topology interface."""
    
    # Force field parameters (class constants)
    KBOND = 4184
    THETA0 = 180
    L0_PRO = 0.38
    L0_RNA = 0.5
    KANGLE_PRO = 4.184
    KANGLE_RNA = 5.021
    CATIONPI_PROPRO = 0.30
    CATIONPI_PRORNA = 0.20
    PIPI_PROPRO = 0.10
    EPS_POLAR = 0.176
    EPS_NOPOL = 0.295
    AZERO_POLAR = 0
    AZERO_HYDRO = 0.0002
    KAPPA = 1
    
    # Elastic network parameters
    FORCE_CONSTANT = 500.0
    CUTOFF_DISTANCE = 0.9 * nanometer
    
    # Residue sets
    PROTEIN_RESIDUES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                        'THR', 'TRP', 'TYR', 'VAL']
    RNA_RESIDUES = ['ADE', 'CYT', 'GUA', 'URA']
    
    # Amino acid 1-letter to 3-letter mapping
    AA_1TO3 = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    
    def __init__(self, box_size: float, topology_info: Dict, positions: np.ndarray,
                 surf: float = 0.7, resources: str = 'CUDA'):
        """
        Initialize COCOMO system builder using topology info from CGSimulator.
        
        Args:
            box_size: Box length in nanometers
            topology_info: Dictionary containing topology information from CGSimulator
                - 'global_sequence': Global sequence string (e.g., "ACD...")
                - 'chain_ids': List of chain IDs for each residue (1-based)
                - 'folded_domains': List of fold domain flags (0/1) for each residue
                - 'component_names': List of component names for each residue
                - 'local_residue_indices': List of local residue indices (1-based)
            positions: Numpy array of atom positions (from PDB)
            surf: Surface scaling factor (default 0.7)
            resources: Compute resource ('CUDA' or 'CPU')
        """
        self.box_size = box_size
        self.surf = surf
        self.resources = resources.upper()
        
        # Validate resources
        if self.resources not in ['CUDA', 'CPU']:
            raise ValueError("Invalid resource specified. Use 'CUDA' or 'CPU'.")
        
        # Store topology information
        self.global_sequence = topology_info['global_sequence']
        self.chain_ids = topology_info['chain_ids']
        self.folded_domains = topology_info.get('folded_domains', [0] * len(self.global_sequence))
        self.component_names = topology_info.get('component_names', [''] * len(self.global_sequence))
        self.local_residue_indices = topology_info.get('local_residue_indices', list(range(1, len(self.global_sequence) + 1)))
        self.sasa_values = topology_info.get('sasa_values', [5.0] * len(self.global_sequence))
        
        # Store positions
        self.positions = positions
        
        # Derived properties
        self.n_residues = len(self.global_sequence)
        self.n_chains = max(self.chain_ids) if self.chain_ids else 1
        
        # Initialize force field parameters
        self._init_forcefield_params()
    
    def _init_forcefield_params(self):
        """Initialize force field parameters for all amino acids."""
        self.ff_param = {
            'ALA': {'mass': 71.079, 'charge': 0.0, 'radius': 0.2845, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 0.796},
            'ARG': {'mass': 157.197, 'charge': 1.0, 'radius': 0.3567, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 1.921},
            'ASN': {'mass': 114.104, 'charge': 0.0, 'radius': 0.3150, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 1.281},
            'ASP': {'mass': 114.080, 'charge': -1.0, 'radius': 0.3114, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 1.162},
            'CYS': {'mass': 103.139, 'charge': 0.0, 'radius': 0.3024, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 1.074},
            'GLN': {'mass': 128.131, 'charge': 0.0, 'radius': 0.3311, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 1.575},
            'GLU': {'mass': 128.107, 'charge': -1.0, 'radius': 0.3279, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 1.462},
            'GLY': {'mass': 57.052, 'charge': 0.0, 'radius': 0.2617, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 0.544},
            'HIS': {'mass': 137.142, 'charge': 0.0, 'radius': 0.3338, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 1.634},
            'ILE': {'mass': 113.160, 'charge': 0.0, 'radius': 0.3360, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 1.410},
            'LEU': {'mass': 113.160, 'charge': 0.0, 'radius': 0.3363, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 1.519},
            'LYS': {'mass': 129.183, 'charge': 1.0, 'radius': 0.3439, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 1.923},
            'MET': {'mass': 131.193, 'charge': 0.0, 'radius': 0.3381, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 1.620},
            'PHE': {'mass': 147.177, 'charge': 0.0, 'radius': 0.3556, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 1.869},
            'PRO': {'mass': 98.125, 'charge': 0.0, 'radius': 0.3187, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 0.974},
            'SER': {'mass': 87.078, 'charge': 0.0, 'radius': 0.2927, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 0.933},
            'THR': {'mass': 101.105, 'charge': 0.0, 'radius': 0.3108, 
                    'epsilon': self.EPS_POLAR, 'azero': self.AZERO_POLAR, 'surface': 1.128},
            'TRP': {'mass': 186.214, 'charge': 0.0, 'radius': 0.3754, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 2.227},
            'TYR': {'mass': 163.176, 'charge': 0.0, 'radius': 0.3611, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 2.018},
            'VAL': {'mass': 99.133, 'charge': 0.0, 'radius': 0.3205, 
                    'epsilon': self.EPS_NOPOL, 'azero': self.AZERO_HYDRO, 'surface': 1.232}
        }
    
    @staticmethod
    def surface_calc(a: float, surf: float) -> float:
        """Calculate surface scaling factor."""
        if surf > 0:
            return np.min([a, surf]) * 1 / surf
        else:
            return 1
    
    def _get_residue_name(self, seq_idx: int) -> str:
        """Get 3-letter residue name from sequence index."""
        aa_1letter = self.global_sequence[seq_idx]
        return self.AA_1TO3.get(aa_1letter, 'ALA')
    
    def _get_chain_letter(self, chain_id: int) -> str:
        """Convert chain ID (1-based int) to chain letter."""
        # Chain ID 1 -> 'A', 2 -> 'B', etc.
        return chr(ord('A') + (chain_id - 1) % 26)
    
    def _get_residues_by_chain(self) -> Dict[int, List[int]]:
        """
        Get list of residue (atom) indices for each chain.
        Returns global atom indices grouped by chain.
        """
        residues_by_chain = {}
        for i, chain_id in enumerate(self.chain_ids):
            if chain_id not in residues_by_chain:
                residues_by_chain[chain_id] = []
            residues_by_chain[chain_id].append(i)
        return residues_by_chain
    
    def _get_topology_atoms_by_chain(self, top: topology.Topology) -> Dict[int, List]:
        """
        Get list of Topology Atom objects grouped by chain.
        Uses chain.id from the topology.
        """
        atoms_by_chain = {}
        for atom in top.atoms():
            chain_id = atom.residue.chain.id
            if chain_id not in atoms_by_chain:
                atoms_by_chain[chain_id] = []
            atoms_by_chain[chain_id].append(atom)
        return atoms_by_chain
    
    def _get_folded_regions(self) -> List[tuple]:
        """Get list of (start, end) tuples for folded regions (1-based)."""
        regions = []
        n = len(self.folded_domains)
        i = 0
        while i < n:
            if self.folded_domains[i] == 1:
                start = i
                while i < n and self.folded_domains[i] == 1:
                    i += 1
                end = i - 1
                # Convert to 1-based for consistency with original code
                regions.append((start, end))
            else:
                i += 1
        return regions
    
    def _build_topology(self) -> topology.Topology:
        """Build molecular topology from global sequence and chain IDs."""
        top = topology.Topology()
        
        # Get residues grouped by chain
        residues_by_chain = self._get_residues_by_chain()
        
        # Track atoms by chain for bond creation
        atoms_by_chain = {}
        
        # Build topology for each chain
        for chain_id in sorted(residues_by_chain.keys()):
            chain_letter = self._get_chain_letter(chain_id)
            chain = top.addChain(chain_letter)
            
            chain_atoms = []
            for seq_idx in residues_by_chain[chain_id]:
                residue_name = self._get_residue_name(seq_idx)
                local_idx = self.local_residue_indices[seq_idx]
                residue = top.addResidue(residue_name, chain, id=str(local_idx))
                atom = top.addAtom('CA', element=element.carbon, residue=residue)
                chain_atoms.append(atom)
            
            atoms_by_chain[chain_id] = chain_atoms
            
            # Add bonds between consecutive atoms in the chain
            for i in range(len(chain_atoms) - 1):
                top.addBond(chain_atoms[i], chain_atoms[i + 1])
        
        return top
    
    def _create_box_vectors(self):
        """Create periodic box vectors."""
        a = Quantity(np.zeros([3]), nanometers)
        a[0] = self.box_size * nanometers
        b = unit.Quantity(np.zeros([3]), nanometers)
        b[1] = self.box_size * nanometers
        c = unit.Quantity(np.zeros([3]), nanometers)
        c[2] = self.box_size * nanometers
        return a, b, c
    
    def _add_bond_force(self, system, top):
        """Add harmonic bond force using sequence information."""
        f_bond = openmm.HarmonicBondForce()
        
        # Get atoms grouped by chain from topology
        atoms_by_chain = self._get_topology_atoms_by_chain(top)
        
        for chain_id in sorted(atoms_by_chain.keys()):
            chain_atoms = atoms_by_chain[chain_id]
            
            for i in range(len(chain_atoms) - 1):
                atom1 = chain_atoms[i]
                atom2 = chain_atoms[i + 1]
                
                res_name1 = self._get_residue_name(atom1.index)
                res_name2 = self._get_residue_name(atom2.index)
                
                if res_name1 in self.PROTEIN_RESIDUES:
                    l0 = self.L0_PRO
                elif res_name1 in self.RNA_RESIDUES:
                    l0 = self.L0_RNA
                else:
                    l0 = self.L0_PRO
                
                f_bond.addBond(atom1.index, atom2.index, l0 * nanometer, 
                              self.KBOND * kilojoules_per_mole / (nanometer**2))
        
        system.addForce(f_bond)
    
    def _add_angle_force(self, system, top):
        """Add harmonic angle force using chain IDs from topology info."""
        f_angle = openmm.HarmonicAngleForce()
        
        # Get atoms grouped by chain from topology
        atoms_by_chain = self._get_topology_atoms_by_chain(top)
        
        for chain_id in sorted(atoms_by_chain.keys()):
            chain_atoms = atoms_by_chain[chain_id]
            
            if chain_atoms:
                # Get residue name of first atom to determine force constant
                first_res_name = self._get_residue_name(chain_atoms[0].index)
                if first_res_name in self.PROTEIN_RESIDUES:
                    kangle = self.KANGLE_PRO
                elif first_res_name in self.RNA_RESIDUES:
                    kangle = self.KANGLE_RNA
                else:
                    kangle = self.KANGLE_PRO
                
                # Add angles for consecutive triplets
                for i in range(len(chain_atoms) - 2):
                    atom1 = chain_atoms[i]
                    atom2 = chain_atoms[i + 1]
                    atom3 = chain_atoms[i + 2]
                    
                    f_angle.addAngle(atom1.index, atom2.index, atom3.index,
                                    self.THETA0 * degrees, 
                                    kangle * kilojoule_per_mole / (radian**2))
        
        system.addForce(f_angle)
    
    def _add_electrostatic_force(self, system, top):
        """Add electrostatic force - match old implementation exactly."""
        k0 = self.KAPPA * nanometer
        
        equation = "S*(A+Z)/r*exp(-r/K0); "
        equation += "A=A1*A2; "
        equation += "Z=Z1+Z2; "
        equation += "S=(S1*S2)^(1/2)"
        
        force1 = CustomNonbondedForce(equation)
        force1.addGlobalParameter("K0", k0)
        force1.addPerParticleParameter("A")
        force1.addPerParticleParameter("Z")
        force1.addPerParticleParameter("S")
        force1.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        force1.setCutoffDistance(3.0 * nanometer)
        
        # Match old implementation: iterate over top.atoms()
        # Use global sequence index (atom.index) for sasa_values
        for atom in top.atoms():
            ff = self.ff_param[atom.residue.name]
            
            # Use actual SASA value from surface file (atom.index is global sequence index)
            sasa_value = self.sasa_values[atom.index]
            surface_factor = self.surface_calc(sasa_value / ff['surface'], self.surf)
            
            force1.addParticle([
                (np.sqrt(0.75 * np.abs(ff['charge'])) * 
                 np.sign(ff['charge'])) * nanometer * kilojoule / mole,
                ff['azero'] * (nanometer * kilojoule / mole) ** (1/2),
                surface_factor
            ])
        
        # Create exclusions for bonded pairs using topology atoms
        force1.createExclusionsFromBonds([(i[0].index, i[1].index) for i in top.bonds()], 1)
        
        system.addForce(force1)
    
    def _add_vdw_force(self, system, top):
        """Add van der Waals force - match old implementation exactly."""
        equation = "S*4*epsilon*((sigma/r)^10-(sigma/r)^5); "
        equation += "sigma=0.5*(sigma1+sigma2); "
        equation += "epsilon=sqrt(epsilon1*epsilon2); "
        equation += "S=(S1*S2)^(1/2)"
        
        force2 = CustomNonbondedForce(equation)
        force2.addPerParticleParameter("sigma")
        force2.addPerParticleParameter("epsilon")
        force2.addPerParticleParameter("S")
        force2.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        force2.setCutoffDistance(3.0 * nanometer)
        
        # Match old implementation: iterate over top.atoms()
        # Use global sequence index (atom.index) for sasa_values
        for atom in top.atoms():
            ff = self.ff_param[atom.residue.name]
            
            # Use actual SASA value from surface file (atom.index is global sequence index)
            sasa_value = self.sasa_values[atom.index]
            surface_factor = self.surface_calc(sasa_value / ff['surface'], self.surf)
            
            force2.addParticle([
                ff['radius'] * 2 * 2**(-1/6) * nanometer,
                ff['epsilon'] * kilojoule / mole,
                surface_factor
            ])
        
        # Create exclusions for bonded pairs
        force2.createExclusionsFromBonds([(i[0].index, i[1].index) for i in top.bonds()], 1)
        
        system.addForce(force2)
    
    def _add_special_forces(self, system, top):
        """Add special interactions (cation-pi, pi-pi) - match old implementation exactly."""
        # Check for required amino acids using topology
        residue_names = [atom.residue.name for atom in top.atoms()]
        has_cation = any(name in ['ARG', 'LYS'] for name in residue_names)
        has_aromatic = any(name in ['PHE', 'TRP', 'TYR'] for name in residue_names)
        
        # Get indices for each amino acid type using global sequence index (atom.index)
        arg_lys_indices = [atom.index for atom in top.atoms() if atom.residue.name in ['ARG', 'LYS']]
        aromatic_indices = [atom.index for atom in top.atoms() if atom.residue.name in ['PHE', 'TRP', 'TYR']]
        
        # Cation-pi interaction
        if has_cation and has_aromatic:
            equation = "S*4*epsilon*((sigma/r)^10-(sigma/r)^5); "
            equation += "sigma=0.5*(sigma1+sigma2); "
            equation += "epsilon=sqrt(epsilon1*epsilon2); "
            equation += "S=(S1*S2)^(1/2)"
            
            force3 = CustomNonbondedForce(equation)
            force3.addPerParticleParameter("sigma")
            force3.addPerParticleParameter("epsilon")
            force3.addPerParticleParameter("S")
            force3.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
            force3.setCutoffDistance(3.0 * nanometer)
            
            # Match old implementation: iterate over top.atoms()
            for atom in top.atoms():
                ff = self.ff_param[atom.residue.name]
                
                # Use actual SASA value from surface file (atom.index is global sequence index)
                sasa_value = self.sasa_values[atom.index]
                surface_factor = self.surface_calc(sasa_value / ff['surface'], self.surf)
                
                force3.addParticle([
                    ff['radius'] * 2 * 2**(-1/6) * nanometer,
                    self.CATIONPI_PROPRO * kilojoule / mole,
                    surface_factor
                ])
            
            # Create exclusions and add interaction group
            force3.createExclusionsFromBonds([(i[0].index, i[1].index) for i in top.bonds()], 1)
            force3.addInteractionGroup(arg_lys_indices, aromatic_indices)
            
            # Add force to system
            system.addForce(force3)
        
        # Pi-pi interaction
        if has_aromatic:
            equation = "S*4*epsilon*((sigma/r)^10-(sigma/r)^5); "
            equation += "sigma=0.5*(sigma1+sigma2); "
            equation += "epsilon=sqrt(epsilon1*epsilon2); "
            equation += "S=(S1*S2)^(1/2)"
            
            force5 = CustomNonbondedForce(equation)
            force5.addPerParticleParameter("sigma")
            force5.addPerParticleParameter("epsilon")
            force5.addPerParticleParameter("S")
            force5.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
            force5.setCutoffDistance(3.0 * nanometer)
            
            # Match old implementation: iterate over top.atoms()
            for atom in top.atoms():
                ff = self.ff_param[atom.residue.name]
                
                # Use actual SASA value from surface file (atom.index is global sequence index)
                sasa_value = self.sasa_values[atom.index]
                surface_factor = self.surface_calc(sasa_value / ff['surface'], self.surf)
                
                force5.addParticle([
                    ff['radius'] * 2 * 2**(-1/6) * nanometer,
                    self.PIPI_PROPRO * kilojoule / mole,
                    surface_factor
                ])
            
            # Create exclusions
            force5.createExclusionsFromBonds([(i[0].index, i[1].index) for i in top.bonds()], 1)
            
            force5.addInteractionGroup(aromatic_indices, aromatic_indices)
            system.addForce(force5)
    
    def _add_elastic_network(self, system, top):
        """Add elastic network for folded domains using folded_domains."""
        # Get folded regions
        folded_regions = self._get_folded_regions()
        
        if not folded_regions:
            print("No folded regions found, skipping elastic network")
            return  # No elastic network needed
        
        elastic_force = CustomBondForce("0.5*k*(r-r0)^2")
        elastic_force.addGlobalParameter("k", self.FORCE_CONSTANT / nanometer ** 2)
        elastic_force.addPerBondParameter("r0")
        
        # Build elastic network for each folded region in each chain
        atoms_by_chain = self._get_topology_atoms_by_chain(top)
        
        # Get residues grouped by chain (using global sequence indices)
        residues_by_chain = self._get_residues_by_chain()
        
        for chain_id in sorted(atoms_by_chain.keys()):
            chain_atoms = atoms_by_chain[chain_id]
            # Convert chain letter to integer index for residues_by_chain
            # chain_id is like 'A', 'B', etc. from topology
            # residues_by_chain uses integers 1, 2, 3, etc.
            chain_int_id = ord(chain_id) - ord('A') + 1
            
            # Get global residue indices for this chain
            chain_global_indices = residues_by_chain.get(chain_int_id, [])
            
            if not chain_global_indices:
                print(f"  Warning: No global indices found for chain {chain_id} (int_id={chain_int_id})")
                continue
            
            for region_start, region_end in folded_regions:
                # Get global indices in this region that belong to this chain
                # region_start and region_end are 0-based global indices
                region_indices = [idx for idx in chain_global_indices 
                                if region_start <= idx <= region_end]
                
                print(f"  Chain {chain_id}: region [{region_start+1}, {region_end+1}] has {len(region_indices)} residues in chain")
                
                # Add bonds between all pairs within cutoff distance
                n_region = len(region_indices)
                bonds_added = 0
                for i in range(n_region):
                    for j in range(i + 3, n_region):  # i+3 to skip nearby residues
                        idx1 = region_indices[i]
                        idx2 = region_indices[j]
                        
                        # Use positions from self.positions (global indices)
                        pos1 = self.positions[idx1]
                        pos2 = self.positions[idx2]
                        distance = unit.norm(pos1 - pos2)
                        
                        if distance < self.CUTOFF_DISTANCE:
                            elastic_force.addBond(idx1, idx2, [distance])
                            bonds_added += 1
                
                print(f"  Added {bonds_added} elastic bonds in region [{region_start+1}, {region_end+1}]")
            
            print(f"Building ENM: chain {chain_id}/{self.n_chains}")
        
        system.addForce(elastic_force)
    
    def create_system(self) -> openmm.System:
        """
        Create and return OpenMM System.
        
        Returns:
            openmm.System: Complete OpenMM system with all forces
        """
        system = openmm.System()
        
        # Build topology
        top = self._build_topology()
        
        # Add particles using sequence information
        for i in range(self.n_residues):
            res_name = self._get_residue_name(i)
            mass = self.ff_param[res_name]['mass']
            system.addParticle(mass * amu)
        
        # Set box vectors
        a, b, c = self._create_box_vectors()
        system.setDefaultPeriodicBoxVectors(a, b, c)
        
        # Add all forces (CMMotionRemover added early to match OLD COCOMO order)
        self._add_bond_force(system, top)
        self._add_angle_force(system, top)
        
        # Add CMotion remover (right after bond/angle forces, matching OLD)
        system.addForce(openmm.CMMotionRemover())
        
        self._add_electrostatic_force(system, top)
        self._add_vdw_force(system, top)
        self._add_special_forces(system, top)
        
        # Add elastic network if there are folded domains
        if any(self.folded_domains):
            self._add_elastic_network(system, top)
        
        return system


def main():
    """Main function to run COCOMO simulation from config file."""
    # Get the directory of the config file
    config_path = os.path.join(os.path.dirname(__file__), 'config_idp.yaml')
    
    # Load configuration from YAML
    config = CGSimulationConfig.from_yaml(config_path)
    
    # Create simulator and setup output
    sim = CGSimulator(config)
    sim.setup(config.output_dir)
    
    # ============================================================
    # 获取体系拓扑信息 (使用新增的接口方法)
    # ============================================================
    print("\n" + "=" * 60)
    print("体系拓扑信息")
    print("=" * 60)
    
    # 1. 获取体系构成
    print("\n1. 体系构成 (get_composition):")
    composition = sim.get_composition()
    for comp in composition:
        print(f"   - {comp['name']}: nmol={comp['nmol']}, nres={comp['nres']}, type={comp['type']}")
    
    # 2. 获取全局序列
    print("\n2. 全局序列 (get_global_sequence):")
    global_seq = sim.get_global_sequence()
    print(f"   - 长度: {len(global_seq)}")
    print(f"   - 前50字符: {global_seq[:50]}...")
    print(f"   - 后50字符: ...{global_seq[-50:]}")
    
    # 3. 获取链ID
    print("\n3. 链ID (get_chain_ids):")
    chain_ids = sim.get_chain_ids()
    print(f"   - 长度: {len(chain_ids)}")
    print(f"   - 链ID范围: 1-{max(chain_ids)} (共 {max(chain_ids)} 条链)")
    
    # 4. 获取fold domain信息
    print("\n4. Folded Domain (get_folded_domains):")
    folded = sim.get_folded_domains()
    n_folded = sum(folded)
    n_unfolded = len(folded) - n_folded
    print(f"   - 长度: {len(folded)}")
    print(f"   - Folded区域: {n_folded} ({n_folded/len(folded)*100:.1f}%)")
    print(f"   - Unfolded区域: {n_unfolded} ({n_unfolded/len(folded)*100:.1f}%)")
    
    # 5. 构建 topology_info 字典
    print("\n5. 构建 COCOMO 系统:")
    topology_info = {
        'global_sequence': sim.get_global_sequence(),
        'chain_ids': sim.get_chain_ids(),
        'folded_domains': sim.get_folded_domains(),
        'component_names': [],  # 可以从 get_composition 构建
        'local_residue_indices': list(range(1, len(sim.get_global_sequence()) + 1))
    }
    
    # 从 composition 构建 component_names
    comp_names = []
    for comp in sim.get_composition():
        for _ in range(comp['nmol']):
            comp_names.extend([comp['name']] * comp['nres'])
    topology_info['component_names'] = comp_names
    
    print(f"   - topology_info 构建完成")
    print(f"   - global_sequence 长度: {len(topology_info['global_sequence'])}")
    print(f"   - chain_ids 长度: {len(topology_info['chain_ids'])}")
    
    print("=" * 60 + "\n")


if __name__ == '__main__':
    main()
