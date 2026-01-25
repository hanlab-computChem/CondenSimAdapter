"""
Biomolecule classes for OpenMpipi.

Provides classes for coarse-grained biomolecule representations:
- CGBiomolecule: Base class
- IDP: Intrinsically Disordered Protein
- MDP: Multi-domain (globular) Protein
- RNA: Single-stranded RNA
"""

import numpy as np
import openmm.app as app
import openmm as mm
import openmm.unit as unit
import mdtraj as md
import os
import hashlib
import json
import pickle

from .coordinate_building import generate_spiral_coords, parse_pdb
from .system_building import get_mpipi_system
from .constants import PLATFORM, PROPERTIES


class CGBiomolecule:
    """
    Base class representing a coarse-grained biomolecule. 

    This class provides the basic structure for a biomolecule, including methods
    for sequence validation, topology creation, and system setup in OpenMM. 
    Subclasses must implement methods for generating specific topologies and initial coordinates.

    Attributes:
        chain_id (str): Identifier for the biomolecular chain.
        sequence (str): Sequence of the biomolecule (protein or DNA).
        globular_indices (list): List of residue indices that form the globular (structured) part of the molecule.
        dyad_positions (list, optional): Positions of dyads in the molecule (used for nucleosome/DNA structures).

        topology (Topology): OpenMM Topology object for the biomolecule.
        initial_coords (ndarray): Initial 3D coordinates for the biomolecule.
        
        chain_mass (Quantity): Total mass of the chain (in daltons or equivalent).
        min_rg_coords (ndarray): Coordinates of the chain's minimum Rg conformation (after relaxation).
        min_rg (float): Minimum radius of gyration observed during relaxation.
        max_rg (float): Maximum radius of gyration observed during relaxation.
    """

    def __init__(self, chain_id, sequence, valid_residues='', globular_indices=[], dyad_positions=None):
        """
        Initializes a CGBiomolecule object with the given chain ID, sequence, and optional attributes.

        Args:
            chain_id (str): Identifier for the biomolecular chain.
            sequence (str): Sequence of the biomolecule (amino acids or nucleotides).
            valid_residues (str, optional): A string containing valid residues. Defaults to an empty string.
            globular_indices (list, optional): Indices of residues that form the globular region. Defaults to an empty list.
            dyad_positions (list, optional): Positions of dyads in the molecule (used for DNA structures). Defaults to None.
        """
        self.chain_id = chain_id
        self.sequence = sequence
        self.globular_indices = globular_indices
        self.validate_sequence(set(valid_residues))  # Validate the sequence with provided valid residues
        self.dyad_positions = dyad_positions
        
        # Attributes set later (after building topology and relaxation):
        self.topology = None
        self.initial_coords = None
        self.chain_mass = None
        self.min_rg_coords = None
        self.min_rg = None
        self.max_rg = None

    def validate_sequence(self, valid_entries):
        """
        Validates the biomolecule's sequence to ensure it contains only valid residue entries.

        Args:
            valid_entries (set): A set of valid residue or nucleotide symbols.

        Raises:
            ValueError: If the sequence contains symbols not present in the set of valid entries.
        """
        if not set(self.sequence).issubset(valid_entries):
            raise ValueError(f"Invalid sequence. Sequence must only contain valid entries: {valid_entries}")

    def create_monomer_topology(self):
        """
        Creates the topology for the biomolecule. 
        
        This is an abstract method that must be implemented by subclasses to define the
        specific topology for the biomolecule type.

        Returns:
            Topology: An OpenMM Topology object for the biomolecule.
        """
        raise NotImplementedError("Subclasses must implement create_monomer_topology method")

    def generate_initial_coords(self):
        """
        Generates the initial coordinates for the biomolecule.

        This is an abstract method that must be implemented by subclasses to define
        the specific coordinate generation strategy for the biomolecule type.

        Returns:
            ndarray: Initial 3D coordinates for the biomolecule.
        """
        raise NotImplementedError("Subclasses must implement generate_initial_coords method")

    def create_system(self, T=280*unit.kelvin, csx=150):
        """
        Creates an OpenMM System for the monomer using get_mpipi_system, 
        based on the current topology, initial_coords, and globular indices.

        Args:
            T (Quantity, optional): Temperature in Kelvin (default 280 K).
            csx (float, optional): Ionic strength in mM (default 150).

        Returns:
            System: The OpenMM System representing this biomolecule.
        """
        # For system building, we need positions as a NumPy array
        # and a dict of {chain_id: globular_indices}.
        glob_dict = {self.chain_id: self.globular_indices}
        
        system = get_mpipi_system(
            np.array(self.initial_coords),    # positions
            self.topology,                    # topology
            glob_dict, 
            T.value_in_unit(unit.kelvin),     # pass float, no units
            csx,
            CM_remover=False, 
            periodic=False
        )
        
        # Store chain mass (sum of all particle masses in the system)
        self.chain_mass = sum([system.getParticleMass(i) for i in range(system.getNumParticles())],
                              0*unit.dalton)
        return system

    def _get_cache_filename(self, simulation_time, T, csx, cache_dir='.relaxed_monomers'):
        """
        Generate a cache filename based on molecule properties and simulation parameters.
        
        Args:
            simulation_time (Quantity): Duration of the simulation.
            T (Quantity): Temperature.
            csx (float): Ionic strength in mM.
            cache_dir (str): Directory to store cache files.
        
        Returns:
            tuple: (cache_dir, cache_file) paths
        """
        # Create a unique identifier from sequence, globular_indices, and simulation params
        cache_key = {
            'chain_id': self.chain_id,
            'sequence': self.sequence,
            'globular_indices': self.globular_indices if hasattr(self, 'globular_indices') else None,
            'simulation_time_ns': simulation_time.value_in_unit(unit.nanosecond),
            'T_kelvin': T.value_in_unit(unit.kelvin),
            'csx_mM': csx
        }
        
        # Create hash from cache key
        cache_str = json.dumps(cache_key, sort_keys=True)
        cache_hash = hashlib.md5(cache_str.encode()).hexdigest()[:12]
        
        cache_filename = f"{self.chain_id}_{cache_hash}.pkl"
        
        # Create cache directory if it doesn't exist
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        
        return cache_dir, os.path.join(cache_dir, cache_filename)
    
    def get_compact_model(self, 
                          simulation_time=5*unit.nanosecond, 
                          T=280*unit.kelvin, 
                          csx=150,
                          use_cache=True,
                          cache_dir='.relaxed_monomers'):
        """
        Relax (compactify) this biomolecule by simulating it in isolation. 

        1. Builds a System via `create_system`.
        2. Minimizes, then runs a short MD simulation at temperature T.
        3. Tracks the radius of gyration (Rg), storing coords of the minimum-Rg conformation.

        Args:
            simulation_time (Quantity, optional): Duration of the simulation (default 5 ns).
            T (Quantity, optional): Temperature (default 280 K).
            csx (float, optional): Ionic strength in mM (default 150).
            use_cache (bool, optional): Whether to use cached results if available (default True).
            cache_dir (str, optional): Directory to store cache files (default '.relaxed_monomers').
        """
        # Check cache first
        if use_cache:
            cache_dir_path, cache_file = self._get_cache_filename(simulation_time, T, csx, cache_dir)
            if os.path.exists(cache_file):
                print(f"  Loading cached relaxed structure from {cache_file}...", flush=True)
                try:
                    with open(cache_file, 'rb') as f:
                        cache_data = pickle.load(f)
                    self.min_rg_coords = cache_data['min_rg_coords']
                    self.min_rg = cache_data['min_rg']
                    self.max_rg = cache_data['max_rg']
                    print(f"  ✓ Loaded cached structure (Rg: {self.min_rg:.3f} nm)", flush=True)
                    return
                except Exception as e:
                    print(f"  Warning: Failed to load cache ({e}), regenerating...", flush=True)
        
        # Build system
        system = self.create_system(T, csx)

        # Create integrator
        integrator = mm.LangevinMiddleIntegrator(T, 0.01/unit.picosecond, 10*unit.femtosecond)
        
        # Setup Simulation
        simulation = app.Simulation(
            self.topology, 
            system, 
            integrator, 
            platform=PLATFORM, 
            platformProperties=PROPERTIES
        )
        simulation.context.setPositions(self.initial_coords)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(T)

        # Write out a temporary trajectory so we can compute Rg
        temp_traj_file = f".temp_{self.chain_id}_Rg_traj.pdb"
        simulation.reporters.append(app.PDBReporter(temp_traj_file, 5000, enforcePeriodicBox=False))

        # Run short simulation
        simulation_steps = int(simulation_time / (10*unit.femtosecond))
        simulation.step(simulation_steps)

        # Load the temporary trajectory in MDTraj, compute Rg
        traj = md.load_pdb(temp_traj_file)
        rg_values = md.compute_rg(traj)

        # Identify the frame with minimum Rg
        min_rg_index = np.argmin(rg_values)
        self.min_rg_coords = traj.xyz[min_rg_index]
        self.min_rg = rg_values[min_rg_index]
        self.max_rg = np.max(rg_values)

        # Save to cache
        if use_cache:
            cache_dir_path, cache_file = self._get_cache_filename(simulation_time, T, csx, cache_dir)
            try:
                cache_data = {
                    'min_rg_coords': self.min_rg_coords,
                    'min_rg': self.min_rg,
                    'max_rg': self.max_rg
                }
                with open(cache_file, 'wb') as f:
                    pickle.dump(cache_data, f)
                print(f"  ✓ Saved relaxed structure to cache: {cache_file}", flush=True)
            except Exception as e:
                print(f"  Warning: Failed to save cache ({e})", flush=True)

        # Clean up
        os.remove(temp_traj_file)


class IDP(CGBiomolecule):
    """
    Class representing an intrinsically disordered protein (IDP).

    Inherits from CGBiomolecule and implements methods specific to IDPs,
    including topology creation and coordinate generation. IDPs are
    modelled as fully-flexible polymers.
    """

    def __init__(self, chain_id, sequence):
        """
        Initializes an IDP object with the given chain ID and sequence.

        Args:
            chain_id (str): Identifier for the protein chain.
            sequence (str): Amino acid sequence of the intrinsically disordered protein.
        """
        super().__init__(chain_id, sequence, valid_residues='ACDEFGHIKLMNPQRSTVWY')
        
        # Build topology and initial coords right away
        self.topology = self.create_monomer_topology()
        self.initial_coords = self.generate_initial_coords()

    def create_monomer_topology(self):
        """
        Creates the topology for the IDP using the provided sequence and chain ID.

        Returns:
            Topology: The OpenMM Topology object for the IDP.
        """
        return create_monomer_topology(self.sequence, self.chain_id, chain_type='prt')

    def generate_initial_coords(self):
        """
        Generates initial coordinates for the IDP, modeled as a spiral structure.

        The coordinates are spaced using a default spacing value of 0.381 nm
        between each residue.

        Returns:
            ndarray: Initial 3D coordinates for the IDP in a spiral configuration.
        """
        return generate_spiral_coords(len(self.sequence), spacing=0.381)


class MDP(CGBiomolecule):
    """
    Class representing a multi-domain (globular) protein (MDP). Globular domains 
    are constrained using an elastic-network model, whilst intrinsically-disordered 
    regions are modelled as fully-flexible polymers.

    Inherits from CGBiomolecule and implements methods specific to globular proteins,
    including topology creation and coordinate generation from a PDB file.
    """

    def __init__(self, chain_id, sequence, globular_indices, pdb_file):
        """
        Initializes an MDP object with the given chain ID, sequence, globular indices, and PDB file.

        Args:
            chain_id (str): Identifier for the protein chain.
            sequence (str): Amino acid sequence of the multi-domain protein.
            globular_indices (list): List of residue indices corresponding to globular (structured) regions.
            pdb_file (str): Path to the PDB file containing the initial structure of the protein.
        """
        self.pdb_file = pdb_file
        super().__init__(chain_id, sequence, 
                         valid_residues='ACDEFGHIKLMNPQRSTVWY', 
                         globular_indices=globular_indices)

        # Build topology and initial coords right away
        self.topology = self.create_monomer_topology()
        self.initial_coords = self.generate_initial_coords()

    def create_monomer_topology(self):
        """
        Creates the topology for the multi-domain protein using the provided sequence,
        chain ID, and globular indices.

        Returns:
            Topology: The OpenMM Topology object for the multi-domain protein.
        """
        return create_monomer_topology(
            self.sequence, 
            self.chain_id, 
            chain_type='prt', 
            globular_indices=self.globular_indices
        )

    def generate_initial_coords(self):
        """
        Parses the PDB file to generate initial coordinates for the multi-domain protein.

        If the sequence in the PDB file does not match the provided sequence, a warning is raised.

        Returns:
            ndarray: Initial 3D coordinates for the multi-domain protein from the PDB file.

        Raises:
            Warning: If there is a mismatch between the provided sequence and the sequence from the PDB file.
        """
        coords, pdb_sequence = parse_pdb(self.pdb_file, self.globular_indices)
        if self.sequence != pdb_sequence:
            print(self.sequence, pdb_sequence)
            raise Warning("Warning! Mismatch between the provided sequence and the sequence read from the PDB file.")
        return coords


class RNA(CGBiomolecule):
    """
    Class representing an unstructured single-stranded RNA molecule.

    Inherits from CGBiomolecule and implements methods specific to RNAs,
    including topology creation and coordinate generation. RNAs are
    modelled as fully-flexible polymers.

    At present, only parameters for U are available so only polyU sequences are valid.
    """

    def __init__(self, chain_id, sequence):
        """
        Initializes an RNA object with the given chain ID and sequence.

        Args:
            chain_id (str): Identifier for the RNA chain.
            sequence (str): Nucleotide sequence of the RNA.
        """
        super().__init__(chain_id, sequence, valid_residues='U')
        
        # Build topology and initial coords right away
        self.topology = self.create_monomer_topology()
        self.initial_coords = self.generate_initial_coords()

    def create_monomer_topology(self):
        """
        Creates the topology for the RNA using the provided sequence and chain ID.

        Returns:
            Topology: The OpenMM Topology object for the RNA.
        """
        return create_monomer_topology(self.sequence, self.chain_id, chain_type='RNA')

    def generate_initial_coords(self):
        """
        Generates initial coordinates for the RNA, modeled as a spiral structure.

        The coordinates are spaced using a default spacing value of 0.50 nm
        between each residue.

        Returns:
            ndarray: Initial 3D coordinates for the RNA in a spiral configuration.
        """
        return generate_spiral_coords(len(self.sequence), spacing=0.50)


def create_monomer_topology(sequence, chain_id, chain_type, globular_indices=None, topology=None):
    """
    Creates an OpenMM Topology object for the given biomolecular sequence.

    This function builds the topology of a biomolecule by adding residues 
    and atoms to the topology. It can either create a new Topology object or add to an existing one.

    Args:
        sequence (str): The sequence of the biomolecule (e.g., protein or nucleic acids).
        chain_id (str): Identifier for the biomolecular chain in the topology.
        chain_type (str): Type of the chain ('prt' for protein, 'RNA' for RNA).
        globular_indices (list, optional): List of indices that represent globular regions of the sequence. Defaults to None.
        topology (Topology, optional): An existing OpenMM Topology object to add the chain to. If not provided, a new Topology will be created.

    Returns:
        Topology: The OpenMM Topology object representing the biomolecule.
    """
    chain_dict = {
        'prt': ['Cu', 'p'],   # "Cu" for coil/unstructured, "p" prefix for IDP/MDP
        'RNA': ['Au', 'r']    # "Au" for single-stranded RNA, "r" prefix
    }

    element_symbol = chain_dict[chain_type][0]
    prefix = chain_dict[chain_type][1]

    if not topology:
        topology = app.Topology()

    chain = topology.addChain(id=chain_id)

    if globular_indices is None:
        globular_indices = []
    
    # Flatten nested lists (if needed)
    flattened_globular_indices = [index for sublist in globular_indices for index in sublist] \
                                 if any(isinstance(el, list) for el in globular_indices) \
                                 else globular_indices

    for res_id, residue_name in enumerate(sequence):
        # If this residue is in the globular domain, use 'Pt' (platinum) to label it
        # (matching your previous code's usage)
        if res_id in flattened_globular_indices:
            residue = topology.addResidue(prefix + residue_name, chain)
            topology.addAtom(prefix + residue_name, app.Element.getBySymbol('Pt'), residue)
        else:
            residue = topology.addResidue(prefix + residue_name, chain)
            topology.addAtom(prefix + residue_name, app.Element.getBySymbol(element_symbol), residue)

    return topology

