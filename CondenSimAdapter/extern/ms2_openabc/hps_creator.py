"""
HPS System Builder for CGSimulator

This module provides functions to build HPS (Hierarchical Physics-based Scaling) 
coarse-grained protein systems for integration with CGSimulator.

Design follows the same pattern as COCOMO and Mpipi-Recharged:
1. Pre-equilibration (using CALVADOS)
2. Build OpenMM System
3. Run Simulation
4. Post-processing
"""

import numpy as np
import pandas as pd
import os
from typing import Dict, List, Optional, Any

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit

from openmm import Platform, LangevinMiddleIntegrator


def get_hps_system(
    positions: np.ndarray,
    topology: app.Topology,
    folded_domains: List[int],
    global_sequence: str,
    chain_ids: List[int],
    box_size: Optional[List[float]] = None,
    use_pbc: bool = True,
    **kwargs
) -> mm.System:
    """
    Build HPS OpenMM System from positions and topology.
    
    This function follows the same pattern as get_mpipi_system() and COCOMO class.
    
    Parameters
    ----------
    positions : np.ndarray
        Atom positions from pre-equilibrated structure.
    
    topology : app.Topology
        OpenMM Topology object.
    
    folded_domains : List[int]
        List of fold domain flags (0/1) for each residue.
        0 = IDP/unfolded region, 1 = folded domain.
    
    global_sequence : str
        Global sequence string (1-letter amino acid codes).
    
    chain_ids : List[int]
        List of chain IDs for each residue (1-based).
    
    box_size : List[float], optional
        Box dimensions [Lx, Ly, Lz] in nm.
        Required if use_pbc=True.
    
    use_pbc : bool
        Whether to use periodic boundary conditions. Default True.
    
    **kwargs : dict
        Additional parameters:
        - hydropathy_scale: 'Urry' or 'KR' (default 'Urry')
        - mu: Hydropathy scale factor (default 1)
        - delta: Hydropathy shift factor (default 0.08)
        - epsilon: Contact strength (default 0.2 kcal/mol)
        - ldby: Debye length in nm (default 1.0)
        - dielectric_water: Water dielectric constant (default 80.0)
        - cutoff: Electrostatic cutoff in nm (default 3.5)
        - en_force_constant: Elastic network force constant (default 700 kJ/mol/nm^2)
    
    Returns
    -------
    system : mm.System
        Complete OpenMM System with all HPS forces.
    
    Notes
    -----
    This function builds the system by:
    1. Building HPSParser objects for each chain
    2. Aggregating into HPSModel
    3. Adding all forces (bonds, contacts, electrostatics, elastic network)
    """
    from CondenSimAdapter.extern.ms2_openabc.forcefields.parsers.hps_parser import HPSParser
    from CondenSimAdapter.extern.ms2_openabc.forcefields.hps_model import HPSModel
    from CondenSimAdapter.extern.ms2_openabc.utils.helper_functions import parse_pdb
    from CondenSimAdapter.extern.ms2_openabc.lib import _amino_acids, _kcal_to_kj
    
    # Default parameters
    hydropathy_scale = kwargs.get('hydropathy_scale', 'Urry')
    mu = kwargs.get('mu', 1)
    delta = kwargs.get('delta', 0.08)
    epsilon = kwargs.get('epsilon', 0.2*_kcal_to_kj)
    ldby = kwargs.get('ldby', 1.0*unit.nanometer)
    dielectric_water = kwargs.get('dielectric_water', 80.0)
    cutoff = kwargs.get('cutoff', 3.5*unit.nanometer)
    en_force_constant = kwargs.get('en_force_constant', 700.0*unit.kilojoule_per_mole/unit.nanometer**2)
    
    # Build HPS system
    model = HPSModel()
    
    # Group residues by chain
    n_residues = len(global_sequence)
    residues_by_chain = {}
    for i, chain_id in enumerate(chain_ids):
        if chain_id not in residues_by_chain:
            residues_by_chain[chain_id] = []
        residues_by_chain[chain_id].append(i)
    
    # Process each chain
    for chain_id in sorted(residues_by_chain.keys()):
        chain_residues = residues_by_chain[chain_id]
        
        # Build chain sequence
        chain_sequence = ''.join([global_sequence[i] for i in chain_residues])
        
        # Check if chain has folded domains
        chain_folded = [folded_domains[i] for i in chain_residues]
        has_folded = any(chain_folded)
        
        # Create temporary CA PDB for this chain
        from CondenSimAdapter.extern.ms2_openabc.utils.helper_functions import build_straight_CA_chain, write_pdb
        
        chain_ca_pdb = f'_temp_chain_{chain_id}.pdb'
        ca_atoms = build_straight_CA_chain(chain_sequence, r0=0.38)
        
        # Set positions from reference
        for j, res_idx in enumerate(chain_residues):
            ca_atoms.loc[j, 'x'] = positions[res_idx, 0] * 10  # nm to Angstrom
            ca_atoms.loc[j, 'y'] = positions[res_idx, 1] * 10
            ca_atoms.loc[j, 'z'] = positions[res_idx, 2] * 10
        
        write_pdb(ca_atoms, chain_ca_pdb)
        
        # Create parser
        if has_folded:
            # Find folded domains for this chain
            # Convert global indices to chain-local indices
            chain_folded_ranges = []
            i = 0
            n = len(chain_folded)
            while i < n:
                if chain_folded[i] == 1:
                    start = i
                    while i < n and chain_folded[i] == 1:
                        i += 1
                    end = i - 1
                    # Convert to 1-based for domain file
                    chain_folded_ranges.append([start + 1, end + 1])
                else:
                    i += 1
            
            # Create domain YAML for this chain
            chain_name = f'chain_{chain_id}'
            domains_yaml = f'_temp_{chain_name}_domains.yaml'
            domains_content = f'{chain_name}:\n'
            for start, end in chain_folded_ranges:
                domains_content += f'  - [{start}, {end}]\n'
            
            with open(domains_yaml, 'w') as f:
                f.write(domains_content)
            
            parser = HPSParser(chain_ca_pdb, fdomains=domains_yaml)
            
            # Cleanup domain file
            os.remove(domains_yaml)
        else:
            parser = HPSParser(chain_ca_pdb)
        
        # Append to model
        model.append_mol(parser)
        
        # Cleanup chain PDB
        os.remove(chain_ca_pdb)
    
    # Create system
    model.create_system(
        top=topology,
        use_pbc=use_pbc,
        box_a=box_size[0] if box_size else 50,
        box_b=box_size[1] if box_size else 50,
        box_c=box_size[2] if box_size else 50
    )
    
    # Add forces
    model.add_protein_bonds(force_group=1)
    model.add_contacts(
        hydropathy_scale=hydropathy_scale,
        epsilon=epsilon,
        mu=mu,
        delta=delta,
        force_group=2
    )
    model.add_dh_elec(
        ldby=ldby,
        dielectric_water=dielectric_water,
        cutoff=cutoff,
        force_group=3
    )
    model.add_elastic_network(
        force_constant=en_force_constant,
        force_group=4
    )
    
    return model.system


def build_hps_topology_info(
    global_sequence: str,
    chain_ids: List[int],
    folded_domains: List[int],
    component_names: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Build topology_info dictionary for HPS system.
    
    This is similar to what COCOMO uses.
    
    Parameters
    ----------
    global_sequence : str
        Global sequence string.
    
    chain_ids : List[int]
        List of chain IDs for each residue.
    
    folded_domains : List[int]
        List of fold domain flags.
    
    component_names : List[str], optional
        List of component names for each residue.
    
    Returns
    -------
    topology_info : Dict
        Dictionary containing topology information.
    
    """
    if component_names is None:
        # Build component names from chain_ids
        component_names = []
        current_chain = None
        for chain_id in chain_ids:
            if chain_id != current_chain:
                current_chain = chain_id
                comp_idx = len(set(chain_ids[:chain_ids.index(chain_id) + 1]))
                component_names.append(f'component_{comp_idx}')
            else:
                component_names.append(component_names[-1])
    
    topology_info = {
        'global_sequence': global_sequence,
        'chain_ids': chain_ids,
        'folded_domains': folded_domains,
        'component_names': component_names,
        'local_residue_indices': list(range(1, len(global_sequence) + 1)),
    }
    
    return topology_info

