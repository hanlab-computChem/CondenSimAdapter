"""
ms2_OpenMpipi - Internalized OpenMpipi coarse-grained simulation package

Mpipi-Recharged forcefield implementation for coarse-grained simulations
of intrinsically disordered proteins (IDPs) and multi-domain proteins (MDPs).

Modules:
- biomolecules: Classes for IDP, MDP, and RNA representations
- system_building: Core system building functions (get_mpipi_system)
- coordinate_building: Initial coordinate generation tools
- model_building: Slab configuration and equilibration tools

Author: Kieran Russell (kor20@cam.ac.uk)
Collepardo Lab, University of Cambridge
"""

# Core modules
from .mpipi_recharged import biomolecules
from .mpipi_recharged import system_building
from .mpipi_recharged import coordinate_building
from .mpipi_recharged import model_building
from .mpipi_recharged import constants

# Expose key classes and functions
from .mpipi_recharged.biomolecules import IDP, MDP, RNA, CGBiomolecule
from .mpipi_recharged.system_building import get_mpipi_system
from .mpipi_recharged.model_building import (
    build_mpipi_recharged_system_from_chains,
    calculate_target_box_vectors,
    build_model,
    equilibrate_slab,
    build_and_equilibrate_model,
)
from .mpipi_recharged.constants import PLATFORM, PROPERTIES, DATA_DIR

__all__ = [
    'biomolecules',
    'system_building',
    'coordinate_building',
    'model_building',
    'constants',
    'IDP',
    'MDP',
    'RNA',
    'CGBiomolecule',
    'get_mpipi_system',
    'build_mpipi_recharged_system_from_chains',
    'calculate_target_box_vectors',
    'build_model',
    'equilibrate_slab',
    'build_and_equilibrate_model',
    'PLATFORM',
    'PROPERTIES',
    'DATA_DIR',
]

__version__ = "0.1.0"
__author__ = "Kieran Russell"
__email__ = "kor20@cam.ac.uk"

