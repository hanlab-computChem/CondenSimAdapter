"""
CondenSimAdapter CG Engine.

Unified OpenMM-native coarse-grained simulation engine for protein condensates.
Supports CALVADOS2/3, HPS-Urry, COCOMO2, and Mpipi-Recharged force fields.

Quick start:
    from CondenSimAdapter.core import CGConfig, CGSimulation

    config = CGConfig.from_yaml("FUS_LC.yaml")
    result = CGSimulation(config).run("FUS_LC_CG/")
"""

from .config     import CGConfig, Component, ComponentType, TopologyType, SimulationResult
from .simulation import CGSimulation
from .forcefield import create_forcefield, list_forcefields

__all__ = [
    "CGConfig",
    "Component",
    "ComponentType",
    "TopologyType",
    "SimulationResult",
    "CGSimulation",
    "create_forcefield",
    "list_forcefields",
]
