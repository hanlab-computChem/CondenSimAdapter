#!/usr/bin/env python3
"""
Legacy src/ shim.

All CG simulation logic has moved to CondenSimAdapter.core.
All backmapping has moved to CondenSimAdapter.backmap.
All minimization has moved to CondenSimAdapter.minimize.

This shim re-exports the new names under their old aliases so that
any remaining code that still imports from ``src`` continues to work.
"""

from ..backmap.backmapper import Backmapper as BackmapSimulator
from ..backmap.backmapper import BackmapResult
from ..core.config import (
    CGConfig as CGSimulationConfig,
)
from ..core.config import (
    Component as CGComponent,
)
from ..core.config import (
    ComponentType,
    SimulationResult,
    TopologyType,
)
from ..core.simulation import CGSimulation as CGSimulator
from ..minimize.config_loader import load_config_from_yaml
from ..minimize.minimizer import MinimizeConfig, MinimizeResult, MinimizeSimulator
from .pdb_tool import ChainLabel, extract_coordinates_from_pdb

# Expose plumed_generator from src (still lives here)
from .plumed_generator import generate_plumed_for_minimize as generate_plumed_dat

__all__ = [
    "CGSimulationConfig",
    "CGComponent",
    "ComponentType",
    "TopologyType",
    "SimulationResult",
    "CGSimulator",
    "BackmapSimulator",
    "BackmapResult",
    "MinimizeSimulator",
    "MinimizeConfig",
    "MinimizeResult",
    "load_config_from_yaml",
    "ChainLabel",
    "extract_coordinates_from_pdb",
    "generate_plumed_dat",
]
