#!/usr/bin/env python3
"""
CondenSimAdapter

Automated multi-stage workflow for protein condensate simulations:
  CG simulation -> CG-to-AA backmapping -> softcore minimization -> production

Architecture:
  core/      - Unified OpenMM CG engine (CALVADOS, HPS, COCOMO, Mpipi)
  backmap/   - cg2all neural-network backmapper
  minimize/  - three-stage softcore AA minimization
  cli/       - adapter command-line interface
  forcefield/- AA force field registry (.ff directories)
"""

from .core import (
    CGConfig,
    Component,
    ComponentType,
    TopologyType,
    SimulationResult,
    CGSimulation,
    create_forcefield,
    list_forcefields,
)

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
