"""
CG force field registry.

Usage:
    from CondenSimAdapter.core.forcefield import create_forcefield
    ff = create_forcefield("calvados")   # version auto-selected by CGConfig
    ff = create_forcefield("calvados2")  # explicit
"""

from __future__ import annotations

from typing import Dict, Type

from .base import CGForceField
from .calvados import CalvadosFF
from .cocomo import CocomoFF
from .hps import HPSFF
from .mpipi import MpipiFF

_REGISTRY: Dict[str, Type[CGForceField]] = {
    "calvados": CalvadosFF,  # resolved to 2 or 3 by CGConfig.resolved_force_field
    "calvados2": CalvadosFF,
    "calvados3": CalvadosFF,
    "hps": HPSFF,
    "cocomo": CocomoFF,
    "mpipi": MpipiFF,
}


def create_forcefield(name: str) -> CGForceField:
    """
    Instantiate a CG force field by name.

    Canonical names: calvados, hps, cocomo, mpipi.
    Explicit version aliases: calvados2, calvados3.
    Pass CGConfig.resolved_force_field to always get the versioned instance.
    """
    key = name.lower().strip()
    if key not in _REGISTRY:
        raise ValueError(f"Unknown CG force field '{name}'. Available: {sorted(_REGISTRY)}")
    if key == "calvados3":
        return CalvadosFF(version=3)
    # both 'calvados' and 'calvados2' default to version 2
    if key in ("calvados", "calvados2"):
        return CalvadosFF(version=2)
    return _REGISTRY[key]()


def list_forcefields():
    return sorted(_REGISTRY.keys())


__all__ = [
    "CGForceField",
    "CalvadosFF",
    "HPSFF",
    "CocomoFF",
    "MpipiFF",
    "create_forcefield",
    "list_forcefields",
]
