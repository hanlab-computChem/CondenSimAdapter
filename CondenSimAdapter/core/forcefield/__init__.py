"""
CG force field registry.

Usage:
    from CondenSimAdapter.core.forcefield import create_forcefield
    ff = create_forcefield("calvados2")
"""

from __future__ import annotations

from typing import Dict, Type

from .base import CGForceField
from .calvados import CalvadosFF
from .hps     import HPSFF
from .cocomo  import CocomoFF
from .mpipi   import MpipiFF

_REGISTRY: Dict[str, Type[CGForceField]] = {
    "calvados2": CalvadosFF,
    "calvados3": CalvadosFF,
    "hps"      : HPSFF,
    "cocomo"   : CocomoFF,
    "mpipi"    : MpipiFF,
}


def create_forcefield(name: str) -> CGForceField:
    """
    Instantiate a CG force field by name.

    Supported names: calvados2, calvados3, hps, cocomo, mpipi.
    """
    key = name.lower().strip()
    if key not in _REGISTRY:
        raise ValueError(
            f"Unknown CG force field '{name}'. "
            f"Available: {sorted(_REGISTRY)}"
        )
    if key == "calvados3":
        return CalvadosFF(version=3)
    if key == "calvados2":
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
