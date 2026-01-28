"""
ms2_calvados - Internalized CALVADOS coarse-grained simulation package
Retain core functions: system construction + OpenMM simulation
"""

# Core modules
from .calvados import cfg
from .calvados import components
from .calvados import build
from .calvados import sim
from .calvados import sequence
from .calvados import interactions

# Expose key classes
from .calvados.cfg import Config, Job, Components
from .calvados.components import Component, Protein, RNA, Lipid, Crowder

__all__ = [
    'cfg',
    'components',
    'build',
    'sim',
    'sequence',
    'interactions',
    'Config',
    'Job',
    'Components',
    'Component',
    'Protein',
    'RNA',
    'Lipid',
    'Crowder',
]

