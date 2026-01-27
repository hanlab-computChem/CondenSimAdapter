#!/usr/bin/env python3
"""
CondenSimAdapter Package

A multi-stage workflow for simulating protein condensates.

Subpackages:
- src: Coarse-grained simulation framework
- extern: External packages (ms2_calvados, etc.)
"""

import sys
import warnings

# Store original stderr
_original_stderr = sys.stderr if hasattr(sys, 'stderr') else None


class FilteredStderr:
    """Custom stderr wrapper that filters out all warning messages."""

    def __init__(self, original_stderr):
        self.original_stderr = original_stderr

    def write(self, text):
        # Filter out all warning messages (case-insensitive)
        if 'warning' in text.lower():
            return
        self.original_stderr.write(text)

    def flush(self):
        self.original_stderr.flush()

    def fileno(self):
        return self.original_stderr.fileno()

    def isatty(self):
        return self.original_stderr.isatty()

    def __getattr__(self, name):
        return getattr(self.original_stderr, name)


# Replace stderr with filtered version BEFORE any imports
if _original_stderr is not None:
    sys.stderr = FilteredStderr(_original_stderr)

# Set up comprehensive warning filters
warnings.filterwarnings('ignore')
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*')

from .src import (
    CGSimulationConfig,
    CGComponent,
    ComponentType,
    TopologyType,
    ComputePlatform,
    SimulationParams,
    SimulationResult,
    CGSimulator,
    CalvadosWrapper,
    run_calvados,
)

__all__ = [
    # Configuration
    'CGSimulationConfig',
    'CGComponent',
    'ComponentType',
    'TopologyType',
    'ComputePlatform',
    'SimulationParams',
    'SimulationResult',

    # Simulator
    'CGSimulator',

    # CALVADOS wrapper
    'CalvadosWrapper',
    'run_calvados',
]

