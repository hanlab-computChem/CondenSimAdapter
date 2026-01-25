#!/usr/bin/env python3
"""
CLI Entry Point

Allows running the CLI with: python -m CondenSimAdapter.cli
"""

import sys
import os
import warnings


class FilteredStderr:
    """Custom stderr wrapper that filters out OpenMM deprecation warnings."""
    
    def __init__(self, original_stderr):
        self.original_stderr = original_stderr
        self.buffer = []
        
    def write(self, text):
        # Filter out OpenMM deprecation warnings
        if 'simtk.openmm' in text.lower() and 'deprecated' in text.lower():
            return
        self.original_stderr.write(text)
        
    def flush(self):
        self.original_stderr.flush()
        
    def fileno(self):
        return self.original_stderr.fileno()


# Replace stderr with filtered version
sys.stderr = FilteredStderr(sys.stderr)

# Suppress deprecation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

from . import cli

if __name__ == '__main__':
    cli()

