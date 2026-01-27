#!/usr/bin/env python3
"""
CLI Entry Point

Allows running the CLI with: python -m CondenSimAdapter.cli
"""

import sys
import os
import warnings


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


# Replace stderr with filtered version BEFORE any imports
sys.stderr = FilteredStderr(sys.stderr)

# Set up comprehensive warning filters BEFORE any other imports
warnings.filterwarnings('ignore')
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*')
warnings.filterwarnings('ignore', category=Warning)

from . import cli

if __name__ == '__main__':
    cli()

