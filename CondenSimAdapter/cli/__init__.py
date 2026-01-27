#!/usr/bin/env python3
"""
Adapter CLI

A workflow for CG and AA protein condensate simulation.

    Core Commands:
        cg               Run coarse-grained simulation
        backmap          Backmap CG structure to all-atom representation
        minimize         Energy minimization with AMBER/CHARMM force fields

    Utility Commands:
        init             Initialize a new configuration template
        droplet-density  Estimate protein density in droplet geometry
        info             Display system and environment information
"""

import sys
import warnings
import os

# CRITICAL: Set up warning filters BEFORE any other imports
# This must be done at the very beginning to catch warnings from all modules

# Suppress all deprecation warnings globally
warnings.filterwarnings('ignore', category=DeprecationWarning)

# Suppress specific warnings from common simulation libraries
warnings.filterwarnings('ignore', message='.*simtk\\.openmm.*')
warnings.filterwarnings('ignore', message='.*xdrlib.*')
warnings.filterwarnings('ignore', message='.*MDAnalysis.*')
warnings.filterwarnings('ignore', message='.*Bio\\..*')
warnings.filterwarnings('ignore', message='.*NumPy.*')
warnings.filterwarnings('ignore', message='.*Pandas.*')

# Suppress UserWarnings and FutureWarnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Also redirect stderr for any warnings that bypass the warning module
_stderr_fileno = None
try:
    _stderr_fileno = sys.stderr.fileno()
except (AttributeError, ValueError):
    pass

# Create a filter function for stderr
class StderrWarningFilter:
    """Filter out warning messages from stderr."""

    def __init__(self):
        self.original_stderr = sys.stderr

    def write(self, text):
        # Skip warning messages
        if 'Warning:' in text or 'warning:' in text.lower():
            return
        if 'deprecated' in text.lower():
            return
        self.original_stderr.write(text)

    def flush(self):
        self.original_stderr.flush()

    def fileno(self):
        return self.original_stderr.fileno()

    def isatty(self):
        return self.original_stderr.isatty()

# Only replace stderr if we're in a non-TTY context (to avoid issues in some environments)
if _stderr_fileno is not None and not os.environ.get('FORCE_ADAPTER_STDERR', ''):
    pass  # Keep original stderr in normal use

import click

from .commands import init_command, cg_command, backmap_command, pace_opt_command, minimize_command, info_command, droplet_density_command


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    """
    Adapter: a workflow for CG and AA protein condensate simulation.
    
    CORE COMMANDS:
        cg               Run coarse-grained simulation
        backmap          Backmap CG structure to all-atom representation
        minimize         Energy minimization with AMBER/CHARMM force fields
    
    UTILITY COMMANDS:
        init             Initialize a new configuration template
        droplet-density  Estimate protein density in droplet geometry
        info             Display system and environment information
    
    Available CG force fields:
        calvados, hps_urry, cocomo, mpipi_recharged
    
    Available all-atom force fields:
        1-a99SBdisp, 2-amber03wsc, 3-amber99sbws-stqp, 4-amber99sbws-stq,
        5-des-amber, 6-des-amber-sf1.0, 7-amber99sb-ildn, 8-amber14sb,
        9-charmm36m
    
    Typical workflow:
        1. adapter init my_project              # Create configuration template
        2. adapter cg -f config.yaml            # Run CG simulation
        3. adapter backmap -i output_CG -f config.yaml  # Backmap to all-atom
        4. adapter minimize -i output_backmap -f config.yaml  # Minimize structure
    
    Additional examples:
        adapter init --topol droplet -c IIIMII  # Multi-component droplet
        adapter droplet-density -f config.yaml -r 15  # Estimate density
        adapter info                            # Check environment
    """
    pass


# Add commands (ordered: core commands first, then utility commands)
# Core commands
main.add_command(cg_command, 'cg')
main.add_command(backmap_command, 'backmap')
main.add_command(minimize_command, 'minimize')

# Utility commands
main.add_command(init_command, 'init')
main.add_command(droplet_density_command, 'droplet-density')
main.add_command(info_command, 'info')

# Hidden/experimental commands
pace_opt_command.hidden = True
main.add_command(pace_opt_command, 'pace-opt')


def cli():
    """Entry point for CLI."""
    sys.exit(main())


if __name__ == '__main__':
    cli()

