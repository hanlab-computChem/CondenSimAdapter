#!/usr/bin/env python3
"""
Shared utilities for CLI commands.

Contains common constants, validation functions, and helper utilities
used across multiple CLI commands.
"""

import os
import sys
import warnings

# CRITICAL: Set up warning filters BEFORE any other imports
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

from pathlib import Path
from typing import Optional, List

import click

from ..src import CGSimulationConfig, CGComponent, ComponentType, TopologyType
from ..src.minimize import MinimizeSimulator, MinimizeConfig
from ..forcefield.registry import REGISTRY, list_force_fields, BUILTIN_FORCE_FIELDS


# Available force fields (for CG simulation)
CG_FORCE_FIELDS = ['calvados', 'hps_urry', 'cocomo', 'mpipi_recharged']

# Mapping from CLI force field names to internal runner method names
FORCE_FIELD_TO_RUNNER = {
    'calvados': 'calvados',
    'hps_urry': 'hps',  # CLI uses hps_urry, but method is run_hps
    'cocomo': 'cocomo',
    'mpipi_recharged': 'mpipi_recharged',
}

# Geometry defaults for init command
GEOMETRY_DEFAULTS = {
    'grid': {
        'box': [20.0, 20.0, 20.0],
        'description': 'Continuous dense phase with periodic boundaries in x, y, z',
    },
    'slab': {
        'box': [10.0, 10.0, 40.0],
        'description': 'Slab geometry with periodic boundaries in x, y and interfaces with dilute phase along z',
    },
    'droplet': {
        'box': [15.0, 15.0, 15.0],  # default radius
        'description': 'Spherical droplet confined within radius r, surrounded by dilute phase',
    },
}


def validate_cg_force_field(ctx, param, value):
    """Validate CG force field name."""
    if value and value not in CG_FORCE_FIELDS:
        raise click.BadParameter(
            f"Invalid force field: {value}. "
            f"Available: {', '.join(CG_FORCE_FIELDS)}"
        )
    return value.lower() if value else None


# For minimize command: available force fields from registry
MINIMIZE_FORCE_FIELDS = list_force_fields()


def validate_minimize_force_field(ctx, param, value):
    """Validate minimize force field name using registry.
    
    Supports:
    - Short number: "1", "2", etc. (converts to CLI name format)
    - CLI name: "1-a99SBdisp", "2-amber03wsc", etc.
    - pdb2gmx name: "a99SBdisp", "amber03wsc", "charmm36m", etc.
    """
    if value:
        # If user enters just a number, convert to CLI name format
        if value.isdigit():
            # Find the corresponding CLI name
            for ff in BUILTIN_FORCE_FIELDS:
                if ff.name.startswith(f"{value}-"):
                    return ff.name
            # Number not found, will fail validation below
        
        is_valid, message = REGISTRY.validate(value)
        if not is_valid:
            raise click.BadParameter(message)
        
        # Return the canonical CLI name (e.g., "1-a99SBdisp")
        ff = REGISTRY.get_force_field(value)
        return ff.name if ff else value
    return '7-amber99sb-ildn'  # Default: amber99sb-ildn


def parse_component_pattern(pattern: str, nmol: int) -> List[dict]:
    """Parse component pattern string like 'IIIMII' -> 3 IDP + 1 MDP + 2 IDP"""
    components = []
    idp_count = 0
    mdp_count = 0
    
    for char in pattern.upper():
        if char == 'I':
            idp_count += 1
            comp = {
                'name': f'IDP_{idp_count}',
                'type': 'IDP',
                'nmol': nmol,
                'ffasta': f'input/IDP_{idp_count}.fasta',
            }
            components.append(comp)
        elif char == 'M':
            mdp_count += 1
            comp = {
                'name': f'MDP_{mdp_count}',
                'type': 'MDP',
                'nmol': nmol,
                'fpdb': f'input/MDP_{mdp_count}.pdb',
                'restraint': True,
                'restraint_type': 'harmonic',
                'charge_termini': 'both',
            }
            components.append(comp)
        else:
            raise ValueError(f"Invalid character '{char}' in component pattern. Use 'I' for IDP, 'M' for MDP.")
    
    if not components:
        raise ValueError("Component pattern cannot be empty. Use 'I' for IDP, 'M' for MDP.")
    
    return components


def parse_box(box_tuple: Optional[tuple], topol: str, default_box: List[float]) -> List[float]:
    """Parse box parameter (now accepts tuple from nargs=3)"""
    if box_tuple is None:
        return default_box
    
    # Convert tuple to list
    values = list(box_tuple)
    
    if topol == 'droplet':
        # Droplet: use x value as radius for all dimensions
        return [values[0], values[0], values[0]]
    else:
        # grid/slab: use as-is
        return values


def generate_yaml_with_comments(config: 'CGSimulationConfig', topol: str, 
                                  geom_description: str, time_ns: float, 
                                  component_list: List[dict], force_field: str = 'calvados') -> str:
    """Generate YAML with detailed comments"""
    import yaml
    
    # Header comments
    lines = [
        "# CG Simulation Configuration",
        "# Generated by: adapter init",
        "#",
        "# Available CG force fields:",
        "#   - calvados: CALVADOS force field (default)",
        "#   - hps_urry: HPS Urry model",
        "#   - cocomo: COCOMO2",
        "#   - mpipi_recharged: Mpipi_recharged",
        "#",
        "# Geometry types:",
        "#   - grid: Continuous dense phase with periodic boundaries in x, y, z",
        "#   - slab: Periodic boundaries in x, y and interfaces with dilute phase along z",
        "#   - droplet: Spherical droplet confined within radius r, surrounded by dilute phase",
        "",
        f"# System information",
        f"system_name: {config.system_name}",
        "",
        "# CG force field",
        f"force_field: {force_field}   # calvados | hps_urry | cocomo | mpipi_recharged",
        "",
        "# Environment parameters",
        f"box: [{', '.join(f'{v:.1f}' for v in config.box)}]   # nm (x, y, z)",
        f"temperature: {config.temperature}         # Kelvin",
        f"ionic: {config.ionic}                # Molar (ionic strength)",
        "",
        "# Topology type:",
        "#   - grid: Continuous dense phase with periodic boundaries in x, y, z.",
        "#   - droplet: Spherical droplet confined within radius r, surrounded by dilute phase",
        "#   - slab: Geometry with periodic boundaries in x, y and interfaces with dilute phase along z.",
        f"topol: {topol}",
        "",
        "# CG simulation parameters",
        "simulation:",
        f"  steps: {config.simulation.steps}   # {time_ns} ns (1 step = 10 fs)",
        f"  wfreq: {config.simulation.wfreq}        # write frequency - save per 50 ps",
        f"  verbose: {str(config.simulation.verbose).lower()}",
        "",
        "# Component definitions",
        "components:",
    ]
    
    # Add components with comments
    for i, comp in enumerate(component_list):
        if i > 0:
            lines.append("")
        
        comp_type = comp['type']
        lines.append(f"  - name: {comp['name']}")
        lines.append(f"    type: {comp_type}          # IDP or MDP")
        lines.append(f"    nmol: {comp['nmol']}           # number of molecules (can be adjusted per component)")
        
        if comp_type == 'IDP':
            lines.append(f"    ffasta: {comp['ffasta']}")
        elif comp_type == 'MDP':
            lines.append(f"    fpdb: {comp['fpdb']}")
            lines.append("    # Domain definitions (required for MDP with restraints):")
            lines.append("    fdomains: |")
            lines.append(f"      {comp['name']}:")
            lines.append("        - [1, 50]    # Domain 1: residues 1-50")
            lines.append("        - [51, 100]  # Domain 2: residues 51-100")
            lines.append("    # Alternative options:")
            lines.append("    # Option 1 - Inline YAML (as above, recommended)")
            lines.append("    # Option 2 - File path:")
            lines.append(f"    #   fdomains: {comp['name']}_domains.yaml")
            if comp.get('restraint'):
                lines.append(f"    restraint: {str(comp['restraint']).lower()}")
                lines.append(f"    restraint_type: {comp.get('restraint_type', 'harmonic')}  # harmonic | go")
            if comp.get('charge_termini'):
                lines.append(f"    charge_termini: {comp['charge_termini']}  # both | n | c | none")
    
    return '\n'.join(lines) + '\n'
