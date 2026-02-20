#!/usr/bin/env python3
"""
CLI Commands Module

Implements the command-line interface for CondenSimAdapter.

This module imports commands from the commands package for backward compatibility.
"""

# Import shared utilities for backward compatibility
from .shared import (
    CG_FORCE_FIELDS,
    FORCE_FIELD_TO_RUNNER,
    GEOMETRY_DEFAULTS,
    MINIMIZE_FORCE_FIELDS,
    validate_cg_force_field,
    validate_minimize_force_field,
    parse_component_pattern,
    parse_box,
    generate_yaml_with_comments,
)

# Import commands from the commands package
from .commands_refactored import (
    init_command,
    cg_command,
    backmap_command,
    pace_opt_command,
    minimize_command,
    info_command,
    droplet_density_command,
    to_run_command,
)

__all__ = [
    # Shared utilities
    'CG_FORCE_FIELDS',
    'FORCE_FIELD_TO_RUNNER',
    'GEOMETRY_DEFAULTS',
    'MINIMIZE_FORCE_FIELDS',
    'validate_cg_force_field',
    'validate_minimize_force_field',
    'parse_component_pattern',
    'parse_box',
    'generate_yaml_with_comments',
    # Commands
    'init_command',
    'cg_command',
    'backmap_command',
    'pace_opt_command',
    'minimize_command',
    'info_command',
    'droplet_density_command',
    'to_run_command',
]
