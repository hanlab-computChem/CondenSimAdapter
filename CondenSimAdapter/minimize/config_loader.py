"""
YAML config loader shared by CLI and MinimizeSimulator.

Replaces the load_config_from_yaml / load_components_from_yaml functions
that were previously buried inside pdb2gmx_utils.py.
"""

from __future__ import annotations

import yaml
from pathlib import Path
from typing import List, Optional, Tuple


def load_config_from_yaml(yaml_path: str) -> Tuple[str, list]:
    """
    Load a simulation config YAML and return (system_name, components).

    The components list contains raw dicts suitable for passing to the
    Component.from_dict() constructor.

    Args:
        yaml_path: Path to the YAML configuration file.

    Returns:
        (system_name, components_list)
    """
    path = Path(yaml_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {yaml_path}")

    with open(path) as f:
        raw = yaml.safe_load(f)

    if raw is None:
        raise ValueError(f"Empty or invalid YAML: {yaml_path}")

    # Support both 'system_name' and legacy 'sysname'
    system_name = raw.get("system_name") or raw.get("sysname") or path.stem

    components = raw.get("components", [])
    return system_name, components


def load_components_from_yaml(yaml_path: str) -> list:
    """Return only the components list from a config YAML."""
    _, components = load_config_from_yaml(yaml_path)
    return components


def get_system_name(yaml_path: str) -> str:
    """Return just the system_name from a config YAML."""
    name, _ = load_config_from_yaml(yaml_path)
    return name
