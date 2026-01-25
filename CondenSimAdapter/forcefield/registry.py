#!/usr/bin/env python3
"""
Force Field Registry

Central registry for all available force fields with their metadata including
pdb2gmx names, water models, and GBSA atom type mappings.

This registry supports numbered selection:
1. a99SBdisp    4. amber99sbws-stq  7. amber99sb-ildn
2. amber03wsc   5. des-amber        8. amber14sb
3. amber99sbws-stqp  6. des-amber-sf1.0  9. charmm36m
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Dict


# =============================================================================
# Force Field Info Data Class
# =============================================================================

@dataclass
class ForceFieldInfo:
    """Information about a force field.
    
    Attributes:
        name: CLI name for the force field (short identifier with number)
        family: Force field family (AMBER or CHARMM)
        pdb2gmx_name: Name used by GROMACS pdb2gmx
        water_model: Water model for pdb2gmx
        solvate_cs: Water coordinate model for gmx solvate -cs
        gbsa_mapping: GBSA atom type mapping to use
        description: Human-readable description
    """
    name: str
    family: str
    pdb2gmx_name: str
    water_model: str
    solvate_cs: str
    gbsa_mapping: str
    description: str = ""


# =============================================================================
# Built-in Force Field Definitions
# =============================================================================

# Force fields ordered by number for CLI selection
# Format: "N-name" where N is the selection number
BUILTIN_FORCE_FIELDS = [
    # AMBER family (numbered 1-8)
    ForceFieldInfo(
        name="1-a99SBdisp",
        family="AMBER",
        pdb2gmx_name="a99SBdisp",
        water_model="a99SBdisp_water",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="a99SB-disp with custom water model"
    ),
    ForceFieldInfo(
        name="2-amber03wsc",
        family="AMBER",
        pdb2gmx_name="amber03wsc",
        water_model="tip4p2005s",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER03 with WSC correction and tip4p2005s water"
    ),
    ForceFieldInfo(
        name="3-amber99sbws-stqp",
        family="AMBER",
        pdb2gmx_name="amber99sbws-STQp",
        water_model="tip4p2005s",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER99SB-WS with STQp correction and tip4p2005s water"
    ),
    ForceFieldInfo(
        name="4-amber99sbws-stq",
        family="AMBER",
        pdb2gmx_name="amber99sbws-stq",
        water_model="tip4p2005s",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER99SB-WS with stq correction and tip4p2005s water"
    ),
    ForceFieldInfo(
        name="5-des-amber",
        family="AMBER",
        pdb2gmx_name="des-amber",
        water_model="tip4pd",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="DES-AMBER with tip4pd water"
    ),
    ForceFieldInfo(
        name="6-des-amber-sf1.0",
        family="AMBER",
        pdb2gmx_name="des-amber-SF1.0",
        water_model="tip4pd",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="DES-AMBER SF1.0 with tip4pd water"
    ),
    ForceFieldInfo(
        name="7-amber99sb-ildn",
        family="AMBER",
        pdb2gmx_name="amber99sb-ildn",
        water_model="tip3p",
        solvate_cs="spc216",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER99SB-ILDN with tip3p water (default)"
    ),
    ForceFieldInfo(
        name="8-amber14sb",
        family="AMBER",
        pdb2gmx_name="amber14sb_parmbsc1",
        water_model="tip3p",
        solvate_cs="spc216",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER14SB with PARMBSC1 correction and tip3p water"
    ),
    
    # CHARMM family (numbered 9)
    ForceFieldInfo(
        name="9-charmm36m",
        family="CHARMM",
        pdb2gmx_name="charmm36-jul2021",
        water_model="tip3p",
        solvate_cs="spc216",
        gbsa_mapping="CHARMM36",
        description="CHARMM36m-jul2021 with tip3p water"
    ),
]


# =============================================================================
# Force Field Registry
# =============================================================================

class ForceFieldRegistry:
    """Registry for managing available force fields.
    
    This class provides a centralized way to:
    - List all available force fields
    - Get force field information by name (supports both "N-name" and "pdb2gmx_name")
    - Get water model for a force field
    - Validate force field names
    - Filter force fields by family
    
    Examples:
        >>> registry = ForceFieldRegistry()
        >>> registry.list_force_fields()
        ['1-a99SBdisp', '2-amber03wsc', '3-amber99sbws-stqp', ...]
        >>> ff = registry.get_force_field('1-a99SBdisp')
        >>> ff.water_model
        'a99SBdisp_water'
        >>> ff = registry.get_force_field('a99SBdisp')
        >>> ff.water_model
        'a99SBdisp_water'
    """
    
    def __init__(self):
        """Initialize the registry with built-in force fields."""
        self._force_fields: Dict[str, ForceFieldInfo] = {}
        self._pdb2gmx_index: Dict[str, str] = {}  # Map pdb2gmx_name to CLI name
        
        for ff in BUILTIN_FORCE_FIELDS:
            # Store CLI name (e.g., "1-a99SBdisp") in lowercase for case-insensitive lookup
            self._force_fields[ff.name.lower()] = ff
            # Also map by pdb2gmx_name for convenience
            self._pdb2gmx_index[ff.pdb2gmx_name.lower()] = ff.name.lower()
            # Add charmm36m as an alias for charmm36-jul2021
            if ff.pdb2gmx_name == "charmm36-jul2021":
                self._pdb2gmx_index["charmm36m"] = ff.name.lower()
    
    def list_force_fields(self) -> List[str]:
        """List all available force field names.
        
        Returns:
            List of force field names in order
        """
        return [ff.name for ff in BUILTIN_FORCE_FIELDS]
    
    def list_by_family(self, family: str) -> List[str]:
        """List force field names by family.
        
        Args:
            family: Family name (AMBER or CHARMM)
        
        Returns:
            List of force field names in the family (in order)
        """
        family = family.upper()
        return [
            ff.name for ff in BUILTIN_FORCE_FIELDS
            if ff.family.upper() == family
        ]
    
    def get_force_field(self, name: str) -> Optional[ForceFieldInfo]:
        """Get force field information by name.
        
        Supports:
        - CLI name: "1-a99SBdisp", "2-amber03wsc", etc.
        - Short number: "1", "2", etc.
        - pdb2gmx name: "a99SBdisp", "amber03wsc", "charmm36m", etc.
        - Alias: "charmm36m" maps to CHARMM36
        
        Args:
            name: Force field name (case-insensitive)
        
        Returns:
            ForceFieldInfo if found, None otherwise
        """
        name_lower = name.lower()
        
        # Handle short number format (e.g., "1" -> "1-a99SBdisp")
        if name_lower.isdigit():
            for ff in BUILTIN_FORCE_FIELDS:
                if ff.name.startswith(f"{name_lower}-"):
                    return ff
            return None
        
        # Try direct lookup first (for "1-a99SBdisp")
        if name_lower in self._force_fields:
            return self._force_fields[name_lower]
        
        # Try pdb2gmx name lookup (for "a99SBdisp", "charmm36m", etc.)
        if name_lower in self._pdb2gmx_index:
            cli_name = self._pdb2gmx_index[name_lower]
            return self._force_fields[cli_name]
        
        return None
    
    def get_pdb2gmx_name(self, name: str) -> Optional[str]:
        """Get pdb2gmx name for a force field.
        
        Args:
            name: Force field name (CLI name or pdb2gmx name)
        
        Returns:
            pdb2gmx name if found, None otherwise
        """
        ff = self.get_force_field(name)
        return ff.pdb2gmx_name if ff else None
    
    def get_water_model(self, name: str) -> Optional[str]:
        """Get water model for a force field.
        
        Args:
            name: Force field name
        
        Returns:
            Water model name if found, None otherwise
        """
        ff = self.get_force_field(name)
        return ff.water_model if ff else None

    def get_solvate_cs(self, name: str) -> Optional[str]:
        """Get gmx solvate -cs water model for a force field.
        
        Args:
            name: Force field name
        
        Returns:
            gmx solvate -cs model name if found, None otherwise
        """
        ff = self.get_force_field(name)
        return ff.solvate_cs if ff else None
    
    def get_gbsa_mapping(self, name: str) -> Optional[str]:
        """Get GBSA atom type mapping for a force field.
        
        Args:
            name: Force field name
        
        Returns:
            GBSA mapping name if found, None otherwise
        """
        ff = self.get_force_field(name)
        return ff.gbsa_mapping if ff else None
    
    def get_family(self, name: str) -> Optional[str]:
        """Get family for a force field.

        Args:
            name: Force field name

        Returns:
            Family name (AMBER or CHARMM) if found, None otherwise
        """
        ff = self.get_force_field(name)
        return ff.family if ff else None

    def get_force_field_path(self, name: str) -> Optional[Path]:
        """Get the path to the force field folder.

        This returns the path to the force field folder in the package
        (e.g., CondenSimAdapter/forcefield/amber99sb-ildn.ff).

        Args:
            name: Force field name

        Returns:
            Path to force field folder if found, None otherwise
        """
        import importlib.resources as resources
        from .. import forcefield

        ff = self.get_force_field(name)
        if not ff:
            return None

        # Force field folder name matches pdb2gmx_name (e.g., "amber99sb-ildn" -> "amber99sb-ildn.ff")
        ff_folder_name = f"{ff.pdb2gmx_name}.ff"

        # Try to get the path using importlib.resources
        try:
            # For Python 3.9+
            ff_path = resources.files(forcefield) / ff_folder_name
            # Convert to Path if possible
            if hasattr(ff_path, '__fspath__'):
                return Path(str(ff_path))
            else:
                # It's a Traversable, try to resolve
                for p in resources.files(forcefield).iterdir():
                    if p.name == ff_folder_name:
                        return Path(str(p))
        except Exception:
            pass

        # Fallback: construct path directly
        package_dir = Path(__file__).parent.parent
        ff_path = package_dir / "forcefield" / ff_folder_name
        return ff_path if ff_path.exists() else None

    def validate(self, name: str) -> tuple[bool, str]:
        """Validate a force field name.
        
        Args:
            name: Force field name to validate
        
        Returns:
            Tuple of (is_valid, message)
        """
        ff = self.get_force_field(name)
        if ff:
            return True, f"Valid force field: {ff.name} ({ff.family})"
        return False, f"Unknown force field: {name}. Available: {', '.join(self.list_force_fields())}"
    
    def is_amber(self, name: str) -> bool:
        """Check if a force field is from AMBER family.
        
        Args:
            name: Force field name
        
        Returns:
            True if AMBER family, False otherwise
        """
        return self.get_family(name) == "AMBER"
    
    def is_charmm(self, name: str) -> bool:
        """Check if a force field is from CHARMM family.
        
        Args:
            name: Force field name
        
        Returns:
            True if CHARMM family, False otherwise
        """
        return self.get_family(name) == "CHARMM"
    
    def count(self) -> int:
        """Get total number of registered force fields.
        
        Returns:
            Number of force fields
        """
        return len(self._force_fields)


# =============================================================================
# Global Registry Instance
# =============================================================================

# Global registry instance for easy access
REGISTRY = ForceFieldRegistry()


# =============================================================================
# Convenience Functions
# =============================================================================

def list_force_fields() -> List[str]:
    """List all available force field names.
    
    Returns:
        List of force field names in order
    """
    return REGISTRY.list_force_fields()


def list_amber_force_fields() -> List[str]:
    """List available AMBER force field names.
    
    Returns:
        List of AMBER force field names (in order)
    """
    return REGISTRY.list_by_family("AMBER")


def list_charmm_force_fields() -> List[str]:
    """List available CHARMM force field names.
    
    Returns:
        List of CHARMM force field names (in order)
    """
    return REGISTRY.list_by_family("CHARMM")


def get_force_field(name: str) -> Optional[ForceFieldInfo]:
    """Get force field information by name.
    
    Args:
        name: Force field name (CLI name or pdb2gmx name)
    
    Returns:
        ForceFieldInfo if found, None otherwise
    """
    return REGISTRY.get_force_field(name)


def get_water_model(name: str) -> Optional[str]:
    """Get water model for a force field.

    Args:
        name: Force field name

    Returns:
        Water model name if found, None otherwise
    """
    return REGISTRY.get_water_model(name)


def get_force_field_path(name: str) -> Optional[str]:
    """Get the path to the force field folder.

    This returns the path to the force field folder in the package
    (e.g., CondenSimAdapter/forcefield/amber99sb-ildn.ff).

    Args:
        name: Force field name

    Returns:
        Path string to force field folder if found, None otherwise
    """
    path = REGISTRY.get_force_field_path(name)
    return str(path) if path else None


def validate_force_field(name: str) -> tuple[bool, str]:
    """Validate a force field name.
    
    Args:
        name: Force field name to validate
    
    Returns:
        Tuple of (is_valid, message)
    """
    return REGISTRY.validate(name)


if __name__ == "__main__":
    # Print all available force fields
    print("Available force fields:")
    print("-" * 60)
    for ff in BUILTIN_FORCE_FIELDS:
        print(f"{ff.name:22} | {ff.family:6} | {ff.water_model:15} | {ff.gbsa_mapping}")
    print("-" * 60)
    print(f"Total: {len(BUILTIN_FORCE_FIELDS)} force fields")
