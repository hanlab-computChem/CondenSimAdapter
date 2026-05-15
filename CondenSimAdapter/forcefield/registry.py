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

import json
import logging
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

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
    source: str = "builtin"  # "builtin" or "custom"
    ff_dir: Optional[str] = None  # Relative or absolute path to *.ff folder


CUSTOM_ID_PATTERN = re.compile(r"^a([1-9][0-9]*)$")

# Old source-tree location (for migration)
_OLD_CUSTOM_FORCEFIELD_DIR = Path(__file__).parent / "custom"
_OLD_CUSTOM_FORCEFIELD_INDEX = Path(__file__).parent / "user_forcefields.json"


def _get_user_data_dir() -> Path:
    """Resolve user data directory for custom force fields.

    Priority: CONDENSIMADAPTER_DATA_DIR env var > platformdirs > ~/.config/CondenSimAdapter
    """
    env_dir = os.environ.get("CONDENSIMADAPTER_DATA_DIR")
    if env_dir:
        return Path(env_dir)
    try:
        from platformdirs import user_data_dir

        return Path(user_data_dir("CondenSimAdapter"))
    except ImportError:
        return Path.home() / ".config" / "CondenSimAdapter"


_USER_DATA_DIR = _get_user_data_dir()
CUSTOM_FORCEFIELD_DIR = _USER_DATA_DIR / "forcefields"
CUSTOM_FORCEFIELD_INDEX = _USER_DATA_DIR / "user_forcefields.json"


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
        description="a99SB-disp with custom water model",
    ),
    ForceFieldInfo(
        name="2-amber03wsc",
        family="AMBER",
        pdb2gmx_name="amber03wsc",
        water_model="tip4p2005s",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER03 with WSC correction and tip4p2005s water",
    ),
    ForceFieldInfo(
        name="3-amber99sbws-stqp",
        family="AMBER",
        pdb2gmx_name="amber99sbws-STQp",
        water_model="tip4p2005s",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER99SB-WS with STQp correction and tip4p2005s water",
    ),
    ForceFieldInfo(
        name="4-amber99sbws-stq",
        family="AMBER",
        pdb2gmx_name="amber99sbws-stq",
        water_model="tip4p2005s",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER99SB-WS with stq correction and tip4p2005s water",
    ),
    ForceFieldInfo(
        name="5-des-amber",
        family="AMBER",
        pdb2gmx_name="des-amber",
        water_model="tip4pd",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="DES-AMBER with tip4pd water",
    ),
    ForceFieldInfo(
        name="6-des-amber-sf1.0",
        family="AMBER",
        pdb2gmx_name="des-amber-SF1.0",
        water_model="tip4pd",
        solvate_cs="tip4p",
        gbsa_mapping="AMBER99SB-ILDN",
        description="DES-AMBER SF1.0 with tip4pd water",
    ),
    ForceFieldInfo(
        name="7-amber99sb-ildn",
        family="AMBER",
        pdb2gmx_name="amber99sb-ildn",
        water_model="tip3p",
        solvate_cs="spc216",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER99SB-ILDN with tip3p water (default)",
    ),
    ForceFieldInfo(
        name="8-amber14sb",
        family="AMBER",
        pdb2gmx_name="amber14sb_parmbsc1",
        water_model="tip3p",
        solvate_cs="spc216",
        gbsa_mapping="AMBER99SB-ILDN",
        description="AMBER14SB with PARMBSC1 correction and tip3p water",
    ),
    # CHARMM family (numbered 9)
    ForceFieldInfo(
        name="9-charmm36m",
        family="CHARMM",
        pdb2gmx_name="charmm36-jul2021",
        water_model="tip3p",
        solvate_cs="spc216",
        gbsa_mapping="CHARMM36",
        description="CHARMM36m-jul2021 with tip3p water",
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
        self._custom_force_fields: Dict[str, ForceFieldInfo] = {}

        self._migrate_old_custom_forcefields()

        for ff in BUILTIN_FORCE_FIELDS:
            ff.source = "builtin"
            # Store CLI name (e.g., "1-a99SBdisp") in lowercase for case-insensitive lookup
            self._force_fields[ff.name.lower()] = ff
            # Also map by pdb2gmx_name for convenience
            self._pdb2gmx_index[ff.pdb2gmx_name.lower()] = ff.name.lower()
            # Add charmm36m as an alias for charmm36-jul2021
            if ff.pdb2gmx_name == "charmm36-jul2021":
                self._pdb2gmx_index["charmm36m"] = ff.name.lower()

        self._load_custom_force_fields()

    def _migrate_old_custom_forcefields(self) -> None:
        """Migrate custom force fields from old source-tree location to user data dir."""
        if CUSTOM_FORCEFIELD_INDEX.exists():
            return  # Already migrated or fresh install
        if not _OLD_CUSTOM_FORCEFIELD_INDEX.exists():
            return  # Nothing to migrate
        try:
            CUSTOM_FORCEFIELD_DIR.mkdir(parents=True, exist_ok=True)
            shutil.copy2(_OLD_CUSTOM_FORCEFIELD_INDEX, CUSTOM_FORCEFIELD_INDEX)
            if _OLD_CUSTOM_FORCEFIELD_DIR.exists():
                for item in _OLD_CUSTOM_FORCEFIELD_DIR.iterdir():
                    dest = CUSTOM_FORCEFIELD_DIR / item.name
                    if not dest.exists():
                        if item.is_dir():
                            shutil.copytree(item, dest)
                        else:
                            shutil.copy2(item, dest)
        except OSError:
            pass  # Migration is best-effort

    def _custom_id_key(self, ff_name: str) -> int:
        """Return numeric part for custom id sorting (a1, a2, ...)."""
        match = CUSTOM_ID_PATTERN.match(ff_name.lower())
        return int(match.group(1)) if match else 0

    def _to_index_data(self, ff: ForceFieldInfo) -> Dict[str, Any]:
        """Serialize custom force field for JSON index."""
        return {
            "id": ff.name,
            "family": ff.family,
            "pdb2gmx_name": ff.pdb2gmx_name,
            "water_model": ff.water_model,
            "solvate_cs": ff.solvate_cs,
            "gbsa_mapping": ff.gbsa_mapping,
            "description": ff.description,
            "ff_dir": ff.ff_dir,
        }

    def _from_index_data(self, data: Dict[str, Any]) -> Optional[ForceFieldInfo]:
        """Deserialize one custom force field record from JSON index."""
        ff_id = str(data.get("id", "")).strip().lower()
        if not CUSTOM_ID_PATTERN.match(ff_id):
            return None

        pdb2gmx_name = str(data.get("pdb2gmx_name", "")).strip()
        if not pdb2gmx_name:
            return None

        ff_dir = data.get("ff_dir")
        if ff_dir is None:
            return None

        ff = ForceFieldInfo(
            name=ff_id,
            family=str(data.get("family", "AMBER")).upper(),
            pdb2gmx_name=pdb2gmx_name,
            water_model=str(data.get("water_model", "tip3p")),
            solvate_cs=str(data.get("solvate_cs", "spc216")),
            gbsa_mapping=str(data.get("gbsa_mapping", "AMBER99SB-ILDN")),
            description=str(data.get("description", "")),
            source="custom",
            ff_dir=str(ff_dir),
        )
        return ff

    def _save_custom_force_fields(self) -> None:
        """Persist custom force fields to source-tree JSON index.

        Uses atomic write (temp file + rename) to avoid corrupting the index on
        partial writes.
        """
        CUSTOM_FORCEFIELD_INDEX.parent.mkdir(parents=True, exist_ok=True)
        records = [
            self._to_index_data(ff)
            for ff in sorted(
                self._custom_force_fields.values(), key=lambda item: self._custom_id_key(item.name)
            )
        ]
        fd, tmp_path = tempfile.mkstemp(
            suffix=".json", prefix=".tmp_forcefield_", dir=CUSTOM_FORCEFIELD_INDEX.parent
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                json.dump({"force_fields": records}, f, indent=2, sort_keys=True)
                f.write("\n")
            os.replace(tmp_path, CUSTOM_FORCEFIELD_INDEX)
        except Exception:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
            raise

    def _load_custom_force_fields(self) -> None:
        """Load custom force fields from source-tree JSON index."""
        if not CUSTOM_FORCEFIELD_INDEX.exists():
            return

        try:
            with open(CUSTOM_FORCEFIELD_INDEX, "r", encoding="utf-8") as f:
                payload = json.load(f)
        except Exception:
            logger.warning(
                "Failed to load custom force field index from %s (file may be corrupted)",
                CUSTOM_FORCEFIELD_INDEX,
                exc_info=True,
            )
            return

        records = payload.get("force_fields", []) if isinstance(payload, dict) else []
        if not isinstance(records, list):
            return

        for record in records:
            if not isinstance(record, dict):
                continue
            ff = self._from_index_data(record)
            if ff is None:
                continue
            if ff.name in self._force_fields:
                continue
            if ff.pdb2gmx_name.lower() in self._pdb2gmx_index:
                continue
            self._force_fields[ff.name] = ff
            self._pdb2gmx_index[ff.pdb2gmx_name.lower()] = ff.name
            self._custom_force_fields[ff.name] = ff

    def _next_custom_id(self) -> str:
        """Allocate next custom force field id (a1, a2, ...)."""
        if not self._custom_force_fields:
            return "a1"
        max_id = max(self._custom_id_key(ff.name) for ff in self._custom_force_fields.values())
        return f"a{max_id + 1}"

    def list_force_fields(self) -> List[str]:
        """List all available force field names.

        Returns:
            List of force field names in order
        """
        builtin_names = [ff.name for ff in BUILTIN_FORCE_FIELDS]
        custom_names = [
            ff.name
            for ff in sorted(
                self._custom_force_fields.values(), key=lambda item: self._custom_id_key(item.name)
            )
        ]
        return builtin_names + custom_names

    def list_by_family(self, family: str) -> List[str]:
        """List force field names by family.

        Args:
            family: Family name (AMBER or CHARMM)

        Returns:
            List of force field names in the family (in order)
        """
        family = family.upper()
        builtin_names = [ff.name for ff in BUILTIN_FORCE_FIELDS if ff.family.upper() == family]
        custom_names = [
            ff.name
            for ff in sorted(
                self._custom_force_fields.values(), key=lambda item: self._custom_id_key(item.name)
            )
            if ff.family.upper() == family
        ]
        return builtin_names + custom_names

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

        # Handle custom short id format (e.g., "a1")
        if CUSTOM_ID_PATTERN.match(name_lower) and name_lower in self._force_fields:
            return self._force_fields[name_lower]

        # Try direct lookup first (for "1-a99SBdisp")
        if name_lower in self._force_fields:
            return self._force_fields[name_lower]

        # Try pdb2gmx name lookup (for "a99SBdisp", "charmm36m", etc.)
        if name_lower in self._pdb2gmx_index:
            cli_name = self._pdb2gmx_index[name_lower]
            return self._force_fields[cli_name]

        return None

    def list_custom_force_fields(self) -> List[str]:
        """List user-registered custom force field ids (a1, a2, ...)."""
        return [
            ff.name
            for ff in sorted(
                self._custom_force_fields.values(), key=lambda item: self._custom_id_key(item.name)
            )
        ]

    def register_custom_force_field(
        self,
        ff_dir: str,
        pdb2gmx_name: Optional[str] = None,
        family: str = "AMBER",
        water_model: str = "tip3p",
        solvate_cs: str = "spc216",
        gbsa_mapping: Optional[str] = None,
        description: str = "",
    ) -> ForceFieldInfo:
        """Register a custom all-atom force field and persist metadata.

        Args:
            ff_dir: Path to source *.ff directory
            pdb2gmx_name: Name used by gmx pdb2gmx (defaults to ff_dir stem)
            family: Force field family (AMBER or CHARMM)
            water_model: Water model for pdb2gmx
            solvate_cs: Water model for gmx solvate -cs
            gbsa_mapping: GBSA mapping identifier (optional, inferred from family if omitted)
            description: Human-readable description

        Returns:
            Registered ForceFieldInfo with allocated id (aN)
        """
        src_dir = Path(ff_dir).expanduser().resolve()
        if not src_dir.exists() or not src_dir.is_dir():
            raise ValueError(f"Force field directory does not exist: {src_dir}")
        if src_dir.suffix != ".ff":
            raise ValueError("Force field directory must end with '.ff'")
        if not (src_dir / "forcefield.itp").exists():
            raise ValueError(f"Missing required file: {(src_dir / 'forcefield.itp')}")

        family_normalized = family.upper()
        if family_normalized not in {"AMBER", "CHARMM"}:
            raise ValueError(f"Unsupported family: {family}. Must be AMBER or CHARMM.")
        if gbsa_mapping is None:
            gbsa_mapping = "AMBER99SB-ILDN" if family_normalized == "AMBER" else "CHARMM36"

        inferred_name = src_dir.name[:-3] if src_dir.name.endswith(".ff") else src_dir.name
        ff_pdb2gmx_name = (pdb2gmx_name or inferred_name).strip()
        if not ff_pdb2gmx_name:
            raise ValueError("pdb2gmx_name cannot be empty")
        if ff_pdb2gmx_name.lower() in self._pdb2gmx_index:
            existing = self._pdb2gmx_index[ff_pdb2gmx_name.lower()]
            raise ValueError(f"pdb2gmx_name '{ff_pdb2gmx_name}' already exists as '{existing}'")

        custom_id = self._next_custom_id()
        CUSTOM_FORCEFIELD_DIR.mkdir(parents=True, exist_ok=True)
        target_dir = CUSTOM_FORCEFIELD_DIR / f"{ff_pdb2gmx_name}.ff"
        if target_dir.exists():
            raise ValueError(f"Target custom force field directory already exists: {target_dir}")

        shutil.copytree(src_dir, target_dir)
        ff_info = ForceFieldInfo(
            name=custom_id,
            family=family_normalized,
            pdb2gmx_name=ff_pdb2gmx_name,
            water_model=water_model,
            solvate_cs=solvate_cs,
            gbsa_mapping=gbsa_mapping,
            description=description,
            source="custom",
            ff_dir=str(target_dir),
        )

        self._force_fields[ff_info.name.lower()] = ff_info
        self._pdb2gmx_index[ff_info.pdb2gmx_name.lower()] = ff_info.name.lower()
        self._custom_force_fields[ff_info.name.lower()] = ff_info
        self._save_custom_force_fields()
        return ff_info

    def remove_custom_force_field(self, name: str) -> ForceFieldInfo:
        """Remove a user-registered force field.

        Only custom aN ids are removable. Built-in force fields are read-only.
        """
        ff = self.get_force_field(name)
        if not ff:
            raise ValueError(f"Unknown force field: {name}")
        if ff.source != "custom":
            raise ValueError(f"Built-in force field cannot be removed: {ff.name}")

        ff_path = self.get_force_field_path(ff.name)
        key = ff.name.lower()
        self._force_fields.pop(key, None)
        self._custom_force_fields.pop(key, None)
        self._pdb2gmx_index.pop(ff.pdb2gmx_name.lower(), None)

        if ff_path and ff_path.exists():
            shutil.rmtree(ff_path)

        self._save_custom_force_fields()
        return ff

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

        # Custom force field path from explicit metadata
        if ff.ff_dir:
            ff_dir = Path(ff.ff_dir)
            if not ff_dir.is_absolute():
                ff_dir = Path(__file__).parent / ff_dir
            return ff_dir if ff_dir.exists() else None

        # Force field folder name matches pdb2gmx_name (e.g., "amber99sb-ildn" -> "amber99sb-ildn.ff")
        ff_folder_name = f"{ff.pdb2gmx_name}.ff"

        # Try to get the path using importlib.resources
        try:
            # For Python 3.9+
            ff_path = resources.files(forcefield) / ff_folder_name
            # Convert to Path if possible
            if hasattr(ff_path, "__fspath__"):
                return Path(str(ff_path))
            else:
                # It's a Traversable, try to resolve
                for p in resources.files(forcefield).iterdir():
                    if p.name == ff_folder_name:
                        return Path(str(p))
        except (ModuleNotFoundError, AttributeError):
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
        return (
            False,
            f"Unknown force field: {name}. Available: {', '.join(self.list_force_fields())}",
        )

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


def list_custom_force_fields() -> List[str]:
    """List user-registered custom force fields (aN ids)."""
    return REGISTRY.list_custom_force_fields()


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


def register_custom_force_field(
    ff_dir: str,
    pdb2gmx_name: Optional[str] = None,
    family: str = "AMBER",
    water_model: str = "tip3p",
    solvate_cs: str = "spc216",
    gbsa_mapping: Optional[str] = None,
    description: str = "",
) -> ForceFieldInfo:
    """Register a custom force field in global registry."""
    return REGISTRY.register_custom_force_field(
        ff_dir=ff_dir,
        pdb2gmx_name=pdb2gmx_name,
        family=family,
        water_model=water_model,
        solvate_cs=solvate_cs,
        gbsa_mapping=gbsa_mapping,
        description=description,
    )


def remove_custom_force_field(name: str) -> ForceFieldInfo:
    """Remove a custom force field from global registry."""
    return REGISTRY.remove_custom_force_field(name)


if __name__ == "__main__":
    # Print all available force fields
    print("Available force fields:")
    print("-" * 60)
    for ff in BUILTIN_FORCE_FIELDS:
        print(f"{ff.name:22} | {ff.family:6} | {ff.water_model:15} | {ff.gbsa_mapping}")
    print("-" * 60)
    print(f"Total: {len(BUILTIN_FORCE_FIELDS)} force fields")
