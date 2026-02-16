#!/usr/bin/env python3
"""
Coarse-Grained Simulation Module

A unified interface for coarse-grained molecular dynamics simulations
supporting multiple force fields (CALVADOS, HPS, MOFF, COCOMO, OpenMpipi).

Architecture:
- CGSimulationConfig: Configuration data class
- CGComponent: Individual component specification
- CGSimulator: Main simulator class with runner methods

Usage:
    from CondenSimAdapter.src import CGSimulationConfig, CGSimulator
    
    config = CGSimulationConfig.from_yaml("config.yaml")
    sim = CGSimulator(config)
    sim.setup("output/")
    sim.run_calvados(gpu_id=0)
"""

import os
import yaml
import shutil
import warnings
import numpy as np
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Any
from enum import Enum
from pathlib import Path

from openmm import Platform, LangevinIntegrator, LangevinMiddleIntegrator, XmlSerializer, Vec3
import openmm as mm
import openmm.unit as unit
from openmm.app import (
    Simulation, StateDataReporter, PDBFile
)
from openmm.unit import (
    picoseconds, nanometers, kilojoule, mole,
    Quantity, kilocalorie, amu, kelvin
)

# mdtraj reporters for trajectory saving with PBC support
from mdtraj.reporters import XTCReporter

from .pdb_tool import ChainLabel


# =============================================================================
# Enums
# =============================================================================

class ComponentType(Enum):
    """Component type."""
    IDP = "idp"   # Intrinsically disordered protein (sequence-based)
    MDP = "mdp"   # Folded protein (structure-based)


class TopologyType(Enum):
    """Topology type."""
    GRID = "grid"       # Continuous dense phase with periodic boundaries in x, y, z
    SLAB = "slab"       # Periodic boundaries in x, y and interfaces with dilute phase along z  
    DROPLET = "droplet" # Droplet geometry surrounded by dilute phase
    # Backward compatibility
    CUBIC = "grid"      # Alias for GRID


class ComputePlatform(Enum):
    """Compute platform."""
    CPU = "CPU"
    CUDA = "CUDA"


# =============================================================================
# Droplet Confinement Force
# =============================================================================

DROPLET_FORCE_K = 1.0
DROPLET_FORCE_STRIDE = 10


def _is_droplet_topology(topol: Any) -> bool:
    if isinstance(topol, TopologyType):
        return topol == TopologyType.DROPLET
    return str(topol).lower() == "droplet"


def _get_droplet_params(config: "CGSimulationConfig") -> tuple:
    box = config.box
    center = (box[0] / 2.0, box[1] / 2.0, box[2] / 2.0)
    radius = config.droplet_radius if config.droplet_radius is not None else box[0] / 2.0
    return radius, center


def add_droplet_force(
    system: mm.System,
    radius: float,
    center: tuple,
    k: float = DROPLET_FORCE_K,
    stride: int = DROPLET_FORCE_STRIDE,
) -> mm.System:
    """
    Add spherical confinement (droplet) force to the system.

    Energy expression:
        E = k * step(r - r0) * (r - r0)
        where r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
    """
    n_particles = system.getNumParticles()
    energy_expression = (
        "k * step(r - r0) * (r - r0); "
        "r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
    )

    confinement_force = mm.CustomExternalForce(energy_expression)
    confinement_force.addGlobalParameter("k", k)
    confinement_force.addGlobalParameter("r0", radius)
    confinement_force.addGlobalParameter("x0", center[0])
    confinement_force.addGlobalParameter("y0", center[1])
    confinement_force.addGlobalParameter("z0", center[2])

    for idx in range(0, n_particles, max(stride, 1)):
        confinement_force.addParticle(idx, [])

    system.addForce(confinement_force)
    return system


# =============================================================================
# Configuration Classes
# =============================================================================

@dataclass
class SimulationParams:
    """
    Core simulation parameters.

    Notes:
        dt and friction use fixed defaults and are not user inputs:
        - dt = 0.01 ps (10 fs) - common default across force fields
        - friction = 0.01 - OpenMM LangevinMiddleIntegrator default
    """
    # Fixed defaults (not user-editable)
    _DT: float = 0.01       # Time step (ps) - common 10 fs default
    _FRICTION: float = 0.01 # Friction - OpenMM default

    steps: int = 100000          # Total integration steps
    wfreq: int = 1000            # Write frequency (save every N steps)
    platform: ComputePlatform = ComputePlatform.CUDA
    verbose: bool = True

    def to_dict(self) -> Dict:
        d = {
            'steps': self.steps,
            'wfreq': self.wfreq,
            'platform': self.platform.value,
            'verbose': self.verbose,
        }
        # Drop None values
        return {k: v for k, v in d.items() if v is not None}

    @classmethod
    def from_dict(cls, d: Dict) -> 'SimulationParams':
        if 'platform' in d and isinstance(d['platform'], str):
            d['platform'] = ComputePlatform(d['platform'])
        # Drop None values
        d = {k: v for k, v in d.items() if v is not None}
        return cls(**d)


@dataclass
class BackmapConfig:
    """Backmap config section."""
    device: str = "cpu"  # Always use CPU
    model_type: Optional[str] = None  # None = auto, or "ResidueBasedModel"/"CalphaBasedModel"
    output_dir: Optional[str] = None  # Output directory, default {system_name}_backmap
    
    def to_dict(self) -> Dict:
        """Convert to dict."""
        d = {
            'device': self.device,
        }
        if self.model_type is not None:
            d['model_type'] = self.model_type
        if self.output_dir is not None:
            d['output_dir'] = self.output_dir
        return d
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'BackmapConfig':
        """Create from dict."""
        return cls(
            device=d.get('device', 'cpu'),
            model_type=d.get('model_type'),
            output_dir=d.get('output_dir'),
        )


@dataclass
class CGComponent:
    """
    Single component specification.
    
    Attributes:
        name: Component unique ID
        type: Component type (IDP or MDP)
        nmol: Molecule count for this component
        ffasta: FASTA file path (IDP)
        fpdb: PDB file path (MDP)
        fdomains: Domain definition file path (MDP)
        restraint: Whether to apply structural restraints (MDP)
        restraint_type: Restraint type (harmonic, go)
        use_com: Apply restraints to center of mass
        k_harmonic: Harmonic restraint force constant
        colabfold: PAE format (0=EBI, 1&2=Colabfold)
        charge_termini: Terminal charges (both, n, c, none)
    """
    name: str
    type: ComponentType
    nmol: int = 1
    
    # Input files
    ffasta: Optional[str] = None       # FASTA file (IDP)
    fpdb: Optional[str] = None         # PDB file (MDP)
    fdomains: Optional[str] = None     # Domain definition file (MDP)
    fpae: Optional[str] = None         # PAE JSON (Go-like potential)
    
    # Restraint settings
    restraint: bool = False
    restraint_type: str = "harmonic"
    use_com: bool = False
    k_harmonic: float = 700.0
    colabfold: int = 1
    
    # Terminal charge settings
    charge_termini: str = "both"
    
    # Derived attributes
    seq: Optional[str] = None
    nres: int = 0
    
    def to_dict(self) -> Dict:
        d = {
            'name': self.name,
            'type': self.type.value,
            'nmol': self.nmol,
            'restraint': self.restraint,
            'restraint_type': self.restraint_type,
            'use_com': self.use_com,
            'k_harmonic': self.k_harmonic,
            'colabfold': self.colabfold,
            'charge_termini': self.charge_termini,
        }
        if self.ffasta:
            d['ffasta'] = self.ffasta
        if self.fpdb:
            d['fpdb'] = self.fpdb
        if self.fdomains:
            d['fdomains'] = self.fdomains
        if self.fpae:
            d['fpae'] = self.fpae
        return d
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'CGComponent':
        comp_type = d.get('type', 'idp')
        if isinstance(comp_type, str):
            # Normalize to lower-case to match enum values
            comp_type = ComponentType(comp_type.lower())
        
        # For MDP type, default restraint to True if not explicitly specified
        # This is the expected behavior when restraint fields are hidden in YAML
        is_mdp = comp_type == ComponentType.MDP
        
        return cls(
            name=d['name'],
            type=comp_type,
            nmol=d.get('nmol', 1),
            ffasta=d.get('ffasta'),
            fpdb=d.get('fpdb'),
            fdomains=d.get('fdomains'),
            fpae=d.get('fpae'),
            restraint=d.get('restraint', True) if is_mdp else d.get('restraint', False),
            restraint_type=d.get('restraint_type', 'harmonic'),
            use_com=d.get('use_com', False),
            k_harmonic=d.get('k_harmonic', 700.0),
            colabfold=d.get('colabfold', 1),
            charge_termini=d.get('charge_termini', 'both'),
            # Derived attributes (if present in dict)
            seq=d.get('seq'),
            nres=d.get('nres', 0),
        )
    
    def validate(self) -> List[str]:
        """Validate configuration.

        fdomains supports two formats:
        1) File path: 'domains.yaml' (check file existence)
        2) Inline YAML: 'TDP43:\\n  - [3, 76]\\n...' (skip file check)
        """
        errors = []
        
        def _is_inline_yaml(text: str) -> bool:
            """Check whether fdomains is inline YAML (not a file path)."""
            if not text:
                return False
            stripped = text.strip()
            # Starts with { or [, or contains YAML-like structure with newlines
            if stripped.startswith('{') or stripped.startswith('['):
                return True
            if '\n' in stripped and (':' in stripped or stripped.startswith('-')):
                return True
            return False
        
        if self.type == ComponentType.IDP:
            if not self.ffasta:
                errors.append(f"Component '{self.name}': IDP requires ffasta file")
            elif not os.path.exists(self.ffasta):
                errors.append(f"Component '{self.name}': FASTA file not found: {self.ffasta}")
        elif self.type == ComponentType.MDP:
            if not self.fpdb:
                errors.append(f"Component '{self.name}': MDP requires fpdb file")
            elif not os.path.exists(self.fpdb):
                errors.append(f"Component '{self.name}': PDB file not found: {self.fpdb}")
            if self.restraint and self.fdomains:
                # Inline YAML does not require file existence check
                if not _is_inline_yaml(self.fdomains) and not os.path.exists(self.fdomains):
                    errors.append(f"Component '{self.name}': Domains file not found: {self.fdomains}")
        return errors


@dataclass
class CGSimulationConfig:
    """
    Full simulation configuration.
    
    Example YAML structure:
        system_name: my_simulation
        box: [25.0, 25.0, 30.0]
        temperature: 310.0
        ionic: 0.15
        topol: cubic  # or slab
        
        simulation:
          steps: 100000
          wfreq: 1000
          platform: CUDA
        
        components:
          - name: protein_A
            type: IDP
            nmol: 20
            ffasta: input/protein_A.fasta
    """
    # System info
    system_name: str = "cg_simulation"
    
    # CG force field
    force_field: Optional[str] = 'calvados'  # calvados | hps_urry | cocomo | mpipi_recharged
    
    # Environment
    box: List[float] = field(default_factory=lambda: [25.0, 25.0, 30.0])
    temperature: float = 310.0       # Kelvin
    ionic: float = 0.15              # Molar (ionic strength)
    
    # Topology
    topol: TopologyType = TopologyType.SLAB
    droplet_radius: Optional[float] = None
    
    # Simulation parameters
    simulation: SimulationParams = field(default_factory=SimulationParams)
    
    # Component list
    components: List[CGComponent] = field(default_factory=list)
    
    # Output
    output_dir: str = "output_cg"
    
    # Backmap config (optional)
    backmap: Optional[BackmapConfig] = None
    
    # Metadata
    config_path: Optional[str] = None
    created_at: str = field(default_factory=lambda: str(__import__('datetime').datetime.now()))
    
    def add_component(self, component: CGComponent):
        """Add a component."""
        self.components.append(component)
    
    def get_component(self, name: str) -> Optional[CGComponent]:
        """Get component by name."""
        for comp in self.components:
            if comp.name == name:
                return comp
        return None
    
    def total_molecules(self) -> int:
        """Compute total molecule count."""
        return sum(comp.nmol for comp in self.components)
    
    def validate(self) -> List[str]:
        """Validate configuration."""
        errors = []
        if not self.system_name:
            errors.append("system_name is required")
        if len(self.box) != 3:
            errors.append("box must be a list of 3 values [x, y, z]")
        if self.droplet_radius is not None and self.droplet_radius <= 0:
            errors.append("droplet radius must be positive")
        if not self.components:
            errors.append("At least one component is required")
        for comp in self.components:
            errors.extend(comp.validate())
        return errors
    
    def to_dict(self) -> Dict:
        """Convert to dict."""
        d = {
            'system_name': self.system_name,
            'force_field': self.force_field,
            'box': self.box,
            'temperature': self.temperature,
            'ionic': self.ionic,
            'topol': self.topol.value if isinstance(self.topol, TopologyType) else self.topol,
            'radius': self.droplet_radius,
            'simulation': self.simulation.to_dict(),
            'components': [c.to_dict() for c in self.components],
            'output_dir': self.output_dir,
        }
        if d.get('radius') is None:
            d.pop('radius', None)
        if self.backmap is not None:
            d['backmap'] = self.backmap.to_dict()
        return d
    
    def to_yaml(self, path: str = None):
        """Save to YAML file."""
        d = self.to_dict()
        if path:
            with open(path, 'w') as f:
                yaml.dump(d, f, default_flow_style=False, sort_keys=False)
        else:
            return yaml.dump(d, default_flow_style=False, sort_keys=False)
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'CGSimulationConfig':
        """Create from dict."""
        topol = d.get('topol', 'slab')
        if isinstance(topol, str):
            # Normalize to lower-case to match enum values
            topol = TopologyType(topol.lower())
        
        sim_dict = d.get('simulation', {})
        if isinstance(sim_dict, dict):
            simulation = SimulationParams.from_dict(sim_dict)
        else:
            simulation = SimulationParams()
        
        components = []
        for comp_dict in d.get('components', []):
            components.append(CGComponent.from_dict(comp_dict))
        
        # Handle backmap config
        backmap = None
        if 'backmap' in d and d['backmap']:
            backmap = BackmapConfig.from_dict(d['backmap'])
        
        droplet_radius = d.get('radius', d.get('droplet_radius'))
        
        return cls(
            system_name=d.get('system_name', 'cg_simulation'),
            force_field=d.get('force_field', 'calvados'),
            box=d.get('box', [25.0, 25.0, 30.0]),
            temperature=d.get('temperature', 310.0),
            ionic=d.get('ionic', 0.15),
            topol=topol,
            droplet_radius=droplet_radius,
            simulation=simulation,
            components=components,
            output_dir=d.get('output_dir', 'output_cg'),
            backmap=backmap,
            config_path=d.get('config_path'),
        )
    
    @classmethod
    def from_yaml(cls, path: str) -> 'CGSimulationConfig':
        """Load from YAML.

        Automatically convert relative paths to absolute paths based on the
        config file directory.
        """
        import os

        # Resolve config directory for relative paths
        config_dir = os.path.dirname(os.path.abspath(path))

        with open(path, 'r', encoding='utf-8') as f:
            d = yaml.safe_load(f)
        d['config_path'] = path

        config = cls.from_dict(d)

        # Convert component relative paths to absolute paths
        for comp in config.components:
            # Handle ffasta path
            if comp.ffasta and not os.path.isabs(comp.ffasta):
                comp.ffasta = os.path.join(config_dir, comp.ffasta)

            # Handle fpdb path
            if comp.fpdb and not os.path.isabs(comp.fpdb):
                comp.fpdb = os.path.join(config_dir, comp.fpdb)

            # Handle fdomains path (skip inline YAML)
            if comp.fdomains and not cls._is_inline_yaml(comp.fdomains):
                if not os.path.isabs(comp.fdomains):
                    comp.fdomains = os.path.join(config_dir, comp.fdomains)

            # Handle fpae path
            if comp.fpae and not os.path.isabs(comp.fpae):
                comp.fpae = os.path.join(config_dir, comp.fpae)

        return config

    @staticmethod
    def _is_inline_yaml(text: str) -> bool:
        """Check whether text is inline YAML (not a file path)."""
        if not text:
            return False
        stripped = text.strip()
        # Starts with { or [, or contains YAML-like structure with newlines
        if stripped.startswith('{') or stripped.startswith('['):
            return True
        if '\n' in stripped and (':' in stripped or stripped.startswith('-')):
            return True
        return False


# =============================================================================
# Simulation Result
# =============================================================================

@dataclass
class SimulationResult:
    """Simulation result."""
    success: bool = False
    output_dir: str = ""
    trajectory: Optional[str] = None
    structure: Optional[str] = None
    checkpoint: Optional[str] = None
    log: Optional[str] = None
    metrics: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)


@dataclass
class TopologyInfo:
    """
    Topology info cache.

    Stores full system topology metadata for analysis interfaces.

    Attributes:
        global_sequence: Global sequence (all residues concatenated)
        chain_ids: Chain ID per residue (numeric, starting at 1)
        is_folded: Folded-domain flag per residue (0=unfolded/IDP, 1=folded)
        molecule_indices: Molecule index per residue (1, 2, 3, ...)
        component_names: Component name per residue
        local_residue_indices: Residue index within a chain (1-based)
        sasa_values: SASA per residue (from all-atom structure)
    """
    global_sequence: str = ""
    chain_ids: List[int] = field(default_factory=list)
    is_folded: List[int] = field(default_factory=list)
    molecule_indices: List[int] = field(default_factory=list)
    component_names: List[str] = field(default_factory=list)
    local_residue_indices: List[int] = field(default_factory=list)
    sasa_values: List[float] = field(default_factory=list)


# =============================================================================
# CG Simulator with Multiple Runners
# =============================================================================

class CGSimulator:
    """
    Coarse-grained simulator.

    Unified interface across multiple force fields with runner methods.

    Attributes:
        config: Simulation configuration
        output_dir: Output directory
        is_setup: Whether setup is complete
        is_running: Whether a simulation is running
    """
    
    def __init__(self, config: CGSimulationConfig):
        """
        Initialize simulator.

        Args:
            config: CGSimulationConfig instance
        """
        self.config = config
        self.output_dir: Optional[str] = None
        self.is_setup: bool = False
        self.is_running: bool = False
        self._result: Optional[SimulationResult] = None
        self._topology_info: Optional[TopologyInfo] = None  # Topology cache

        # Validate configuration
        errors = self.config.validate()
        if errors:
            raise ValueError(
                f"Configuration validation failed:\n" +
                "\n".join(f"  - {e}" for e in errors)
            )

        print(f"[CGSimulator] Initialized")
        print(f"  System: {config.system_name}")
        print(f"  Components: {len(config.components)}")
        print(f"  Total molecules: {config.total_molecules()}")
    
    # -------------------------------------------------------------------------
    # Setup Methods
    # -------------------------------------------------------------------------

    def setup(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
        """
        Set up the simulation environment.

        Create output directory and prepare inputs.

        Args:
            output_dir: Output directory
            overwrite: Whether to overwrite existing directory

        Returns:
            Dict of generated file paths
        """
        self._ensure_not_running()

        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            if not overwrite:
                raise FileExistsError(
                    f"Output directory exists: {output_dir}\n"
                    f"Use overwrite=True to replace."
                )
            # If overwrite=True, keep directory; force-field-specific methods handle backups

        self.output_dir = output_dir

        print(f"\n[CGSimulator] Setting up...")
        print(f"  Output directory: {output_dir}")
        print(f"  System: {self.config.system_name}")

        self.is_setup = True
        print(f"  ✓ Setup complete")

        return {
            'output_dir': output_dir,
            'config': os.path.join(output_dir, 'config.yaml'),
        }

    # -------------------------------------------------------------------------
    # Topology / Sequence Interface Methods
    # -------------------------------------------------------------------------

    def get_composition(self) -> List[Dict[str, Any]]:
        """
        Return system composition summary.

        Returns:
            List of dicts, each containing:
            - name: Component name
            - nmol: Molecule count
            - sequence: Single-chain sequence
            - nres: Sequence length
            - type: Component type (IDP or MDP)
        """
        composition = []
        for comp in self.config.components:
            # Get sequence
            seq = self._get_component_sequence(comp)
            composition.append({
                'name': comp.name,
                'nmol': comp.nmol,
                'sequence': seq,
                'nres': len(seq),
                'type': comp.type.value,
            })
        return composition

    def get_global_sequence(self) -> str:
        """
        Return the global sequence (concatenate all chains).

        Returns:
            Global sequence string
        """
        info = self._get_topology_info()
        return info.global_sequence

    def get_chain_ids(self) -> List[int]:
        """
        Return chain IDs aligned with the global sequence.

        Chain IDs are numeric and start at 1.

        Returns:
            Chain ID list (numeric, 1-based)
        """
        info = self._get_topology_info()
        return info.chain_ids

    def get_folded_domains(self) -> List[int]:
        """
        Return folded-domain flags aligned with the global sequence.

        Returns:
            Integer list: 0 = IDP/unfolded, 1 = folded domain
        """
        info = self._get_topology_info()
        return info.is_folded

    def get_chain_identifiers(self) -> List[str]:
        """
        Return chain identifiers aligned with the global sequence.

        These identifiers are suitable for OpenMM topology.
        Format: 'A', 'B', ..., 'Z', 'A1', 'B1', ..., 'Z1', 'A2', ...
        Supports unlimited chains.

        Returns:
            Chain identifier list aligned with the global sequence
        """
        info = self._get_topology_info()
        chain_ids = info.chain_ids

        def chain_id_to_identifier(chain_id: int) -> str:
            """Convert chain ID to unique identifier."""
            if chain_id <= 26:
                return chr(ord('A') + chain_id - 1)
            else:
                letter_idx = (chain_id - 1) % 26
                suffix = (chain_id - 1) // 26
                return chr(ord('A') + letter_idx) + str(suffix)

        return [chain_id_to_identifier(cid) for cid in chain_ids]

    def get_unique_chain_identifiers(self) -> List[str]:
        """
        Return unique chain identifiers.

        Returns:
            Unique identifiers, sorted alphabetically
        """
        identifiers = set(self.get_chain_identifiers())
        return sorted(list(identifiers))

    def _build_component_names(self) -> List[str]:
        """
        Build component name list per residue.

        Returns:
            Component name list aligned with global sequence
        """
        comp_names = []
        for comp in self.config.components:
            for _ in range(comp.nmol):
                comp_names.extend([comp.name] * len(self._get_component_sequence(comp)))
        return comp_names

    def _compute_sasa_for_component(self, comp: CGComponent) -> List[float]:
        """
        Compute SASA values for an MDP component.

        Uses mdsim to compute per-residue SASA from all-atom PDB.
        For IDP components, returns default SASA values.

        Args:
            comp: CGComponent instance

        Returns:
            SASA list (one value per residue)
        """
        # IDP uses default SASA values
        if comp.type == ComponentType.IDP:
            nres = len(self._get_component_sequence(comp))
            return [5.0] * nres

        # MDP requires all-atom PDB
        if comp.type == ComponentType.MDP:
            if not comp.fpdb:
                raise ValueError(f"Component '{comp.name}' is MDP but no fpdb file specified")

            # Try importing system_handling for SASA computation
            try:
                from CondenSimAdapter.extern.cocomo.src.cocomo.system_handling import ComponentType as SASAComponentType
            except ImportError:
                warnings.warn(
                    f"Could not import cocomo.system_handling for SASA calculation. "
                    f"Using default SASA values for component '{comp.name}'.",
                    UserWarning
                )
                nres = len(self._get_component_sequence(comp))
                return [5.0] * nres

            try:
                # Create ComponentType and compute SASA automatically
                ctype = SASAComponentType(
                    name=comp.name,
                    pdb=comp.fpdb,
                    getsasa="auto",  # Use mdsim automatic SASA
                    mask_sasa_bydomain=True  # Mask by domain
                )

                # Extract SASA values
                sasa_list = [val for _, val in ctype.sasa]

                # Ensure SASA length matches sequence length
                nres = len(self._get_component_sequence(comp))
                if len(sasa_list) != nres:
                    warnings.warn(
                        f"SASA length ({len(sasa_list)}) doesn't match sequence length ({nres}) "
                        f"for component '{comp.name}'. Padding with default values.",
                        UserWarning
                    )
                    # Pad or truncate to correct length
                    if len(sasa_list) < nres:
                        sasa_list.extend([5.0] * (nres - len(sasa_list)))
                    else:
                        sasa_list = sasa_list[:nres]

                return sasa_list

            except Exception as e:
                warnings.warn(
                    f"Failed to compute SASA for component '{comp.name}': {e}. "
                    f"Using default SASA values.",
                    UserWarning
                )
                nres = len(self._get_component_sequence(comp))
                return [5.0] * nres

        # Default case
        nres = len(self._get_component_sequence(comp))
        return [5.0] * nres

    def _compute_all_sasa_values(self) -> List[float]:
        """
        Compute SASA values for the entire system.

        Compute per-component SASA and repeat by molecule count.

        Returns:
            SASA list aligned with the global sequence
        """
        all_sasa = []

        for comp in self.config.components:
            # Compute SASA for a single chain
            single_sasa = self._compute_sasa_for_component(comp)

            # Repeat by molecule count
            for _ in range(comp.nmol):
                all_sasa.extend(single_sasa)

        return all_sasa

    def _get_sasa_values(self) -> Optional[np.ndarray]:
        """
        Try to load or compute SASA values.

        Prefer loading from file; compute if not found.
        For MDP components, use mdsim on all-atom PDB.
        For IDP components, use default SASA values.

        Returns:
            SASA array, or None if unavailable
        """
        # Try cached file (sasa_mdsim.txt from previous mdsim computation)
        sasa_mdsim_file = os.path.join(self.output_dir or self.config.output_dir, 'sasa_mdsim.txt')
        if os.path.exists(sasa_mdsim_file):
            print(f"  Loading SASA data: {sasa_mdsim_file}")
            return np.loadtxt(sasa_mdsim_file)

        # Try legacy cached file
        surface_file = os.path.join(self.output_dir or self.config.output_dir, 'surface')
        if os.path.exists(surface_file):
            print(f"  Loading SASA data: {surface_file}")
            return np.loadtxt(surface_file)

        # Try alternative path
        alt_surface_file = os.path.join(self.output_dir or self.config.output_dir, 'sasa_values.txt')
        if os.path.exists(alt_surface_file):
            print(f"  Loading SASA data: {alt_surface_file}")
            return np.loadtxt(alt_surface_file)

        # Compute if no cached file
        try:
            sasa = self._compute_all_sasa_values()
            if sasa:
                print(f"  Calculated {len(sasa)} SASA values")

                # Save mdsim-computed SASA values as plain text
                np.savetxt(sasa_mdsim_file, sasa, fmt='%.4f')
                print(f"  Saved mdsim SASA values to: {sasa_mdsim_file}")

                return np.array(sasa)
        except Exception as e:
            print(f"  SASA calculation failed: {e}")

        print(f"  SASA file not found, using default values")
        return None

    def get_sasa_values(self) -> List[float]:
        """
        Get SASA values for the entire system.

        Returns a list aligned with the global sequence.
        For MDP components, computed from all-atom PDB;
        for IDP components, use default SASA values.

        Returns:
            SASA value list
        """
        info = self._get_topology_info()
        if info.sasa_values:
            return info.sasa_values

        # Compute if topology cache is empty
        sasa = self._compute_all_sasa_values()

        # Cache into topology info
        self._topology_info.sasa_values = sasa

        return sasa

    def _get_topology_info(self) -> TopologyInfo:
        """
        Get topology info (lazy initialization, cached).

        Returns:
            TopologyInfo instance
        """
        if self._topology_info is None:
            self._topology_info = self._build_topology_info()
        return self._topology_info

    def _build_topology_info(self) -> TopologyInfo:
        """
        Build topology info.

        Builds global sequence, chain IDs, folded-domain flags.
        For MDP components, compute SASA from all-atom PDB.

        Returns:
            TopologyInfo instance
        """
        global_sequence_parts = []
        chain_ids = []
        is_folded = []
        molecule_indices = []
        component_names = []
        local_residue_indices = []
        sasa_values = []

        # Current chain ID (1-based)
        current_chain_id = 1

        # Precompute SASA per component (once)
        component_sasa = {}
        for comp in self.config.components:
            component_sasa[comp.name] = self._compute_sasa_for_component(comp)

        for comp in self.config.components:
            # Single-chain sequence
            single_seq = self._get_component_sequence(comp)
            nres = len(single_seq)
            nmol = comp.nmol

            # Folded-domain flags for one chain
            single_folded = self._get_component_folded_domains(comp, nres)

            # SASA values for this component
            single_sasa = component_sasa.get(comp.name, [5.0] * nres)

            # Build per-molecule info
            for mol_idx in range(nmol):
                # Append sequence
                global_sequence_parts.append(single_seq)

                # Append per-residue info
                for res_idx in range(nres):
                    chain_ids.append(current_chain_id)
                    is_folded.append(single_folded[res_idx])
                    molecule_indices.append(current_chain_id)  # chain ID is molecule ID
                    component_names.append(comp.name)
                    local_residue_indices.append(res_idx + 1)  # 1-based
                    sasa_values.append(single_sasa[res_idx])

                current_chain_id += 1

        return TopologyInfo(
            global_sequence="".join(global_sequence_parts),
            chain_ids=chain_ids,
            is_folded=is_folded,
            molecule_indices=molecule_indices,
            component_names=component_names,
            local_residue_indices=local_residue_indices,
            sasa_values=sasa_values,
        )

    def _get_component_sequence(self, comp: CGComponent) -> str:
        """
        Get single-chain sequence for component

        Args:
            comp: CGComponent instance

        Returns:
            Single-chain sequence string
        """
        # If sequence already exists, return directly
        if comp.seq:
            return comp.seq

        # Read from FASTA file (IDP) or PDB file (MDP)
        if comp.type == ComponentType.IDP:
            if comp.ffasta:
                return self._read_fasta(comp.ffasta, component_name=comp.name)
            else:
                raise ValueError(f"Component '{comp.name}' is IDP but no ffasta file specified")
        elif comp.type == ComponentType.MDP:
            if comp.fpdb:
                return self._seq_from_pdb(comp.fpdb)
            else:
                raise ValueError(f"Component '{comp.name}' is MDP but no fpdb file specified")

        return ""

    def _read_fasta(self, fasta_path: str, component_name: str = None) -> str:
        """
        Read sequence from FASTA file

        If component_name is specified, return the sequence matching that name.
        If not specified or no match found, return the first sequence.

        Args:
            fasta_path: FASTA file path
            component_name: Component name (used for sequence selection)

        Returns:
            Sequence string
        """
        from Bio import SeqIO

        # Handle relative paths
        if not os.path.isabs(fasta_path):
            # Relative to current working directory or config file directory
            if hasattr(self.config, 'config_path') and self.config.config_path:
                config_dir = os.path.dirname(os.path.abspath(self.config.config_path))
                fasta_path = os.path.join(config_dir, fasta_path)

        # Read FASTA file
        records = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            raise ValueError(f"Empty or invalid FASTA file: {fasta_path}")

        # If component name is specified, try to match
        if component_name:
            if component_name in records:
                return str(records[component_name].seq)
            else:
                # Try case-insensitive matching
                for name in records:
                    if name.lower() == component_name.lower():
                        return str(records[name].seq)
                # No match found, use first sequence with warning
                print(f"  [WARNING] Component '{component_name}' not found in fasta, using first sequence")

        # Return first sequence
        return str(list(records.values())[0].seq)

    def _seq_from_pdb(self, pdb_path: str) -> str:
        """
        Extract sequence from PDB file

        Args:
            pdb_path: PDB file path

        Returns:
            Sequence string
        """
        # Handle relative paths
        if not os.path.isabs(pdb_path):
            # Relative to current working directory or config file directory
            if hasattr(self.config, 'config_path') and self.config.config_path:
                config_dir = os.path.dirname(os.path.abspath(self.config.config_path))
                pdb_path = os.path.join(config_dir, pdb_path)

        # Use MDAnalysis to extract sequence
        try:
            from MDAnalysis import Universe
        except ImportError:
            raise ImportError("MDAnalysis is required to extract sequences from PDB files")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = Universe(pdb_path)

            # Get unique residues
            residues = u.residues
            n_res = len(residues)

            # 3-letter to 1-letter amino acid mapping
            aa_3to1 = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                'SEC': 'U', 'PYL': 'O',  # Selenocysteine, Pyrrolysine
            }

            fastapdb = ""
            for res in residues:
                resname = res.resname
                fastapdb += aa_3to1.get(resname, 'X')

            return fastapdb

    def _get_component_folded_domains(self, comp: CGComponent, nres: int) -> List[int]:
        """
        Get folded domain information for component

        Args:
            comp: CGComponent instance
            nres: Sequence length

        Returns:
            List of length nres, 0=unfolded/IDP, 1=folded domain
        """
        # IDP defaults to all 0
        if comp.type == ComponentType.IDP:
            return [0] * nres

        # MDP check if fdomains is configured
        if not comp.fdomains:
            return [0] * nres

        # Parse fdomains YAML
        domains = self._parse_fdomains(comp.fdomains)

        # Build folded array
        folded = [0] * nres
        for (start, end) in domains:
            # Ensure within valid range
            start = max(1, start)  # 1-based
            end = min(nres, end)
            if start <= end:
                for i in range(start - 1, end):  # Convert to 0-based
                    folded[i] = 1

        return folded

    def _parse_fdomains(self, fdomains: str) -> List[tuple]:
        """
        Parse fdomains configuration

        Supports two formats:
        1. File path: parse YAML file
        2. Inline YAML: parse string directly

        Args:
            fdomains: fdomains configuration

        Returns:
            List of domains, each domain is (start, end) tuple (1-based)
        """
        # Check if it's inline YAML
        if self.config._is_inline_yaml(fdomains):
            # Parse string directly
            data = yaml.safe_load(fdomains)
        else:
            # Parse file
            fdomains_abs = fdomains
            if not os.path.isabs(fdomains):
                if hasattr(self.config, 'config_path') and self.config.config_path:
                    config_dir = os.path.dirname(os.path.abspath(self.config.config_path))
                    fdomains_abs = os.path.join(config_dir, fdomains)

            if not os.path.exists(fdomains_abs):
                return []

            with open(fdomains_abs, 'r') as f:
                data = yaml.safe_load(f)

        # Parse domain definitions
        domains = []
        if isinstance(data, dict):
            for protein_name, domain_list in data.items():
                if isinstance(domain_list, list):
                    for domain in domain_list:
                        if isinstance(domain, (list, tuple)) and len(domain) == 2:
                            domains.append((domain[0], domain[1]))

        return domains

    def clear_topology_cache(self):
        """
        Clear topology information cache

        Will be recalculated when interface methods are called next time.
        """
        self._topology_info = None

    def prepare_calvados_output(self) -> Dict[str, str]:
        """
        Prepare CALVADOS output directory structure

        Unified output structure:
        {output_dir}/
        ├── {system_name}_CG/
        │   ├── raw/                  # Native output
        │   ├── trajectory.xtc        # Organized trajectory
        │   ├── final.pdb             # Organized final structure
        │   └── simulation.log        # High-level log

        Returns:
            Dictionary containing output paths
        """
        self._ensure_setup()
        self._ensure_not_running()

        # Check if self.output_dir already contains _CG suffix
        expected_suffix = f"{self.config.system_name}_CG"
        if self.output_dir.endswith(expected_suffix):
            # Already contains _CG suffix, use directly
            output_dir = self.output_dir
            task_name = expected_suffix
        else:
            # Add _CG suffix
            task_name = expected_suffix
            output_dir = os.path.join(self.output_dir, task_name)

        raw_dir = os.path.join(output_dir, 'raw')

        # If directory exists, backup and recreate
        import shutil
        from datetime import datetime

        if os.path.exists(output_dir):
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{output_dir}_backup_{timestamp}"
            shutil.move(output_dir, backup_dir)
            print(f"  Backup old results to: {backup_dir}")

        # Create parent directory first, then raw subdirectory
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(raw_dir, exist_ok=True)

        return {
            'output_dir': output_dir,
            'raw_dir': raw_dir,
            'task_name': task_name,
        }
    
    def _copy_input_files(self, output_dir: str):
        """Copy input files to output directory"""
        input_dir = os.path.join(output_dir, 'input')
        os.makedirs(input_dir, exist_ok=True)
        
        for comp in self.config.components:
            if comp.type == ComponentType.IDP and comp.ffasta:
                if os.path.exists(comp.ffasta):
                    shutil.copy2(comp.ffasta, os.path.join(input_dir, os.path.basename(comp.ffasta)))
            
            elif comp.type == ComponentType.MDP:
                if comp.fpdb and os.path.exists(comp.fpdb):
                    shutil.copy2(comp.fpdb, os.path.join(input_dir, os.path.basename(comp.fpdb)))
                if comp.fdomains and os.path.exists(comp.fdomains):
                    shutil.copy2(comp.fdomains, os.path.join(input_dir, os.path.basename(comp.fdomains)))
                if comp.fpae and os.path.exists(comp.fpae):
                    shutil.copy2(comp.fpae, os.path.join(input_dir, os.path.basename(comp.fpae)))
    
    def _ensure_setup(self):
        """Ensure setup is complete"""
        if not self.is_setup:
            raise RuntimeError("Simulation not set up. Call setup() first.")
    
    def _ensure_not_running(self):
        """Ensure not running"""
        if self.is_running:
            raise RuntimeError("Simulation is already running")

    # -------------------------------------------------------------------------
    # Pre-equilibration Methods (Using CALVADOS to Build Initial Structure)
    # -------------------------------------------------------------------------
    # Note: All non-CALVADOS runners automatically call this method for pre-equilibration
    # Pre-equilibration parameters are hardcoded, but each runner can have different defaults
    # -------------------------------------------------------------------------

    def _run_pre_equilibration(
        self,
        gpu_id: int = 0,
        steps: int = 100000,
        mapping: str = "ca",
        k_restraint: float = 10000.0,
        use_com: bool = True,
        platform: ComputePlatform = ComputePlatform.CUDA,
    ) -> Optional[str]:
        """
        Run pre-equilibration (Using CALVADOS to Build Initial Structure)

        Purpose of pre-equilibration:
        1. Build initial CG structure using CALVADOS force field
        2. Run short simulation for MDP protein (with restraints) to adapt to CG representation
        3. Run short simulation for IDP protein (without restraints) to generate initial CG structure
        4. Provide good initial structure for subsequent force field simulations (COCOMO/HPS/MOFF/OpenMpipi)

        Mapping options for MDP protein:
        - CA (alpha carbon): Use alpha-carbon coordinates, suitable for scenarios requiring backbone information
        - COM (center of mass): Use residue center of mass, suitable for smoother mapping

        Args:
            gpu_id: GPU device ID
            steps: Pre-equilibration steps (default 100k steps)
            mapping: Mapping method ('ca' or 'com')
            k_restraint: Restraint force constant (kJ/(mol·nm²))
            use_com: Whether to use COM restraint
            platform: Computing platform (CUDA or CPU)

        Returns:
            Path to final.pdb after pre-equilibration, or None if no components
        """
        from .calvados_wrapper import CalvadosWrapper
        import shutil

        self._ensure_setup()
        self._ensure_not_running()

        # Check if there are components to process
        has_components = len(self.config.components) > 0
        has_mdp = any(comp.type == ComponentType.MDP for comp in self.config.components)
        if not has_components:
            return None

        # Print restraint info for MDP systems, skip for IDP systems
        if has_mdp:
            print(f"\n[Pre-equilibration] Using CALVADOS to Build Initial Structure (MDP protein, with restraints)...")
            print(f"  Mapping method: {mapping.upper()}")
            print(f"  Restraint force constant: {k_restraint} kJ/(mol·nm²)")
            print(f"  Use COM restraint: {'Yes' if use_com else 'No'}")
        else:
            print(f"\n[Pre-equilibration] Using CALVADOS to Build Initial Structure (Pure IDP system)...")
            print(f"  Mapping method: {mapping.upper()}")
            print(f"  Pre-equilibration steps: {steps}")

        # Save original component configuration
        original_components = []
        for comp in self.config.components:
            original_components.append({
                'restraint': comp.restraint,
                'restraint_type': comp.restraint_type,
                'use_com': comp.use_com,
                'k_harmonic': comp.k_harmonic,
            })

        # Create temporary configuration for pre-equilibration
        temp_config = self._create_preequil_config(
            steps=steps,
            mapping=mapping,
            k_restraint=k_restraint,
            use_com=use_com,
            platform=platform,
        )

        # Output directory structure: {output_dir}/{system_name}_CG/equilibration/
        # Consistent with run_calvados
        task_name = f"{self.config.system_name}_CG"
        equilibration_dir = os.path.join(self.output_dir, task_name, 'equilibration')
        raw_dir = os.path.join(equilibration_dir, 'raw')

        # If directory exists, backup first
        if os.path.exists(equilibration_dir):
            import time
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{equilibration_dir}_backup_{timestamp}"
            shutil.move(equilibration_dir, backup_dir)
            print(f"  Backup old pre-equilibration results to: {backup_dir}")

        os.makedirs(raw_dir, exist_ok=True)

        try:
            # GPU selection is controlled by DeviceIndex attribute in sim.py
            # No need to set CUDA_VISIBLE_DEVICES environment variable anymore

            # Call CalvadosWrapper
            wrapper = CalvadosWrapper(temp_config)

            # Write config files to raw directory (pass verbose)
            verbose = self.config.simulation.verbose
            files = wrapper._write_to_dir(raw_dir, gpu_id=gpu_id, verbose=verbose)
            print(f"  Config files written to: {files['components']}")
            print(f"  Platform: {platform.value}")

            # Run CALVADOS simulation
            from CondenSimAdapter.extern.ms2_calvados.calvados import sim as calvados_sim
            try:
                calvados_sim.run(
                    path=raw_dir,
                    fconfig='config.yaml',
                    fcomponents='components.yaml'
                )
            except Exception as calvados_error:
                # If CUDA fails, try using CPU
                if 'CUDA' in str(calvados_error) or 'Platform' in str(calvados_error) or 'no registered Platform' in str(calvados_error):
                    print(f"  CALVADOS pre-equilibration CUDA failed, trying CPU...")
                    # Modify platform in config file to CPU
                    import yaml
                    config_file = os.path.join(raw_dir, 'config.yaml')
                    with open(config_file, 'r') as f:
                        config_dict = yaml.safe_load(f)
                    config_dict['platform'] = 'CPU'
                    with open(config_file, 'w') as f:
                        yaml.dump(config_dict, f)
                    
                    # Re-run
                    calvados_sim.run(
                        path=raw_dir,
                        fconfig='config.yaml',
                        fcomponents='components.yaml'
                    )
                else:
                    raise

            # Find the generated final structure
            final_pdb = os.path.join(equilibration_dir, 'final.pdb')
            if os.path.exists(os.path.join(raw_dir, 'checkpoint.pdb')):
                shutil.copy2(
                    os.path.join(raw_dir, 'checkpoint.pdb'),
                    final_pdb
                )
            else:
                # Find timestamped PDB files
                for f in os.listdir(raw_dir):
                    if f.endswith('.pdb') and f != 'top.pdb':
                        shutil.copy2(
                            os.path.join(raw_dir, f),
                            final_pdb
                        )
                        break

            # Copy to output directory root (convenient for subsequent force fields)
            output_pdb = os.path.join(self.output_dir, 'preequil_final.pdb')
            if os.path.exists(final_pdb):
                shutil.copy2(final_pdb, output_pdb)
                print(f"  Pre-equilibration complete: {output_pdb}")
                print(f"  Pre-equilibration output: {equilibration_dir}")
                return output_pdb
            else:
                print(f"  Pre-equilibration output file not found")
                return None

        except Exception as e:
            print(f"  Pre-equilibration failed: {e}")
            import traceback
            traceback.print_exc()
            return None

        finally:
            # Restore original component configuration
            for i, comp in enumerate(self.config.components):
                orig = original_components[i]
                comp.restraint = orig['restraint']
                comp.restraint_type = orig['restraint_type']
                comp.use_com = orig['use_com']
                comp.k_harmonic = orig['k_harmonic']

    def _create_preequil_config(
        self,
        steps: int = 100000,
        mapping: str = "ca",
        k_restraint: float = 10000.0,
        use_com: bool = True,
        platform: ComputePlatform = ComputePlatform.CUDA,
    ) -> 'CGSimulationConfig':
        """
        Create temporary configuration for pre-equilibration

        For MDP components:
        - Enable restraint
        - Set restraint_type to harmonic
        - Set use_com based on mapping:
          - CA: use_com=False (restrain to CA atom of each residue)
          - COM: use_com=True (restrain to center of mass of each residue)
        - Set k_harmonic to k_restraint
        - Use steps for steps
        - Set platform

        Args:
            steps: Pre-equilibration steps
            mapping: Mapping method ('ca' or 'com')
            k_restraint: Restraint force constant
            use_com: Whether to use COM restraint
            platform: Computing platform (CUDA or CPU)

        Returns:
            Temporary configuration object
        """
        from copy import deepcopy

        # Deep copy configuration
        temp_config = deepcopy(self.config)

        # Temporarily modify MDP component restraint settings
        for comp in temp_config.components:
            if comp.type == ComponentType.MDP:
                comp.restraint = True
                comp.restraint_type = 'harmonic'
                comp.use_com = use_com
                comp.k_harmonic = k_restraint

        # Modify simulation parameters to pre-equilibration parameters
        temp_config.simulation = deepcopy(self.config.simulation)
        temp_config.simulation.steps = steps
        temp_config.simulation.wfreq = min(steps // 10, 1000)
        temp_config.simulation.platform = platform

        return temp_config

    def get_pre_equilibrated_structure(self) -> Optional[str]:
        """
        Get the path to the pre-equilibrated structure file

        Returns:
            Path to pre-equilibration structure file, or None if pre-equilibration was not run
        """
        if self.output_dir is None:
            return None

        preequil_pdb = os.path.join(self.output_dir, 'preequil_final.pdb')
        if os.path.exists(preequil_pdb):
            return preequil_pdb
        return None

    # -------------------------------------------------------------------------
    # Runner Methods
    # -------------------------------------------------------------------------
    
    def run_calvados(self, gpu_id: int = 0, continue_from: str = None, **kwargs) -> SimulationResult:
        """
        Run CALVADOS simulation

        Directly use CALVADOS-generated timestamped PDB files as final output.

        Args:
            gpu_id: GPU device ID
            continue_from: PDB file path for continuing simulation, supports resuming from specified structure
            **kwargs: Additional parameters

        Returns:
            SimulationResult
        """
        from .calvados_wrapper import CalvadosWrapper
        import time
        from datetime import datetime

        self._ensure_setup()
        self._ensure_not_running()

        # Prepare output directory
        dirs = self.prepare_calvados_output()
        output_dir = dirs['output_dir']
        raw_dir = dirs['raw_dir']
        task_name = dirs['task_name']

        # Copy input files (after backup logic)
        self._copy_input_files(output_dir)

        self.is_running = True
        result = SimulationResult()
        result.output_dir = output_dir

        try:
            print(f"\n[CALVADOS] Running simulation via CGSimulator...")
            print(f"  GPU ID: {gpu_id}")
            print(f"  Task: {task_name}")
            print(f"  Raw output: {raw_dir}")
            print(f"  Topology: {self.config.topol.value if hasattr(self.config.topol, 'value') else self.config.topol}")

            # GPU selection is controlled by DeviceIndex attribute in sim.py
            # No need to set CUDA_VISIBLE_DEVICES environment variable anymore

            # Call CalvadosWrapper config writing and simulation run (pass raw_dir, gpu_id, and verbose)
            verbose = self.config.simulation.verbose
            wrapper = CalvadosWrapper(self.config)
            
            # Handle continue_from parameter
            continue_from_file = None
            if continue_from:
                if not os.path.exists(continue_from):
                    raise FileNotFoundError(f"Continue from file not found: {continue_from}")
                # Copy coordinate file to raw directory
                continue_from_file = os.path.basename(continue_from)
                continue_from_dst = os.path.join(raw_dir, continue_from_file)
                import shutil
                shutil.copy2(continue_from, continue_from_dst)
                print(f"  [Continue From] Using coordinates from: {continue_from}")
                print(f"  [Continue From] Copied to: {continue_from_dst}")
            
            wrapper._write_to_dir(raw_dir, gpu_id=gpu_id, verbose=verbose, continue_from=continue_from_file)

            # Run simulation
            from CondenSimAdapter.extern.ms2_calvados.calvados import sim as calvados_sim
            calvados_sim.run(
                path=raw_dir,
                fconfig='config.yaml',
                fcomponents='components.yaml'
            )

            # Copy trajectory file
            self._copy_trajectory(raw_dir, output_dir)

            # Copy timestamped PDB file directly as final.pdb
            self._copy_final_pdb(raw_dir, output_dir)

            # Write log
            elapsed = 0  # TODO: track actual time
            self._write_simulation_log(output_dir, task_name, elapsed, True)

            result.success = True
            print(f"  ✓ CALVADOS simulation completed")

        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            print(f"  ✗ CALVADOS simulation failed: {e}")

        finally:
            self.is_running = False

        # Set result file paths
        result.trajectory = os.path.join(output_dir, 'trajectory.xtc')
        result.structure = os.path.join(output_dir, 'final.pdb')

        for key in ['trajectory', 'structure']:
            path = getattr(result, key)
            if path and not os.path.exists(path):
                setattr(result, key, None)

        self._result = result
        return result

    def _copy_trajectory(self, raw_dir: str, output_dir: str):
        """Copy trajectory file from raw directory"""
        import shutil
        sysname = self.config.system_name
        src_xtc = os.path.join(raw_dir, f'{sysname}.xtc')
        dst_xtc = os.path.join(output_dir, 'trajectory.xtc')
        if os.path.exists(src_xtc):
            shutil.copy2(src_xtc, dst_xtc)
            print(f"  📦 trajectory.xtc")

    def _copy_final_pdb(self, raw_dir: str, output_dir: str):
        """
        Copy timestamped PDB file from raw directory to output_dir as final.pdb

        CALVADOS generates PDB files in format {system_name}_{timestamp}.pdb
        """
        import glob

        # Find timestamped PDB files
        pattern = os.path.join(raw_dir, f'{self.config.system_name}_*.pdb')
        pdb_files = glob.glob(pattern)

        if not pdb_files:
            print(f"  ⚠ No timestamped PDB file found in {raw_dir}")
            return

        # Find latest file (sorted by modification time)
        latest_pdb = max(pdb_files, key=os.path.getmtime)

        # Copy to output_dir/final.pdb
        dst_pdb = os.path.join(output_dir, 'final.pdb')
        shutil.copy2(latest_pdb, dst_pdb)
        print(f"  📦 final.pdb (copied from {os.path.basename(latest_pdb)})")

    def _organize_calvados_output(self, raw_dir: str, output_dir: str, task_name: str):
        """
        Organize CALVADOS output files to unified structure

        Unified naming rules:
        - trajectory.xtc  <- {system_name}.xtc
        - final.pdb       <- checkpoint.pdb or timestamped PDB
        """
        import shutil

        sysname = self.config.system_name

        # 1. Process trajectory file
        src_xtc = os.path.join(raw_dir, f'{sysname}.xtc')
        dst_xtc = os.path.join(output_dir, 'trajectory.xtc')
        if os.path.exists(src_xtc):
            shutil.copy2(src_xtc, dst_xtc)
            print(f"  📦 trajectory.xtc")

        # 2. Find and copy final structure
        src_pdb = os.path.join(raw_dir, 'checkpoint.pdb')
        if not os.path.exists(src_pdb):
            for f in os.listdir(raw_dir):
                if f.endswith('.pdb') and f != 'top.pdb':
                    src_pdb = os.path.join(raw_dir, f)
                    break

        dst_pdb = os.path.join(output_dir, 'final.pdb')
        if os.path.exists(src_pdb):
            shutil.copy2(src_pdb, dst_pdb)
            print(f"  📦 final.pdb")

        # 3. Copy important files
        important_files = [
            (f'{sysname}.xml', 'system.xml'),
            ('top.pdb', 'top.pdb'),
            ('restart.chk', 'restart.chk'),
            ('checkpoint.pdb', 'checkpoint.pdb'),
        ]
        for src_name, dst_name in important_files:
            src = os.path.join(raw_dir, src_name)
            if os.path.exists(src):
                dst = os.path.join(raw_dir, dst_name)
                if src != dst:
                    shutil.copy2(src, dst)

        print(f"  Raw output organized to: {raw_dir}")

    def _write_simulation_log(self, output_dir: str, task_name: str, elapsed: float, success: bool):
        """Write high-level simulation log"""
        from datetime import datetime

        log_file = os.path.join(output_dir, 'simulation.log')

        status = "SUCCESS" if success else "FAILED"
        components_info = []
        for comp in self.config.components:
            comp_info = f"  - {comp.name}: {comp.type.value if hasattr(comp.type, 'value') else comp.type}, nmol={comp.nmol}"
            if comp.type == ComponentType.IDP:
                comp_info += f", seq={comp.ffasta}"
            elif comp.type == ComponentType.MDP:
                comp_info += f", pdb={comp.fpdb}"
            components_info.append(comp_info)

        log_content = f"""# CondenSimAdapter CG Simulation Log
# ============================

Task: {task_name}
Force Field: CALVADOS
Date: {datetime.now().isoformat()}

Status: {status}
Duration: {elapsed:.2f} seconds

System Configuration:
  Box: {self.config.box} nm
  Temperature: {self.config.temperature} K
  Ionic Strength: {self.config.ionic} M
  Topology: {self.config.topol.value if hasattr(self.config.topol, 'value') else self.config.topol}

Components ({len(self.config.components)}):
{chr(10).join(components_info)}

Output Files:
  - final.pdb: Final structure
  - trajectory.xtc: Simulation trajectory
  - raw/: Native simulation output files
"""
        with open(log_file, 'w') as f:
            f.write(log_content)
        print(f"  📝 simulation.log")
    
    def run_hps(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """
        Run HPS-Urry simulation

        Workflow:
        1. Build initial CG structure using CALVADOS force field (pre-equilibration)
        2. Use HPSParser to parse each component (MDP supports domains)
        3. Build HPSModel and add force field
        4. Run OpenMM simulation

        Args:
            gpu_id: GPU device ID
            **kwargs: Additional parameters
                - preequil_steps: Pre-equilibration steps
                - preequil_mapping: Pre-equilibration mapping method ('ca' or 'com')
                - preequil_k_restraint: Pre-equilibration restraint force constant
                - preequil_use_com: Pre-equilibration whether to use COM restraint
                - preequil_platform: Pre-equilibration platform
                - platform: Simulation platform (CPU/CUDA)

        Returns:
            SimulationResult
        """
        self._ensure_setup()
        self._ensure_not_running()

        # Pre-equilibration parameters
        preequil_steps = kwargs.get('preequil_steps', 100000)
        preequil_mapping = kwargs.get('preequil_mapping', 'ca')
        preequil_k_restraint = kwargs.get('preequil_k_restraint', 10000.0)
        preequil_use_com = kwargs.get('preequil_use_com', True)
        preequil_platform = kwargs.get('preequil_platform', self.config.simulation.platform)
        
        # Check CUDA availability, fall back to CPU if not available
        if preequil_platform == ComputePlatform.CUDA:
            try:
                from openmm import Platform
                Platform.getPlatformByName('CUDA')
            except:
                print(f"  ⚠️  CUDA not available, pre-equilibration will use CPU")
                preequil_platform = ComputePlatform.CPU

        # Pre-equilibration (Using CALVADOS to Build Initial Structure)
        preequil_pdb = self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=preequil_steps,
            mapping=preequil_mapping,
            k_restraint=preequil_k_restraint,
            use_com=preequil_use_com,
            platform=preequil_platform,
        )
        if preequil_pdb:
            print(f"  [HPS] Using pre-equilibration structure: {preequil_pdb}")

        self.is_running = True
        result = SimulationResult()

        try:
            print(f"\n[HPS-Urry] Running simulation...")
            print(f"  GPU ID: {gpu_id}")

            # ===== 1. Import HPS related modules =====
            try:
                from CondenSimAdapter.extern.ms2_openabc.forcefields.parsers.hps_parser import HPSParser
                from CondenSimAdapter.extern.ms2_openabc.forcefields import HPSModel
                from CondenSimAdapter.extern.ms2_openabc.lib import _kcal_to_kj
            except ImportError as e:
                raise ImportError(f"ms2_openabc module not available: {e}")

            # ===== 2. Create HPSParser for each component =====
            parsers = []
            print("\n  Building HPSParser...")
            
            for comp in self.config.components:
                if comp.type == ComponentType.MDP:
                    # MDP: Use provided PDB + domain definition
                    if not comp.fpdb:
                        raise ValueError(f"MDP component '{comp.name}' requires fpdb file")
                    
                    # Import tool function for converting atomistic structure to CA
                    from CondenSimAdapter.extern.ms2_openabc.utils.helper_functions import atomistic_pdb_to_ca_pdb
                    
                    # Create temporary CA-only PDB file (if original PDB is not CA-only)
                    ca_pdb_path = comp.fpdb
                    temp_ca_pdb = None
                    
                    # Check if original PDB is CA-only
                    from CondenSimAdapter.extern.ms2_openabc.utils import parse_pdb
                    original_atoms = parse_pdb(comp.fpdb)
                    if not (original_atoms['name'].eq('CA').all()):
                        # Original PDB is all-atom, need to convert to CA-only
                        import os
                        temp_ca_pdb = os.path.join(self.output_dir, f'_temp_{comp.name}_ca.pdb')
                        atomistic_pdb_to_ca_pdb(comp.fpdb, temp_ca_pdb)
                        ca_pdb_path = temp_ca_pdb
                    
                    # Handle fdomains: could be file path or inline YAML content
                    fdomains_path = comp.fdomains
                    temp_domains_file = None
                    
                    if comp.fdomains:
                        # Check if it's inline YAML content (contains newlines and YAML format)
                        if '\n' in comp.fdomains and ('[' in comp.fdomains or '-' in comp.fdomains):
                            # Inline YAML content, need to write to temp file
                            import os
                            import tempfile
                            temp_domains_file = os.path.join(self.output_dir, f'_temp_{comp.name}_domains.yaml')
                            with open(temp_domains_file, 'w') as f:
                                f.write(comp.fdomains)
                            fdomains_path = temp_domains_file
                        elif not os.path.isfile(comp.fdomains):
                            # Neither file nor inline content, might be string format list
                            try:
                                import ast
                                # Try to parse as Python list
                                domains_list = ast.literal_eval(comp.fdomains)
                                import os
                                temp_domains_file = os.path.join(self.output_dir, f'_temp_{comp.name}_domains.yaml')
                                # Convert to YAML format
                                import yaml
                                yaml_content = {comp.name: domains_list}
                                with open(temp_domains_file, 'w') as f:
                                    yaml.dump(yaml_content, f)
                                fdomains_path = temp_domains_file
                            except:
                                # If parsing fails, set to None for HPSParser to handle
                                fdomains_path = None
                    
                    parser = HPSParser(
                        ca_pdb=ca_pdb_path,
                        fdomains=fdomains_path
                    )
                    print(f"    {comp.name}: MDP, {len(parser.atoms)} CA atoms")
                    
                    if parser.enm_pairs:
                        print(f"      → {len(parser.enm_pairs)} ENM pairs")
                    
                    # Clean up temporary files
                    if temp_ca_pdb and os.path.exists(temp_ca_pdb):
                        os.remove(temp_ca_pdb)
                    if temp_domains_file and os.path.exists(temp_domains_file):
                        os.remove(temp_domains_file)
                    
                    parsers.append((comp, parser))
                    
                elif comp.type == ComponentType.IDP:
                    # IDP: Build straight CA chain from FASTA sequence
                    if not comp.ffasta:
                        raise ValueError(f"IDP component '{comp.name}' requires ffasta file")
                    
                    # Read FASTA sequence
                    with open(comp.ffasta, 'r') as f:
                        fasta_content = f.read()
                    
                    # Parse FASTA (extract sequence matching component name)
                    lines = fasta_content.strip().split('\n')
                    sequence = None
                    current_seq_lines = []
                    in_target_sequence = False
                    
                    for line in lines:
                        if line.startswith('>'):
                            # If previous segment was our target sequence, save it
                            if in_target_sequence:
                                sequence = ''.join(current_seq_lines)
                                break
                            # Check if this line is the sequence we're looking for
                            in_target_sequence = (comp.name in line.replace('>', '').strip())
                            current_seq_lines = []
                        elif in_target_sequence:
                            current_seq_lines.append(line.strip())
                    
                    # Handle the last sequence
                    if sequence is None and in_target_sequence:
                        sequence = ''.join(current_seq_lines)
                    
                    if sequence is None:
                        raise ValueError(f"IDP component '{comp.name}': Sequence not found in FASTA")
                    
                    # Use build_straight_CA_chain to build straight CA chain from sequence
                    from CondenSimAdapter.extern.ms2_openabc.utils.helper_functions import build_straight_CA_chain, write_pdb
                    import os
                    
                    # Create temporary PDB file
                    temp_pdb = os.path.join(self.output_dir, f'_temp_idp_{comp.name}_ca.pdb')
                    ca_atoms = build_straight_CA_chain(sequence, r0=0.38)
                    write_pdb(ca_atoms, temp_pdb)
                    
                    # Create HPSParser using temporary PDB
                    parser = HPSParser(
                        ca_pdb=temp_pdb,
                        fdomains=None
                    )
                    
                    print(f"    {comp.name}: IDP, {len(sequence)} residues (built from FASTA sequence)")
                    parsers.append((comp, parser))
            
            if not parsers:
                raise ValueError("No valid components")

            # ===== 3. Build HPSModel =====
            print("\n  Building HPSModel...")
            model = HPSModel()
            
            # Set periodic boundary conditions (all topology types are periodic)
            is_periodic = True
            model.use_pbc = is_periodic
            
            # Add all molecules (considering each component's nmol)
            for comp, parser in parsers:
                n_before = len(model.atoms) if model.atoms is not None else 0
                # Add nmol copies for each component
                for _ in range(comp.nmol):
                    model.append_mol(parser)
                n_after = len(model.atoms)
                print(f"    {comp.name}: added {comp.nmol} copies, {n_after - n_before} total atoms")

            # ===== 4. Create Topology and System =====
            print("\n  Creating Topology and System...")
            
            # Create HPS output subdirectory (for temp files and output files)
            hps_output_dir = os.path.join(self.output_dir, 'HPS')
            os.makedirs(hps_output_dir, exist_ok=True)
            
            # Create temp PDB from model.atoms, then read to create topology (in HPS subdirectory)
            temp_pdb = os.path.join(hps_output_dir, '_temp_hps_model.pdb')
            model.atoms_to_pdb(temp_pdb, reset_serial=True)
            topology = PDBFile(temp_pdb).topology
            
            # Create OpenMM System
            box_size = self.config.box
            is_periodic = True  # All topology types are periodic
            model.create_system(
                top=topology,
                use_pbc=is_periodic,
                box_a=box_size[0],
                box_b=box_size[1],
                box_c=box_size[2]
            )
            print(f"    ✓ System created: {model.system.getNumParticles()} particles")
            
            # ===== 5. Add Force Field =====
            print("\n  Adding Force Field...")
            
            # Protein bonds (Harmonic)
            print("    - Protein bonds (harmonic)")
            model.add_protein_bonds(force_group=1)
            
            # Non-bonded contacts (Ashbaugh-Hatch with Urry scale)
            print("    - Contacts (Ashbaugh-Hatch, Urry scale)")
            model.add_contacts(
                hydropathy_scale='Urry',
                epsilon=0.2 * _kcal_to_kj,  # Convert kcal to kJ
                mu=1,
                delta=0.08,
                force_group=2
            )
            
            # Electrostatics (Debye-Hückel)
            print("    - Electrostatics (Debye-Hückel)")
            model.add_dh_elec(
                ldby=1 * unit.nanometer,
                dielectric_water=80.0,
                cutoff=3.5 * unit.nanometer,
                force_group=3
            )
            
            # Elastic network (for MDP)
            has_enm = any(p[1].enm_pairs for p in parsers)
            if has_enm:
                print("    - Elastic network (for folded domains)")
                model.add_elastic_network(
                    force_constant=700.0 * unit.kilojoule_per_mole / unit.nanometer ** 2,
                    force_group=4
                )

            # ===== 6. Create OpenMM Simulation =====
            print("\n  Creating OpenMM Simulation...")
            
            system = model.system
            
            if _is_droplet_topology(self.config.topol):
                radius, center = _get_droplet_params(self.config)
                system = add_droplet_force(
                    system=system,
                    radius=radius,
                    center=center,
                    k=DROPLET_FORCE_K,
                    stride=DROPLET_FORCE_STRIDE,
                )
                print(f"  Droplet confinement enabled (k={DROPLET_FORCE_K}, stride={DROPLET_FORCE_STRIDE})")
            
            # Read initial coordinates from CALVADOS pre-equilibration structure (consistent with COCOMO mode)
            if preequil_pdb and os.path.exists(preequil_pdb):
                print(f"  Reading initial coordinates from pre-equilibration structure: {preequil_pdb}")
                pdb = PDBFile(preequil_pdb)
                positions = pdb.getPositions(asNumpy=True)
                
                # Verify coordinate count matches
                if len(positions) != len(model.atoms):
                    print(f"  ⚠️  Warning: Pre-equilibration structure atom count ({len(positions)}) does not match model atom count ({len(model.atoms)})")
                    print(f"  Reading coordinates from temporary PDB")
                    positions = PDBFile(temp_pdb).getPositions(asNumpy=True)
                else:
                    print(f"  ✓ Coordinate count matches: {len(positions)} atoms")
            else:
                print(f"  ⚠️  Pre-equilibration structure not found, reading coordinates from temporary PDB")
                positions = PDBFile(temp_pdb).getPositions(asNumpy=True)

            # Select platform
            config_platform = self.config.simulation.platform
            platform_name = config_platform.value if hasattr(config_platform, 'value') else str(config_platform)
            
            # Try GPU, fall back to CPU if unavailable
            if gpu_id >= 0:
                for pname in [platform_name, 'CUDA', 'OpenCL']:
                    try:
                        platform = Platform.getPlatformByName(pname)
                        properties = {'DeviceIndex': str(gpu_id)}
                        print(f"    Platform: {pname}")
                        break
                    except:
                        continue
                else:
                    platform = Platform.getPlatformByName('CPU')
                    properties = {}
                    print(f"    Platform: CPU (GPU unavailable)")
            else:
                platform = Platform.getPlatformByName('CPU')
                properties = {}
                print(f"    Platform: CPU")

            # Create integrator
            temperature = self.config.temperature
            integrator = LangevinMiddleIntegrator(
                temperature * kelvin,
                0.01 / picoseconds,
                0.01 * picoseconds
            )

            simulation = Simulation(
                topology,
                system,
                integrator=integrator,
                platform=platform,
                platformProperties=properties
            )

            # Set initial positions
            simulation.context.setPositions(positions)

            # Set box vectors
            if is_periodic:
                box_size = self.config.box
                box_vecs = [
                    mm.Vec3(x=box_size[0], y=0.0, z=0.0),
                    mm.Vec3(x=0.0, y=box_size[1], z=0.0),
                    mm.Vec3(x=0.0, y=0.0, z=box_size[2])
                ] * unit.nanometer
                simulation.context.setPeriodicBoxVectors(*box_vecs)

            # Energy minimization
            print("  Energy minimization...")
            simulation.minimizeEnergy()
            print("  Energy minimization completed")
            simulation.context.setVelocitiesToTemperature(temperature * kelvin)

            # ===== 6. Run simulation =====
            print(f"\n  Running HPS simulation: {self.config.simulation.steps} steps...")
            
            # Add reporters (using HPS output directory created earlier)
            output_dir = hps_output_dir
            
            wfreq = self.config.simulation.wfreq
            traj_file = os.path.join(output_dir, 'trajectory.xtc')
            log_file = os.path.join(output_dir, 'simulation.log')
            
            simulation.reporters.append(XTCReporter(traj_file, wfreq))
            simulation.reporters.append(StateDataReporter(
                log_file,
                wfreq,
                step=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                volume=True,
            ))

            # Progress bar
            from tqdm import tqdm
            total_steps = self.config.simulation.steps
            batch_size = 1000
            
            for _ in tqdm(range(0, total_steps, batch_size), desc="HPS"):
                simulation.step(min(batch_size, total_steps - simulation.currentStep))

            print("  Simulation completed!")

            # ===== 7. Save results =====
            print("\n  Saving results...")
            
            # Final structure
            state_final = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=True,
                getEnergy=True,
                enforcePeriodicBox=True
            )
            positions_final = state_final.getPositions()
            
            # Get box vectors and set to topology (so PDB will contain PBC info)
            box_vectors = state_final.getPeriodicBoxVectors()
            simulation.topology.setPeriodicBoxVectors(box_vectors)
            
            final_pdb = os.path.join(output_dir, 'final.pdb')
            final_pdb_hps_format = os.path.join(output_dir, 'final_hps_format.pdb')
            with open(final_pdb_hps_format, 'w') as f:
                PDBFile.writeFile(topology, positions_final, f, keepIds=True)
            print(f"    - final_hps_format.pdb")

            # Post-processing: Convert HPS format PDB to calvados format (chain/resSeq unified)
            print(f"\n  Post-processing PDB format conversion (HPS -> Calvados)...")
            try:
                from .backmap import standardize_pdb_with_calvados
                standardize_pdb_with_calvados(
                    pdb_path=final_pdb_hps_format,
                    config=self.config,
                    output_pdb=final_pdb
                )
                print(f"  ✓ PDB format conversion completed: {final_pdb}")
            except Exception as e:
                print(f"  ⚠️  PDB format conversion failed: {e}")
                print(f"  Keeping original HPS format PDB: {final_pdb_hps_format}")
                shutil.copy2(final_pdb_hps_format, final_pdb)
            
            # Copy to root directory (output_dir is HPS/, need to copy to TDP43_CTD_CG/)
            system_name = self.config.system_name
            cg_root_dir = os.path.dirname(output_dir)  # output_dir is .../TDP43_CTD_CG/HPS, so dirname is TDP43_CTD_CG
            final_pdb_root = os.path.join(cg_root_dir, 'final.pdb')
            os.makedirs(cg_root_dir, exist_ok=True)
            shutil.copy2(final_pdb, final_pdb_root)
            print(f"  Copied final structure to: {final_pdb_root}")

            # Save system XML
            system_xml = os.path.join(output_dir, 'system.xml')
            with open(system_xml, 'w') as f:
                f.write(XmlSerializer.serialize(system))
            print(f"    - system.xml")

            # Save checkpoint
            checkpoint_file = os.path.join(output_dir, 'restart.chk')
            simulation.saveCheckpoint(checkpoint_file)
            print(f"    - restart.chk")

            result.success = True
            result.trajectory = traj_file
            result.structure = final_pdb
            result.output_dir = output_dir
            
            print(f"\n  HPS-Urry output directory: {output_dir}")

        except ImportError as e:
            result.success = False
            result.errors.append(f"ms2_openabc module not available: {e}")
            print(f"  ✗ HPS-Urry simulation failed: ms2_openabc module not available")
            import traceback
            traceback.print_exc()
        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            print(f"  ✗ HPS-Urry simulation failed: {e}")
            import traceback
            traceback.print_exc()

        finally:
            self.is_running = False

        self._result = result
        return result

    def run_moff(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """
        Run MOFF simulation

        Automatically perform pre-equilibration (Using CALVADOS to Build Initial Structure):
        1. Build initial CG structure using CALVADOS force field
        2. Run short simulation for MDP protein
        3. Then switch to MOFF force field for formal simulation

        Uses MOFF force field from OpenABC package.

        Pre-equilibration parameters (hardcoded, can be overridden via kwargs):
        - steps: Pre-equilibration steps (default 100000)
        - mapping: Mapping method 'ca' or 'com' (default 'ca')
        - k_restraint: Restraint force constant (default 10000.0)
        - use_com: Whether to use COM restraint (default True)
        - platform: Computing platform (default from config, CUDA)

        Args:
            gpu_id: GPU device ID
            **kwargs: Additional parameters
                - salt_conc: Salt concentration (default 150 mM)
                - preequil_steps: Pre-equilibration steps
                - preequil_mapping: Pre-equilibration mapping method ('ca' or 'com')
                - preequil_k_restraint: Pre-equilibration restraint force constant
                - preequil_use_com: Pre-equilibration whether to use COM restraint
                - preequil_platform: Pre-equilibration platform (CUDA or CPU)

        Returns:
            SimulationResult
        """
        self._ensure_setup()
        self._ensure_not_running()

        # Pre-equilibration parameters (hardcoded, can be overridden)
        preequil_steps = kwargs.get('preequil_steps', 100000)
        preequil_mapping = kwargs.get('preequil_mapping', 'ca')
        preequil_k_restraint = kwargs.get('preequil_k_restraint', 10000.0)
        preequil_use_com = kwargs.get('preequil_use_com', True)
        # Pre-equilibration platform: prefer value from kwargs, otherwise use config (default CUDA)
        preequil_platform = kwargs.get('preequil_platform', self.config.simulation.platform)

        # Pre-equilibration (Using CALVADOS to Build Initial Structure)
        preequil_pdb = self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=preequil_steps,
            mapping=preequil_mapping,
            k_restraint=preequil_k_restraint,
            use_com=preequil_use_com,
            platform=preequil_platform,
        )
        if preequil_pdb:
            print(f"  [MOFF] Using pre-equilibration structure: {preequil_pdb}")

        self.is_running = True
        result = SimulationResult()
        result.output_dir = self.output_dir

        try:
            print(f"\n[MOFF] Running simulation...")
            print(f"  GPU ID: {gpu_id}")

            # TODO: Implement MOFF runner
            # Need to use openabc.forcefields.MOFFModel

            result.success = True
            print(f"  ✓ MOFF simulation completed (placeholder)")

        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            print(f"  ✗ MOFF simulation failed: {e}")

        finally:
            self.is_running = False

        self._result = result
        return result
    
    def run_cocomo(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """
        Run COCOMO simulation

        Automatically perform pre-equilibration (Using CALVADOS to Build Initial Structure):
        1. Build initial CG structure using CALVADOS force field
        2. Run short simulation for MDP protein
        3. Then switch to COCOMO force field for formal simulation

        COCOMO uses COCOMO2 force field with SASA correction enabled by default.
        Ionic strength is fixed at 0.1 M (force field default value).

        Pre-equilibration parameters (hardcoded, can be overridden via kwargs):
        - steps: Pre-equilibration steps (default 100000)
        - mapping: Mapping method 'ca' or 'com' (default 'ca')
        - k_restraint: Restraint force constant (default 10000.0)
        - use_com: Whether to use COM restraint (default True)
        - platform: Computing platform (default from config, CUDA)

        Args:
            gpu_id: GPU device ID
            **kwargs: Additional parameters
                - nstep: Simulation steps (override config value)
                - wfreq: Write frequency
                - tstep: Time step (ps, default 0.01)
                - gamma: Friction coefficient (1/ps, default 0.01)
                - preequil_steps: Pre-equilibration steps
                - preequil_mapping: Pre-equilibration mapping method ('ca' or 'com')
                - preequil_k_restraint: Pre-equilibration restraint force constant
                - preequil_use_com: Pre-equilibration whether to use COM restraint
                - preequil_platform: Pre-equilibration platform (CUDA or CPU)

        Returns:
            SimulationResult
        """
        import warnings as py_warnings

        self._ensure_setup()
        self._ensure_not_running()

        # Warning for ionic strength
        if abs(self.config.ionic - 0.1) > 0.001:
            py_warnings.warn(
                "COCOMO uses default ionic strength of 0.1 M. "
                f"Config value {self.config.ionic} M will be ignored.",
                UserWarning
            )

        # Pre-equilibration parameters (hardcoded, can be overridden)
        preequil_steps = kwargs.get('preequil_steps', 100000)
        preequil_mapping = kwargs.get('preequil_mapping', 'ca')
        preequil_k_restraint = kwargs.get('preequil_k_restraint', 10000.0)
        preequil_use_com = kwargs.get('preequil_use_com', True)
        # Pre-equilibration platform: prefer value from kwargs, otherwise use config (default CUDA)
        preequil_platform = kwargs.get('preequil_platform', self.config.simulation.platform)

        # Pre-equilibration (Using CALVADOS to Build Initial Structure)
        preequil_pdb = self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=preequil_steps,
            mapping=preequil_mapping,
            k_restraint=preequil_k_restraint,
            use_com=preequil_use_com,
            platform=preequil_platform,
        )

        # Prepare output directory
        dirs = self._prepare_cocomo_output()
        output_dir = dirs['output_dir']

        self.is_running = True
        result = SimulationResult()

        try:
            print(f"\n[COCOMO] Running COCOMO2 simulation...")
            print(f"  Pre-equilibration structure: {preequil_pdb}")

            # Run simulation using new COCOMO class
            from .cocomo2_creator import COCOMO

            # Build topology_info dictionary
            topology_info = {
                'global_sequence': self.get_global_sequence(),
                'chain_ids': self.get_chain_ids(),
                'folded_domains': self.get_folded_domains(),
                'component_names': self._build_component_names(),
                'local_residue_indices': list(range(1, len(self.get_global_sequence()) + 1)),
            }

            # Read pre-equilibration structure coordinates
            from openmm.app import PDBFile
            pdb = PDBFile(preequil_pdb)
            positions = pdb.getPositions(asNumpy=True)

            # Get box size (use config box parameter, in nm)
            box_size = self.config.box

            # Get SASA values
            sasa_values = self._get_sasa_values()

            if sasa_values is not None:
                topology_info['sasa_values'] = sasa_values

            # Create COCOMO system
            cocomo = COCOMO(
                box_size=box_size,
                topology_info=topology_info,
                positions=positions,
                surf=0.7,
                resources='CUDA' if gpu_id >= 0 else 'CPU'
            )

            # Set _topology_info for ENM calculation
            # Use SimpleNamespace to create object with attributes
            from types import SimpleNamespace
            cocomo._topology_info = SimpleNamespace(
                global_sequence=topology_info['global_sequence'],
                chain_ids=topology_info['chain_ids'],
                is_folded=topology_info['folded_domains'],
                sasa_values=topology_info.get('sasa_values', [])
            )

            # Create OpenMM system
            system, topology = cocomo.create_system()
            
            if _is_droplet_topology(self.config.topol):
                radius, center = _get_droplet_params(self.config)
                system = add_droplet_force(
                    system=system,
                    radius=radius,
                    center=center,
                    k=DROPLET_FORCE_K,
                    stride=DROPLET_FORCE_STRIDE,
                )
                print(f"  Droplet confinement enabled (k={DROPLET_FORCE_K}, stride={DROPLET_FORCE_STRIDE})")

            # Create Simulation object (prefer platform specified in config)
            config_platform = self.config.simulation.platform
            platform_name = config_platform.value if hasattr(config_platform, 'value') else str(config_platform)
            
            # Platform selection: prefer config-specified platform, support automatic fallback
            if gpu_id >= 0:
                # User wants to use GPU
                try:
                    platform = Platform.getPlatformByName(platform_name)
                    properties = {'DeviceIndex': str(gpu_id)}
                    print(f"  Using {platform_name} platform")
                except Exception:
                    # Try other GPU platforms
                    for fallback in ['CUDA', 'OpenCL']:
                        if fallback != platform_name:
                            try:
                                platform = Platform.getPlatformByName(fallback)
                                properties = {'DeviceIndex': str(gpu_id)}
                                print(f"  {platform_name} unavailable, falling back to {fallback}")
                                break
                            except:
                                continue
                    else:
                        print(f"  GPU unavailable, falling back to CPU")
                        platform = Platform.getPlatformByName('CPU')
                        properties = {}
            else:
                # User wants to use CPU
                platform = Platform.getPlatformByName('CPU')
                properties = {}
                print("  Using CPU platform")

            # Use LangevinIntegrator (consistent with Legacy version)
            # Temperature from config, friction coefficient 0.01/ps, time step 0.01 ps
            temperature = self.config.temperature
            integrator = LangevinIntegrator(
                temperature * kelvin,
                0.01 / picoseconds,
                0.01 * picoseconds
            )

            simulation = Simulation(
                topology,
                system,
                integrator=integrator,
                platform=platform,
                platformProperties=properties
            )

            # Set initial positions
            simulation.context.setPositions(positions)

            # Set initial velocities (consistent with Legacy version)
            simulation.context.setVelocitiesToTemperature(temperature * kelvin)

            # Energy minimization
            print("  Running energy minimization...")
            simulation.minimizeEnergy()
            print("  Energy minimization completed")

            # Run simulation
            wfreq = self.config.simulation.wfreq
            xtc_file = os.path.join(output_dir, 'trajectory.xtc')
            log_file = os.path.join(output_dir, 'simulation.log')

            # Add reporters (use mdtraj's XTCReporter to preserve PBC info)
            simulation.reporters.append(
                XTCReporter(xtc_file, wfreq)
            )
            simulation.reporters.append(
                StateDataReporter(
                    log_file,
                    wfreq,
                    step=True,
                    potentialEnergy=True,
                    kineticEnergy=True,
                    totalEnergy=True,
                    temperature=True,
                    volume=True,
                )
            )

            # Run simulation (use tqdm to show progress)
            total_steps = self.config.simulation.steps
            print(f"  Starting simulation: {total_steps} steps")

            from tqdm import tqdm
            n_batches = 10
            batch_size = total_steps // n_batches

            for _ in tqdm(range(n_batches), desc="COCOMO"):
                simulation.step(batch_size)
                simulation.saveCheckpoint(os.path.join(output_dir, 'restart.chk'))

            # Handle remaining steps
            remaining = total_steps % n_batches
            if remaining > 0:
                simulation.step(remaining)

            print(f"  Simulation completed!")

            # Get final state (contains PBC info)
            state_final = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=True,
                getEnergy=True,
                enforcePeriodicBox=True  # Force periodic boundary conditions
            )

            # Get final positions and box vectors
            positions_final = state_final.getPositions()
            box_vectors = state_final.getPeriodicBoxVectors()

            # Set box vectors on topology (so PDB will contain PBC info)
            simulation.topology.setPeriodicBoxVectors(box_vectors)

            # Save final structure as PDB format (contains PBC)
            from openmm.app import PDBFile
            final_pdb = os.path.join(output_dir, 'final.pdb')
            with open(final_pdb, 'w') as f:
                PDBFile.writeFile(
                    simulation.topology,
                    positions_final,
                    f
                )
            print(f"  Saved final structure: {final_pdb}")

            # final.pdb is already in root directory (output_dir = self.output_dir = {system_name}_CG), no need to copy

            # Save system XML
            system_xml = os.path.join(output_dir, 'system.xml')
            with open(system_xml, 'w') as f:
                f.write(XmlSerializer.serialize(system))
            print(f"  Saved system XML: {system_xml}")

            result.success = True
            result.trajectory = xtc_file
            result.structure = final_pdb
            result.output_dir = output_dir

            print(f"  COCOMO output directory: {output_dir}")

        except Exception as e:
            result = SimulationResult()
            result.success = False
            result.errors.append(str(e))
            result.output_dir = output_dir
            print(f"  ✗ COCOMO simulation failed: {e}")
            import traceback
            traceback.print_exc()

        finally:
            self.is_running = False

        self._result = result
        return result
    
    def _prepare_cocomo_output(self) -> Dict[str, str]:
        """
        Prepare COCOMO output directory

        Output directory structure: {output_dir}/ (CLI already sets self.output_dir = {system_name}_CG)
        Use self.output_dir directly to avoid nested directories
        """
        # CLI already sets self.output_dir = {system_name}_CG, use directly
        output_dir = self.output_dir

        # Backup old results (only when directory contains previous simulation results)
        import shutil
        from datetime import datetime

        # Check if directory exists and contains previous simulation results
        should_backup = False
        if os.path.exists(output_dir):
            # Check for simulation result files (excluding pre-equilibration files)
            result_files = ['final.pdb', 'trajectory.xtc', 'simulation.log', 'system.xml']
            has_results = any(os.path.exists(os.path.join(output_dir, f)) for f in result_files)
            
            # Only backup if simulation results exist (pre-equilibration files don't count)
            if has_results:
                should_backup = True
        
        if should_backup:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{output_dir}_backup_{timestamp}"
            shutil.move(output_dir, backup_dir)
            print(f"  📁 Backing up old results to: {backup_dir}")
            # Recreate output directory
            os.makedirs(output_dir, exist_ok=True)
        else:
            # Directory doesn't exist or only contains pre-equilibration files, create directly if needed
            os.makedirs(output_dir, exist_ok=True)

        return {
            'output_dir': output_dir,
            'task_name': 'COCOMO',
        }

    def run_mpipi_recharged(self, gpu_id: int = 0, return_system: bool = False, **kwargs) -> SimulationResult:
        """
        Run Mpipi-Recharged simulation

        Use the same logic as test_100_molecules.py:
        1. Create biomolecule objects (MDP/IDP)
        2. Use gmx insert-molecules to place molecules
        3. Build Mpipi-Recharged system
        4. Energy minimization
        5. Run MD simulation

        Note: Optionally use CALVADOS pre-equilibration to generate initial structure.

        Args:
            gpu_id: GPU device ID
            return_system: If True, return (system, topology, positions) instead of running simulation (default False)
            **kwargs: Additional parameters
                - use_gmx_insert: Whether to use gmx insert-molecules to place molecules (default True)
                - gmx_radius: Minimum distance for gmx insert-molecules (nm, default 0.35)
                - verbose: Whether to output detailed information (default True)

        Returns:
            SimulationResult (default)
            Or (system, topology, positions) when return_system=True
        """
        from openmm.app import PDBFile
        from openmm import Platform, LangevinMiddleIntegrator, Vec3, LocalEnergyMinimizer
        import openmm as mm
        from openmm import unit
        from openmm.unit import kelvin, picoseconds, nanometer

        self._ensure_setup()
        self._ensure_not_running()

        # Pre-equilibration parameters (optional)
        preequil_steps = kwargs.get('preequil_steps', 100000)
        preequil_mapping = kwargs.get('preequil_mapping', 'ca')
        preequil_k_restraint = kwargs.get('preequil_k_restraint', 10000.0)
        preequil_use_com = kwargs.get('preequil_use_com', True)
        preequil_platform = kwargs.get('preequil_platform', self.config.simulation.platform)

        if preequil_platform == ComputePlatform.CUDA:
            try:
                from openmm import Platform
                Platform.getPlatformByName('CUDA')
            except Exception:
                print(f"  ⚠️  CUDA unavailable, pre-equilibration will use CPU")
                preequil_platform = ComputePlatform.CPU

        preequil_pdb = self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=preequil_steps,
            mapping=preequil_mapping,
            k_restraint=preequil_k_restraint,
            use_com=preequil_use_com,
            platform=preequil_platform,
        )
        if preequil_pdb:
            print(f"  [Mpipi] Using pre-equilibrated structure: {preequil_pdb}")

        # Output directory
        dirs = self._prepare_mpipi_output()
        output_dir = dirs['output_dir']

        self.is_running = True
        result = SimulationResult()

        try:
            print(f"\n[Mpipi-Recharged] Running simulation...")

            # Import biomolecule class and new functions from ms2_OpenMpipi
            from CondenSimAdapter.extern.ms2_OpenMpipi import MDP, IDP, build_mpipi_recharged_system_from_chains

            # Build biomolecule objects from config.components
            print("\n  Building biomolecule objects...")
            chain_objects = []
            for comp in self.config.components:
                # Get sequence
                if comp.type == ComponentType.IDP:
                    # IDP reads sequence from FASTA
                    if comp.ffasta:
                        with open(comp.ffasta, 'r') as f:
                            fasta_content = f.read()
                        # Parse FASTA (only extract sequences matching component names)
                        lines = fasta_content.strip().split('\n')
                        sequence = None
                        current_seq_lines = []
                        in_target_sequence = False
                        
                        for line in lines:
                            if line.startswith('>'):
                                # If the previous entry is the target sequence, save it
                                if in_target_sequence:
                                    sequence = ''.join(current_seq_lines)
                                    break
                                # Check whether this line is the target sequence
                                in_target_sequence = (comp.name in line.replace('>', '').strip())
                                current_seq_lines = []
                            elif in_target_sequence:
                                current_seq_lines.append(line.strip())
                        
                        # Handle the last sequence
                        if sequence is None and in_target_sequence:
                            sequence = ''.join(current_seq_lines)
                        
                        if sequence is None:
                            print(f"    ⚠️  {comp.name}: sequence not found in FASTA, skipping")
                            continue
                    else:
                        continue
                    
                    # Create IDP object
                    idp = IDP(comp.name, sequence)
                    chain_objects.append((comp, idp))
                    print(f"    {comp.name}: IDP, {len(sequence)} residues")

                elif comp.type == ComponentType.MDP:
                    # MDP reads structure and sequence from fpdb
                    if not comp.fpdb:
                        continue
                    
                    # Parse fdomains to get folded-domain indices
                    domains = self._parse_fdomains(comp.fdomains)
                    
                    # Read sequence from PDB (MDP always uses the PDB sequence)
                    # Note: PDB residue names are three-letter codes; convert to one-letter codes
                    pdb_temp = PDBFile(comp.fpdb)
                    three_to_one = {
                        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                    }
                    sequence = []
                    for res in pdb_temp.topology.residues():
                        if res.name in three_to_one:
                            sequence.append(three_to_one[res.name])
                    sequence = ''.join(sequence)
                    
                    # Create MDP object
                    mdp = MDP(comp.name, sequence, domains, comp.fpdb)
                    chain_objects.append((comp, mdp))
                    print(f"    {comp.name}: MDP, {len(sequence)} residues, {len(domains)} domains")

            if not chain_objects:
                raise ValueError("No valid biomolecules could be constructed for Mpipi simulation.")

            # Build chain_info dict
            print("\n  Building chain_info dict...")
            chain_info = {}
            for comp, biomol in chain_objects:
                chain_info[biomol] = comp.nmol
                print(f"    {comp.name}: {comp.nmol} copies")

            # Call build_mpipi_recharged_system_from_chains
            # This function will: 1) relax each monomer 2) build model 3) add force field
            csx = self.config.ionic * 1000  # M -> mM
            is_periodic = True  # All topology types are periodic
            box_size = self.config.box  # Use box size from config
            
            print("\n  Building Mpipi-Recharged system (relaxation + model building + force field)...")
            # Use only coordinates from pre-equilibrated structure; no molecular placement
            if not preequil_pdb or not os.path.exists(preequil_pdb):
                raise FileNotFoundError(
                    f"Pre-equilibrated structure not found: {preequil_pdb}. Please generate preequil_final.pdb first"
                )
            use_gmx_insert = False
            use_grid_placement = False
            gmx_radius = kwargs.get('gmx_radius', 0.35)  # Keep parameter for API compatibility
            verbose = kwargs.get('verbose', True)
            print("  Using pre-equilibrated coordinates; skipping placement")
            
            system, model = build_mpipi_recharged_system_from_chains(
                chain_info=chain_info,
                box_size=box_size,
                topol=self.config.topol.value,
                T=self.config.temperature * unit.kelvin,
                csx=csx,
                CM_remover=True,
                periodic=is_periodic,
                use_gmx_insert=use_gmx_insert,  # Use gmx insert-molecules
                use_grid_placement=use_grid_placement,
                gmx_radius=gmx_radius,  # Minimum distance
                verbose=verbose
            )
            print(f"  System built: {system.getNumParticles()} particles, {system.getNumForces()} forces")

            # Get topology and positions from model
            topology = model.topology
            positions = model.positions

            # Use pre-equilibrated structure coordinates
            print(f"  Reading initial coordinates from pre-equilibrated structure: {preequil_pdb}")
            pdb = PDBFile(preequil_pdb)
            pre_positions = pdb.getPositions(asNumpy=True)
            model_atom_count = sum(1 for _ in model.topology.atoms())
            if len(pre_positions) != model_atom_count:
                raise ValueError(
                    f"Pre-equilibrated atom count ({len(pre_positions)}) does not match model atom count ({model_atom_count})"
                )
            print(f"  ✓ Coordinate count matches: {len(pre_positions)} atoms")
            positions = pre_positions
            
            if _is_droplet_topology(self.config.topol):
                radius, center = _get_droplet_params(self.config)
                system = add_droplet_force(
                    system=system,
                    radius=radius,
                    center=center,
                    k=DROPLET_FORCE_K,
                    stride=DROPLET_FORCE_STRIDE,
                )
                print(f"  Droplet confinement enabled (k={DROPLET_FORCE_K}, stride={DROPLET_FORCE_STRIDE})")

            # If only system is needed, skip simulation
            if return_system:
                print(f"\n  [return_system=True] Returning system, skipping simulation")
                self.is_running = False
                return system, topology, positions

            # Create Simulation object (prefer platform specified in config)
            config_platform = self.config.simulation.platform
            platform_name = config_platform.value if hasattr(config_platform, 'value') else str(config_platform)
            
            # Platform selection: prefer config platform, support fallback
            if gpu_id >= 0:
                # User wants GPU
                try:
                    platform = Platform.getPlatformByName(platform_name)
                    properties = {'DeviceIndex': str(gpu_id)}
                    print(f"  Using {platform_name} platform (GPU {gpu_id})")
                except Exception:
                    # Try other GPU platforms
                    for fallback in ['CUDA', 'OpenCL']:
                        if fallback != platform_name:
                            try:
                                platform = Platform.getPlatformByName(fallback)
                                properties = {'DeviceIndex': str(gpu_id)}
                                print(f"  {platform_name} unavailable, falling back to {fallback}")
                                break
                            except:
                                continue
                    else:
                        print(f"  GPU unavailable, falling back to CPU")
                        platform = Platform.getPlatformByName('CPU')
                        properties = {}
            else:
                # User wants CPU
                platform = Platform.getPlatformByName('CPU')
                properties = {}
                print("  Using CPU platform")

            # Use LangevinMiddleIntegrator (consistent with test_100_molecules.py)
            temperature = self.config.temperature
            integrator = LangevinMiddleIntegrator(
                temperature * kelvin,
                0.01 / picoseconds,
                0.01 * picoseconds
            )

            simulation = mm.app.Simulation(
                topology,
                system,
                integrator=integrator,
                platform=platform,
                platformProperties=properties
            )

            # Set initial coordinates: use pre-equilibrated structure coordinates
            print("  Using pre-equilibrated structure coordinates")
            simulation.context.setPositions(positions)

            # If the system has PBC, set box vectors
            if is_periodic:
                box_size = self.config.box
                box_vecs = [
                    mm.Vec3(x=box_size[0], y=0.0, z=0.0),
                    mm.Vec3(x=0.0, y=box_size[1], z=0.0),
                    mm.Vec3(x=0.0, y=0.0, z=box_size[2])
                ] * unit.nanometer
                simulation.context.setPeriodicBoxVectors(*box_vecs)

            # Save structure before minimization (for debugging)
            debug_pdb_before_min = os.path.join(output_dir, 'before_minimization.pdb')
            with open(debug_pdb_before_min, 'w') as f:
                PDBFile.writeFile(
                    simulation.topology,
                    positions,
                    f
                )
            print(f"  Saving structure before minimization: {debug_pdb_before_min}")

            # Energy minimization (tolerance=500 kJ/mol/nm)
            print("  Running energy minimization...")
            simulation.minimizeEnergy(tolerance=500 * unit.kilojoule_per_mole / unit.nanometer)
            print("  Energy minimization completed")

            # Get minimized state and energy
            state_after_min = simulation.context.getState(getPositions=True, getEnergy=True)
            energy_after_min = state_after_min.getPotentialEnergy()
            print(f"  Minimization energy: {energy_after_min}")

            # Check energy for NaN
            if np.isnan(energy_after_min.value_in_unit(unit.kilojoule_per_mole)):
                raise RuntimeError("Energy minimization failed: NaN energy")

            # Save structure after minimization
            positions_after_min = state_after_min.getPositions()
            debug_pdb_after_min = os.path.join(output_dir, 'after_minimization.pdb')
            with open(debug_pdb_after_min, 'w') as f:
                PDBFile.writeFile(
                    simulation.topology,
                    positions_after_min,
                    f
                )
            print(f"  Saving structure after minimization: {debug_pdb_after_min}")

            # Run simulation
            wfreq = self.config.simulation.wfreq
            xtc_file = os.path.join(output_dir, 'trajectory.xtc')
            log_file = os.path.join(output_dir, 'simulation.log')

            # Add reporters
            simulation.reporters.append(
                XTCReporter(xtc_file, wfreq)
            )
            simulation.reporters.append(
                StateDataReporter(
                    log_file,
                    wfreq,
                    step=True,
                    potentialEnergy=True,
                    kineticEnergy=True,
                    totalEnergy=True,
                    temperature=True,
                    volume=True,
                )
            )

            # Run simulation (show progress with tqdm)
            total_steps = self.config.simulation.steps
            print(f"  Starting simulation: {total_steps} steps")

            from tqdm import tqdm
            n_batches = 10
            batch_size = total_steps // n_batches

            for _ in tqdm(range(n_batches), desc="Mpipi-Recharged"):
                simulation.step(batch_size)
                simulation.saveCheckpoint(os.path.join(output_dir, 'restart.chk'))

            # Handle remaining steps
            remaining = total_steps % n_batches
            if remaining > 0:
                simulation.step(remaining)

            print(f"  Simulation finished!")

            # Get final state (including PBC info)
            state_final = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=True,
                getEnergy=True,
                enforcePeriodicBox=True
            )

            # Get final positions and box vectors
            positions_final = state_final.getPositions()
            box_vectors = state_final.getPeriodicBoxVectors()

            # Set box vectors on topology
            simulation.topology.setPeriodicBoxVectors(box_vectors)

            # Save final structure as PDB (including bonds - CONECT records)
            final_pdb = os.path.join(output_dir, 'final.pdb')
            final_pdb_mpipi_format = os.path.join(output_dir, 'final_mpipi_format.pdb')
            
            # Save mpipi-format PDB first (pA, pG, etc.)
            with open(final_pdb_mpipi_format, 'w') as f:
                PDBFile.writeFile(
                    simulation.topology,
                    positions_final,
                    f,
                    keepIds=True
                )
            print(f"  Saving mpipi-format PDB: {final_pdb_mpipi_format}")

            # Post-process: convert mpipi-format PDB to calvados format (CA + three-letter codes)
            print(f"\n  Post-processing PDB format conversion...")
            try:
                final_pdb_calvados_format = self._convert_mpipi_pdb_to_calvados_format(
                    mpipi_pdb=final_pdb_mpipi_format,
                    output_pdb=final_pdb
                )
                print(f"  ✓ PDB format conversion complete: {final_pdb_calvados_format}")
                print(f"  - Original mpipi format: {final_pdb_mpipi_format}")
                print(f"  - Calvados format: {final_pdb}")
            except Exception as e:
                print(f"  ⚠️  PDB format conversion failed: {e}")
                print(f"  Keeping original mpipi-format PDB: {final_pdb_mpipi_format}")
                # If conversion fails, copy mpipi format as final.pdb
                import shutil
                shutil.copy2(final_pdb_mpipi_format, final_pdb)
                import traceback
                traceback.print_exc()

            # Copy final.pdb to root directory (self.output_dir = {system_name}_CG)
            import shutil
            final_pdb_root = os.path.join(self.output_dir, 'final.pdb')
            shutil.copy2(final_pdb, final_pdb_root)
            print(f"  Copying final structure to root directory: {final_pdb_root}")

            # Copy file to standard location
            shutil.copy(final_pdb, os.path.join(output_dir, 'preequil_final.pdb'))

            # Save system XML
            system_xml = os.path.join(output_dir, 'system.xml')
            with open(system_xml, 'w') as f:
                f.write(mm.XmlSerializer.serialize(system))

            # Save checkpoint
            shutil.copy(os.path.join(output_dir, 'restart.chk'), os.path.join(output_dir, 'final.chk'))

            # Print result summary
            print(f"\n  ✓ Simulation completed!")
            print(f"  Final structure: {final_pdb}")
            print(f"  Trajectory: {xtc_file}")
            
            # Print final energy
            final_energy = state_final.getPotentialEnergy()
            print(f"  Final potential energy: {final_energy}")
            print(f"  Final energy per particle: {final_energy.value_in_unit(unit.kilojoule_per_mole) / system.getNumParticles():.3f} kJ/mol")

            result.success = True
            result.output_dir = output_dir
            result.trajectory = xtc_file
            result.structure = final_pdb
            result.errors = []

        except Exception as e:
            import traceback
            print(f"\n  ✗ Mpipi-Recharged simulation failed: {e}")
            if kwargs.get('verbose', False):
                traceback.print_exc()
            
            result.success = False
            result.output_dir = output_dir
            result.errors = [str(e)]

        finally:
            self.is_running = False

        return result
    
    def _build_globular_indices_dict(self) -> Dict[str, list]:
        """
        Build globular_indices_dict from config.components

        Used for the OpenMpipi get_mpipi_system function.
        OpenMpipi expects format: {chain_id: [[start1, end1], [start2, end2], ...]}
        Each [start, end] is an inclusive index range of a folded domain.

        Note: OpenMpipi expects **local indices** (0-based per chain),
        not global system indices.

        Returns:
            Dictionary mapping chain_id to list of domain ranges [start, end] (local indices)
        """
        globular_indices_dict = {}

        # Get list of chain IDs
        chain_ids = self.get_chain_identifiers()

        # Get folded domain information
        folded_domains = self.get_folded_domains()

        # Build dict: {chain_id: [[start1, end1], [start2, end2], ...]}
        # First collect start positions for each chain (local index baseline)
        chain_local_start = {}  # {chain_id: local_index_offset}
        
        for res_idx, chain_id in enumerate(chain_ids):
            if chain_id not in chain_local_start:
                chain_local_start[chain_id] = res_idx

        # Traverse each residue to detect continuous folded regions (domain ranges)
        current_chain = None
        domain_start = None  # Global index
        domain_start_local = None  # Local index

        for res_idx, (chain_id, is_folded) in enumerate(zip(chain_ids, folded_domains)):
            if chain_id not in globular_indices_dict:
                globular_indices_dict[chain_id] = []

            # New chain starts, reset state
            if current_chain != chain_id:
                if current_chain is not None and domain_start is not None:
                    # Save previous domain (using local indices)
                    local_start = domain_start - chain_local_start[current_chain]
                    local_end = (res_idx - 1) - chain_local_start[current_chain]
                    globular_indices_dict[current_chain].append([local_start, local_end])
                current_chain = chain_id
                domain_start = None

            # If folded domain, record start position
            if is_folded:
                if domain_start is None:
                    domain_start = res_idx
            else:
                # If was in domain, now ended; save domain range
                if domain_start is not None:
                    local_start = domain_start - chain_local_start[chain_id]
                    local_end = (res_idx - 1) - chain_local_start[chain_id]
                    globular_indices_dict[chain_id].append([local_start, local_end])
                    domain_start = None

        # Handle last domain (if chain ends folded)
        if current_chain is not None and domain_start is not None:
            local_start = domain_start - chain_local_start[current_chain]
            # Local end index is the last residue of the chain
            chain_length = sum(1 for cid in chain_ids if cid == current_chain)
            local_end = chain_length - 1
            globular_indices_dict[current_chain].append([local_start, local_end])

        return globular_indices_dict

    def _prepare_mpipi_output(self) -> Dict[str, str]:
        """
        Prepare Mpipi-Recharged output directory

        Output directory structure: {output_dir}/Mpipi-Recharged/
        CLI already sets self.output_dir = {system_name}_CG; use directly

        Expected structure:
        {system_name}_CG/
        ├── Mpipi-Recharged/     # Main simulation output
        │   ├── trajectory.xtc
        │   ├── final.pdb
        │   └── ...
        ├── final.pdb            # Copied to root directory
        └── equilibration/       # Pre-equilibration output (created by _run_pre_equilibration)
            └── raw/
                └── ...
        """
        # CLI already sets self.output_dir = {system_name}_CG; use directly
        # Mpipi-Recharged main simulation output directory
        mpipi_dir = os.path.join(self.output_dir, 'Mpipi-Recharged')

        # Backup old results
        import shutil
        from datetime import datetime

        if os.path.exists(mpipi_dir):
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{mpipi_dir}_backup_{timestamp}"
            shutil.move(mpipi_dir, backup_dir)
            print(f"  📁 Backing up old results to: {backup_dir}")

        os.makedirs(mpipi_dir, exist_ok=True)

        return {
            'output_dir': mpipi_dir,
            'task_name': 'Mpipi-Recharged',
        }
    
    # -------------------------------------------------------------------------
    # Utility Methods
    # -------------------------------------------------------------------------

    def _convert_mpipi_pdb_to_calvados_format(self, mpipi_pdb: str, output_pdb: str) -> str:
        """
        Convert mpipi_recharged-format PDB to calvados format (CA + three-letter codes)

        Use the standardization function in the backmap module.

        Args:
            mpipi_pdb: mpipi_recharged output PDB file path
            output_pdb: output PDB file path (overwrites)

        Returns:
            Output PDB file path
        """
        from .backmap import standardize_pdb_with_calvados
        return standardize_pdb_with_calvados(mpipi_pdb, self.config, output_pdb)

    def get_result(self) -> Optional[SimulationResult]:
        """Get latest simulation results"""
        return self._result

    def cleanup(self):
        """Clean up temporary files"""
        self.is_setup = False
        self._result = None

    # -------------------------------------------------------------------------
    # Pre-equilibration Utility Methods
    # -------------------------------------------------------------------------

    def run_pre_equilibration(
        self,
        gpu_id: int = 0,
        steps: int = 100000,
        mapping: str = "ca",
        k_restraint: float = 10000.0,
        use_com: bool = True,
        platform: Optional[ComputePlatform] = None,
    ) -> Optional[str]:
        """
        Run pre-equilibration only (build initial structure with CALVADOS)

        This method allows users to pre-generate CG structure without running a full simulation.
        The generated `preequil_final.pdb` file can be used for subsequent force-field simulations.

        Args:
            gpu_id: GPU device ID
            steps: pre-equilibration steps (default 100000)
            mapping: mapping method ('ca' or 'com') (default 'ca')
            k_restraint: restraint force constant (kJ/(mol·nm²)) (default 10000.0)
            use_com: whether to use COM restraints (default True)
            platform: compute platform (CUDA or CPU), default from config

        Returns:
            Path to pre-equilibrated structure file, or None if no MDP component
        """
        # If platform not specified, use config value (default CUDA)
        if platform is None:
            platform = self.config.simulation.platform

        return self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=steps,
            mapping=mapping,
            k_restraint=k_restraint,
            use_com=use_com,
            platform=platform,
        )

    def get_pre_equilibrated_structure(self) -> Optional[str]:
        """
        Get path to pre-equilibrated structure file

        Returns:
            Path to pre-equilibrated structure, or None if pre-equil not run
        """
        if self.output_dir is None:
            return None

        preequil_pdb = os.path.join(self.output_dir, 'preequil_final.pdb')
        if os.path.exists(preequil_pdb):
            return preequil_pdb
        return None


# =============================================================================
# Module Exports
# =============================================================================

__all__ = [
    # Configuration
    'CGSimulationConfig',
    'CGComponent',
    'ComponentType',
    'TopologyType',
    'ComputePlatform',
    'SimulationParams',

    # Result
    'SimulationResult',

    # Topology Info
    'TopologyInfo',

    # Simulator
    'CGSimulator',

    # PDB Tools (chain labeling)
    'ChainLabel',
]

