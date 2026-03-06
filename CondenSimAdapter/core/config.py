"""
Data models for the unified CG simulation engine.

Protein-only: IDP (intrinsically disordered) and MDP (multi-domain with folded regions).
"""

from __future__ import annotations

import yaml
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class ComponentType(Enum):
    IDP = "IDP"   # intrinsically disordered protein
    MDP = "MDP"   # multi-domain protein with folded regions


class TopologyType(Enum):
    SLAB    = "slab"
    DROPLET = "droplet"
    CUBIC   = "cubic"


class ForceField(Enum):
    CALVADOS2 = "calvados2"
    CALVADOS3 = "calvados3"
    HPS       = "hps"
    COCOMO    = "cocomo"
    MPIPI     = "mpipi"


@dataclass
class Component:
    """A single protein component (one species, possibly multiple copies)."""
    name: str
    comp_type: ComponentType
    nmol: int = 1

    # Sequence / structure input -- at least one required
    sequence: Optional[str] = None        # one-letter amino acid sequence
    fasta_path: Optional[str] = None      # path to FASTA file
    pdb_path: Optional[str] = None        # required for MDP; CA template

    # Folded domains: list of (start, end) 1-based inclusive residue ranges
    folded_domains: List[Tuple[int, int]] = field(default_factory=list)

    def get_sequence(self) -> str:
        if self.sequence:
            return self.sequence
        if self.fasta_path:
            return _read_fasta(self.fasta_path, self.name)
        if self.pdb_path:
            return _seq_from_pdb(self.pdb_path)
        raise ValueError(f"Component '{self.name}': no sequence, fasta_path, or pdb_path provided.")

    @classmethod
    def from_dict(cls, d: dict) -> Component:
        comp_type = ComponentType(d.get("type", "IDP").upper())
        domains_raw = d.get("fdomains", [])
        folded_domains = _parse_fdomains(domains_raw)
        return cls(
            name=d["name"],
            comp_type=comp_type,
            nmol=int(d.get("nmol", 1)),
            sequence=d.get("sequence"),
            fasta_path=d.get("ffasta"),
            pdb_path=d.get("fpdb"),
            folded_domains=folded_domains,
        )


@dataclass
class SimulationParams:
    steps: int    = 100_000_000
    dt: float     = 0.01          # ps
    wfreq: int    = 10_000
    log_freq: int = 1_000_000
    friction: float = 0.01        # ps^-1
    platform: str = "CUDA"
    gpu_id: int   = 0

    @classmethod
    def from_dict(cls, d: dict) -> SimulationParams:
        return cls(
            steps   = int(d.get("steps", 100_000_000)),
            dt      = float(d.get("dt", 0.01)),
            wfreq   = int(d.get("wfreq", 10_000)),
            log_freq= int(d.get("log_freq", 1_000_000)),
            friction= float(d.get("friction", 0.01)),
            platform= d.get("platform", "CUDA"),
            gpu_id  = int(d.get("gpu_id", 0)),
        )


@dataclass
class CGConfig:
    """Complete configuration for one CG simulation run."""
    system_name: str
    force_field: str                    # ForceField.value string
    components: List[Component]
    box: List[float]                    # [Lx, Ly, Lz] in nm
    topology: TopologyType
    temperature: float     = 300.0     # K
    ionic_strength: float  = 0.15      # M
    simulation: SimulationParams = field(default_factory=SimulationParams)

    # Slab-specific (CALVADOS legacy default: 100 nm)
    slab_width: float = 100.0            # nm, matches CALVADOS default_config.yaml

    # Droplet-specific
    droplet_radius: Optional[float] = None   # nm
    droplet_k: float = 1.0                   # kJ/mol/nm^2 (confinement spring)

    @property
    def n_molecules(self) -> int:
        return sum(c.nmol for c in self.components)

    @property
    def resolved_force_field(self) -> str:
        """Resolve 'calvados' to the appropriate version based on component types.

        CALVADOS2 is parameterized for IDPs only.
        CALVADOS3 extends it to folded/mixed-disorder proteins (MDP).
        If the user simply writes force_field: calvados, the backend picks
        the correct version automatically.
        """
        if self.force_field == "calvados":
            has_mdp = any(c.comp_type == ComponentType.MDP for c in self.components)
            return "calvados3" if has_mdp else "calvados2"
        return self.force_field

    def get_component(self, name: str) -> Optional[Component]:
        for c in self.components:
            if c.name == name:
                return c
        return None

    @classmethod
    def from_dict(cls, d: dict) -> CGConfig:
        components = [Component.from_dict(c) for c in d.get("components", [])]
        topol_str = d.get("topol", d.get("topology", "cubic")).lower()
        # "grid" is CALVADOS's name for a cubic periodic box
        if topol_str == "grid":
            topol_str = "cubic"
        box_raw = d.get("box", [20.0, 20.0, 20.0])
        sim_raw = d.get("simulation", {})
        return cls(
            system_name   = d.get("system_name", d.get("sysname", "system")),
            force_field   = d.get("force_field", d.get("ff", "calvados")).lower(),
            components    = components,
            box           = [float(x) for x in box_raw],
            topology      = TopologyType(topol_str),
            temperature   = float(d.get("temperature", d.get("temp", 300.0))),
            ionic_strength= float(d.get("ionic_strength", d.get("ionic", 0.15))),
            simulation    = SimulationParams.from_dict(sim_raw),
            slab_width    = float(d.get("slab_width", 100.0)),
            droplet_radius= d.get("droplet_radius"),
            droplet_k     = float(d.get("droplet_k", 1.0)),
        )

    @classmethod
    def from_yaml(cls, path: str) -> CGConfig:
        with open(path) as f:
            raw = yaml.safe_load(f)
        return cls.from_dict(raw)


@dataclass
class SimulationResult:
    success: bool
    output_dir: str
    final_pdb: Optional[str] = None
    trajectory: Optional[str] = None
    log_file: Optional[str] = None
    error: Optional[str] = None
    elapsed_seconds: float = 0.0


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _read_fasta(path: str, name: Optional[str] = None) -> str:
    """Read a single protein sequence from a FASTA file."""
    seq_lines = []
    found = name is None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # If we already found and collected our target sequence, return it
                if found and seq_lines:
                    return "".join(seq_lines)
                header = line[1:]
                found = (name is None) or (name in header)
                seq_lines = []
            elif found:
                seq_lines.append(line.replace(" ", ""))
    if not seq_lines:
        raise ValueError(f"Sequence '{name}' not found in {path}")
    return "".join(seq_lines)


def _seq_from_pdb(pdb_path: str) -> str:
    """Extract CA-only sequence from a PDB file using MDAnalysis."""
    import MDAnalysis as mda
    u = mda.Universe(pdb_path)
    ca = u.select_atoms("name CA")
    from MDAnalysis.lib.util import convert_aa_code
    seq = ""
    for res in ca.residues:
        try:
            seq += convert_aa_code(res.resname)
        except Exception:
            seq += "G"
    return seq


def _parse_fdomains(raw) -> List[Tuple[int, int]]:
    """
    Parse folded-domain specifications.

    Accepts:
      - list of [start, end] pairs: [[1, 50], [80, 130]]
      - string: "1-50, 80-130"
      - dict with component name: {"CompName": [[1,50]]}
    """
    if not raw:
        return []
    if isinstance(raw, dict):
        # pick the first (or only) value
        raw = next(iter(raw.values()))
    if isinstance(raw, str):
        domains = []
        for part in raw.replace(";", ",").split(","):
            part = part.strip()
            if "-" in part:
                a, b = part.split("-", 1)
                domains.append((int(a), int(b)))
        return domains
    # list of [start, end]
    return [(int(p[0]), int(p[1])) for p in raw]
