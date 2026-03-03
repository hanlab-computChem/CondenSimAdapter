"""
All-atom energy minimizer.

Workflow:
  1. Copy AA force field to output directory (GROMACS needs it locally).
  2. Generate per-component topology with ``gmx pdb2gmx``.
  3. Merge multi-component topologies.
  4. Generate structure GRO from input PDB.
  5. (Optional) Resize / solvate box.
  6. Run three-stage OpenMM softcore minimization in a subprocess.

Replaces the 950-line src/minimize.py.
Changes vs old code:
  - No click.echo -- uses Python logging throughout.
  - resize_box() is self-contained (no longer delegates to BackmapSimulator).
  - MinimizeConfig and MinimizeResult interface are preserved for CLI compatibility.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class MinimizeConfig:
    """Configuration for the all-atom minimization workflow."""

    # AA force field name (must be registered in forcefield.registry)
    forcefield_type: str = "1-a99SBdisp"

    # Implicit-solvent GB model used in stage 3
    gb_model: str = "GBn2"                  # GBn2 | OBC2

    # Softcore schedule
    softcore_lambdas: List[float] = field(default_factory=lambda: [0.75, 0.85, 0.95])

    # OpenMM platform
    platform: str = "CUDA"
    gpu_id: int = 0

    # Minimisation convergence
    tolerance: float = 100.0               # kJ/(mol·nm)
    max_iterations: int = 5000

    # Non-bonded cutoff
    nonbonded_cutoff: float = 2.0          # nm

    # pdb2gmx options
    disable_disulfide: bool = False
    his_type: Optional[int] = None         # 0 or 1 (None = interactive auto)

    # Solvation (explicit water)
    solvate: bool = False
    ion_concentration: float = 0.15        # M

    # Post-minimization box building for droplet systems
    droplet_box_type: Optional[str] = None # dodecahedron | cubic | octahedron
    droplet_distance: float = 2.0          # nm

    # Box resize (optional pre-processing)
    box_resize: bool = False
    box_resize_dims: Optional[List[float]] = None


@dataclass
class MinimizeResult:
    """Result of the minimization workflow."""
    success: bool
    output_pdb: str = ""
    input_pdb: str = ""
    output_dir: str = ""
    errors: List[str] = field(default_factory=list)
    intermediate_files: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

class MinimizeSimulator:
    """
    Protein all-atom energy minimizer.

    Wraps GROMACS (topology / structure preparation) and OpenMM
    (three-stage softcore minimization).
    """

    def __init__(
        self,
        config: MinimizeConfig,
        components: list,      # list of Component (core.config) or legacy CGComponent
        system_name: str,
    ):
        self.config      = config
        self.components  = components
        self.system_name = system_name
        self._ff_name    = self._resolve_ff_name(config.forcefield_type)

    @classmethod
    def from_yaml(cls, yaml_path: str, config: Optional[MinimizeConfig] = None) -> MinimizeSimulator:
        """Construct from a CGConfig YAML file."""
        from .config_loader import load_config_from_yaml
        from ..core.config import Component
        system_name, raw_components = load_config_from_yaml(yaml_path)
        components = [Component.from_dict(c) for c in raw_components]
        return cls(config or MinimizeConfig(), components, system_name)

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(
        self,
        input_pdb: str,
        output_dir: Optional[str] = None,
    ) -> MinimizeResult:
        """
        Execute the minimization workflow.

        Args:
            input_pdb:   Path to the backmapped all-atom PDB.
            output_dir:  Output directory (defaults to {system_name}_minimize/).
        """
        out = Path(output_dir or f"{self.system_name}_minimize").resolve()
        out.mkdir(parents=True, exist_ok=True)

        topology_dir = out / "topology"
        structure_dir = out / "structure"
        minimize_dir  = out / "minimize"
        for d in (topology_dir, structure_dir, minimize_dir):
            d.mkdir(exist_ok=True)

        result = MinimizeResult(success=False, input_pdb=input_pdb, output_dir=str(out))

        try:
            # 1. Copy AA force field for GROMACS
            ff_path = self._get_ff_path()
            if ff_path and Path(ff_path).exists():
                target_ff = out / Path(ff_path).name
                if not target_ff.exists():
                    shutil.copytree(ff_path, target_ff)
                log.info(f"Using force field: {Path(ff_path).name}")
            else:
                log.warning(f"Force field folder not found for '{self.config.forcefield_type}'")

            # 2. Generate per-component topology
            from .topology_builder import generate_all_atom_topology, merge_topologies
            total_nmol = sum(getattr(c, "nmol", 1) for c in self.components)
            his_repeat = max(total_nmol * 30, 30)

            log.info("Generating per-component topologies ...")
            comp_tops, water_model = generate_all_atom_topology(
                self.components,
                self._ff_name,
                topology_dir,
                disable_disulfide=self.config.disable_disulfide,
                his_type=self.config.his_type,
                his_repeat_count=his_repeat,
            )

            merged_top = merge_topologies(
                comp_tops, topology_dir, self._ff_name, water_model, self.system_name
            )
            log.info(f"Merged topology: {merged_top.name}")

            # 3. Generate structure GRO from PDB
            from .topology_builder import run_pdb2gmx_for_structure
            log.info("Generating structure GRO ...")
            structure_gro = run_pdb2gmx_for_structure(
                Path(input_pdb), structure_dir, self._ff_name,
                water_model="none",
                disable_disulfide=self.config.disable_disulfide,
                his_type=self.config.his_type,
                his_repeat_count=his_repeat,
            )

            # 4. Optional: resize box
            if self.config.box_resize and self.config.box_resize_dims:
                structure_gro = self._resize_box(
                    structure_gro, structure_dir, self.config.box_resize_dims
                )

            # 5. Optional: solvate
            topology_top = merged_top
            if self.config.solvate:
                structure_gro, topology_top = self.solvate_system(
                    str(structure_gro), str(merged_top),
                    out / "solvate", self.config.ion_concentration,
                )

            # 6. OpenMM three-stage minimization
            log.info("Running OpenMM softcore minimization ...")
            min_result = self._run_openmm_minimization(
                str(structure_gro), str(topology_top), minimize_dir
            )

            final_pdb = min_result.get("final_pdb", "")
            result.success     = bool(final_pdb and Path(final_pdb).exists())
            result.output_pdb  = final_pdb
            result.intermediate_files = min_result

            if result.success:
                log.info(f"Minimization complete: {final_pdb}")

                # 7. Optional: build droplet box
                if self.config.droplet_box_type:
                    final_pdb = self._build_droplet_box(
                        final_pdb, out / "droplet",
                        self.config.droplet_box_type,
                        self.config.droplet_distance,
                    )
                    result.output_pdb = final_pdb

        except Exception as exc:
            log.exception("Minimization failed")
            result.errors.append(str(exc))

        return result

    # ------------------------------------------------------------------
    # Step helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _resize_box(
        gro_path: Path,
        out_dir: Path,
        dimensions: List[float],
    ) -> Path:
        """Resize the simulation box using ``gmx editconf``."""
        out_gro = out_dir / "resized.gro"
        dim_str = " ".join(str(d) for d in dimensions)
        subprocess.run(
            ["gmx", "editconf", "-f", str(gro_path), "-o", str(out_gro),
             "-box", *[str(d) for d in dimensions]],
            check=True, capture_output=True, text=True,
        )
        return out_gro

    def solvate_system(
        self,
        structure_gro: str,
        topology_top: str,
        solvate_dir: Path,
        ion_concentration: float = 0.15,
    ):
        """Add explicit solvent and ions (TIP3P)."""
        solvate_dir.mkdir(parents=True, exist_ok=True)

        # editconf -- add 1.2 nm box padding
        box_gro = str(solvate_dir / "box.gro")
        subprocess.run(
            ["gmx", "editconf", "-f", structure_gro, "-o", box_gro,
             "-c", "-d", "1.2", "-bt", "cubic"],
            check=True, capture_output=True, text=True,
        )

        # solvate
        solvated_gro = str(solvate_dir / "solvated.gro")
        subprocess.run(
            ["gmx", "solvate", "-cp", box_gro, "-o", solvated_gro,
             "-p", topology_top],
            check=True, capture_output=True, text=True,
        )

        # Add ions
        em_mdp = solvate_dir / "ions.mdp"
        em_mdp.write_text("integrator = steep\nnsteps = 0\n")
        ions_tpr = str(solvate_dir / "ions.tpr")
        subprocess.run(
            ["gmx", "grompp", "-f", str(em_mdp), "-c", solvated_gro,
             "-p", topology_top, "-o", ions_tpr, "-maxwarn", "5"],
            check=True, capture_output=True, text=True,
        )
        ionized_gro = str(solvate_dir / "ionized.gro")
        subprocess.run(
            ["gmx", "genion", "-s", ions_tpr, "-o", ionized_gro,
             "-p", topology_top, "-pname", "NA", "-nname", "CL",
             "-conc", str(ion_concentration), "-neutral"],
            input="SOL\n", check=True, capture_output=True, text=True,
        )
        return Path(ionized_gro), Path(topology_top)

    def _build_droplet_box(
        self,
        pdb_path: str,
        droplet_dir: Path,
        box_type: str,
        distance: float,
    ) -> str:
        """Build a GROMACS box for a droplet system after minimization."""
        droplet_dir.mkdir(parents=True, exist_ok=True)
        out_gro = str(droplet_dir / "droplet.gro")
        subprocess.run(
            ["gmx", "editconf", "-f", pdb_path, "-o", out_gro,
             "-c", "-d", str(distance), "-bt", box_type],
            check=True, capture_output=True, text=True,
        )
        return out_gro

    def _run_openmm_minimization(
        self,
        structure_gro: str,
        topology_top: str,
        minimize_dir: Path,
    ) -> dict:
        """
        Launch the three-stage softcore minimization as a subprocess.

        The worker script (minimize/worker.py) does the actual OpenMM work
        in a child process so it can be killed cleanly on error.
        """
        import sys
        worker = Path(__file__).parent / "worker.py"
        cmd = [
            sys.executable, str(worker),
            "--gro", structure_gro,
            "--top", topology_top,
            "--out", str(minimize_dir),
            "--ff", self.config.forcefield_type,
            "--gb", self.config.gb_model,
            "--tol", str(self.config.tolerance),
            "--max-iter", str(self.config.max_iterations),
            "--cutoff", str(self.config.nonbonded_cutoff),
            "--platform", self.config.platform,
            "--gpu", str(self.config.gpu_id),
            "--lambdas", ",".join(str(l) for l in self.config.softcore_lambdas),
        ]
        log.debug("minimize worker cmd: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            cwd=str(minimize_dir.parent),
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"OpenMM minimization failed:\n{proc.stderr[-2000:]}"
            )

        final_pdb = str(minimize_dir / "final.pdb")
        return {"final_pdb": final_pdb, "minimize_dir": str(minimize_dir)}

    # ------------------------------------------------------------------
    # Force field resolution
    # ------------------------------------------------------------------

    def _resolve_ff_name(self, forcefield_type: str) -> str:
        """Convert registry key (e.g. '1-a99SBdisp') to pdb2gmx name (e.g. 'a99SBdisp')."""
        # Strip leading numeric prefix if present
        parts = forcefield_type.split("-", 1)
        return parts[-1] if len(parts) > 1 else forcefield_type

    def _get_ff_path(self) -> Optional[str]:
        try:
            from ..forcefield.registry import get_force_field_path
            return get_force_field_path(self.config.forcefield_type)
        except Exception:
            return None
