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
    solvated_top: str = ""
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
            import click

            # Step 1 — topology
            click.echo(f"\n  [1/4] Building GROMACS topology ({self._ff_name}) ...")
            ff_path = self._get_ff_path()
            if ff_path and Path(ff_path).exists():
                target_ff = out / Path(ff_path).name
                if not target_ff.exists():
                    shutil.copytree(ff_path, target_ff)
            else:
                log.warning(f"Force field folder not found for '{self.config.forcefield_type}'")

            from .topology_builder import generate_all_atom_topology, merge_topologies
            total_nmol = sum(getattr(c, "nmol", 1) for c in self.components)
            his_repeat = max(total_nmol * 30, 30)

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

            # Step 2 — structure GRO
            click.echo(f"  [2/4] Processing input structure ...")
            from .topology_builder import run_pdb2gmx_for_structure
            structure_gro = run_pdb2gmx_for_structure(
                Path(input_pdb), structure_dir, self._ff_name,
                water_model="none",
                disable_disulfide=self.config.disable_disulfide,
                his_type=self.config.his_type,
                his_repeat_count=his_repeat,
            )

            if self.config.box_resize and self.config.box_resize_dims:
                structure_gro = self._resize_box(
                    structure_gro, structure_dir, self.config.box_resize_dims
                )

            # Step 3 — OpenMM softcore minimization (implicit solvent)
            click.echo(f"  [3/4] OpenMM softcore minimization (3 stages) ...")
            min_result = self._run_openmm_minimization(
                str(structure_gro), str(merged_top), minimize_dir
            )

            final_pdb = min_result.get("final_pdb", "")
            result.success     = bool(final_pdb and Path(final_pdb).exists())
            result.intermediate_files = min_result

            if result.success:
                # Promote key output files to the root output dir (like old minimize.py)
                root_pdb = out / "minimize_final.pdb"
                root_top = out / "topol.top"
                shutil.copy2(final_pdb, root_pdb)
                shutil.copy2(merged_top, root_top)
                final_pdb = str(root_pdb)
                result.output_pdb = final_pdb
                click.echo(f"  Minimization done  →  {root_pdb.name}")

                # Step 4 — optional droplet box / solvation
                if self.config.droplet_box_type:
                    click.echo(f"  [4/4] Building droplet box ({self.config.droplet_box_type}) ...")
                    droplet_gro = self._build_droplet_box(
                        final_pdb, out / "droplet",
                        self.config.droplet_box_type,
                        self.config.droplet_distance,
                    )
                    # Promote to root dir
                    root_box = out / "minimize_final_box.gro"
                    shutil.copy2(droplet_gro, root_box)
                    final_pdb = str(root_box)
                    result.output_pdb = final_pdb
                    click.echo(f"  Droplet box done  →  {root_box.name}")

                if self.config.solvate:
                    click.echo(f"  [4/4] Explicit solvation ({water_model}, {self.config.ion_concentration} M) ...")
                    solvated_gro, solvated_top = self.solvate_system(
                        final_pdb, str(root_top),
                        out / "solvate", self.config.ion_concentration,
                    )
                    # Promote solvated outputs to root dir
                    root_solvated_gro = out / "minimize_final_solvated.gro"
                    root_solvated_top = out / "topol.top"
                    shutil.copy2(solvated_gro, root_solvated_gro)
                    shutil.copy2(solvated_top, root_solvated_top)
                    result.output_pdb  = str(root_solvated_gro)
                    result.solvated_top = str(root_solvated_top)
                    click.echo(f"  Solvation done  →  {root_solvated_gro.name}")

        except Exception as exc:
            import traceback
            result.errors.append(str(exc))
            result.errors.append(traceback.format_exc())

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
        structure_pdb: str,
        topology_top: str,
        solvate_dir: Path,
        ion_concentration: float = 0.15,
    ):
        """Add explicit solvent and ions to a minimized (or droplet-boxed) structure.

        The input is always a PDB (final.pdb from OpenMM minimization or the
        droplet-box GRO).  The box is already the right size from the CG run
        (slab/cubic) or was rebuilt by _build_droplet_box (droplet).  We must
        NOT call editconf with -d/-bt/-c here — that would create an enormous
        new box for slab/multi-chain systems.
        """
        solvate_dir.mkdir(parents=True, exist_ok=True)

        # Resolve -cs (solvent coordinate file) for this force field.
        # a99SBdisp / DES-AMBER / amber03wsc → tip4p
        # amber99sb-ildn / amber14sb / charmm36m → spc216
        from ..forcefield.registry import get_force_field
        ff_info = get_force_field(self.config.forcefield_type)
        solvate_cs = ff_info.solvate_cs if ff_info else "spc216"
        log.info(f"  Solvating with -cs {solvate_cs} ({self.config.forcefield_type})")

        # Convert PDB → GRO (format only, no box modification).
        system_gro = solvate_dir / "system.gro"
        subprocess.run(
            ["gmx", "editconf", "-f", structure_pdb, "-o", str(system_gro)],
            check=True, capture_output=True, text=True,
        )

        # solvate with the force-field-specific solvent coordinate file
        solvated_gro = str(solvate_dir / "solvated.gro")
        subprocess.run(
            ["gmx", "solvate", "-cp", str(system_gro), "-o", solvated_gro,
             "-cs", solvate_cs, "-p", topology_top],
            check=True, capture_output=True, text=True,
        )

        # Copy topology and force field folder into solvate_dir so that
        # gmx grompp can resolve relative #include paths (e.g. a99SBdisp.ff/...)
        # when run with cwd=solvate_dir.  This mirrors the old minimize.py logic.
        local_top = solvate_dir / "topol.top"
        shutil.copy2(topology_top, local_top)
        try:
            from ..forcefield.registry import get_force_field_path
            ff_path = get_force_field_path(self.config.forcefield_type)
            if ff_path:
                ff_dest = solvate_dir / Path(ff_path).name
                if not ff_dest.exists():
                    shutil.copytree(ff_path, ff_dest)
        except Exception:
            pass

        # Add ions — run grompp in solvate_dir so relative FF #includes resolve
        em_mdp = solvate_dir / "ions.mdp"
        em_mdp.write_text("integrator = steep\nnsteps = 0\n")
        ions_tpr = str(solvate_dir / "ions.tpr")
        subprocess.run(
            ["gmx", "grompp", "-f", "ions.mdp", "-c", solvated_gro,
             "-p", str(local_top), "-o", ions_tpr, "-maxwarn", "5"],
            check=True, capture_output=True, text=True,
            cwd=str(solvate_dir),
        )
        ionized_gro = str(solvate_dir / "ionized.gro")
        subprocess.run(
            ["gmx", "genion", "-s", ions_tpr, "-o", ionized_gro,
             "-p", str(local_top), "-pname", "NA", "-nname", "CL",
             "-conc", str(ion_concentration), "-neutral"],
            input="SOL\n", check=True, capture_output=True, text=True,
            cwd=str(solvate_dir),
        )
        return Path(ionized_gro), Path(local_top)

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
        Run three-stage softcore minimization directly in the current process.

        Runs in-process (no subprocess) to guarantee that the same CUDA/GPU
        environment visible to the parent is used without any re-discovery issues.
        """
        from ..forcefield.registry import get_force_field as _get_ff
        ff_info = _get_ff(self.config.forcefield_type)
        ff_type = ff_info.family.lower() if ff_info else "amber"
        ff_name = ff_info.pdb2gmx_name if ff_info else self.config.forcefield_type
        device  = self.config.platform.lower()

        from .openmm_runner import run_minimization
        run_minimization(
            input_gro=structure_gro,
            input_top=topology_top,
            output_dir=str(minimize_dir),
            device=device,
            gpu_id=self.config.gpu_id,
            max_iterations=self.config.max_iterations,
            ff_type=ff_type,
            ff_name=ff_name,
            gb_model=self.config.gb_model,
            salt_conc=self.config.ion_concentration,
            cutoff=self.config.nonbonded_cutoff,
            tolerance=self.config.tolerance,
        )

        final_pdb = str(minimize_dir / "final.pdb")
        return {"final_pdb": final_pdb, "minimize_dir": str(minimize_dir)}

    # ------------------------------------------------------------------
    # Force field resolution
    # ------------------------------------------------------------------

    def _resolve_ff_name(self, forcefield_type: str) -> str:
        """Convert registry key (e.g. '1-a99SBdisp') to pdb2gmx name (e.g. 'a99SBdisp')."""
        from ..forcefield.registry import get_force_field
        
        ff = get_force_field(forcefield_type)
        if ff:
            return ff.pdb2gmx_name
        
        # Fallback: strip leading numeric prefix
        parts = forcefield_type.split("-", 1)
        return parts[-1] if len(parts) > 1 else forcefield_type

    def _get_ff_path(self) -> Optional[str]:
        try:
            from ..forcefield.registry import get_force_field_path
            return get_force_field_path(self.config.forcefield_type)
        except Exception:
            return None
