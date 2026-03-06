"""
Simplified CG -> all-atom backmapper.

The entire logic reduces to:
  1. Validate input PDB exists.
  2. Call convert_cg2all().
  3. Optionally re-centre slab structures in z.

Supported model types (from cg2all checkpoints):
  CalphaBasedModel  -- standard CA-based CG (CALVADOS, HPS, COCOMO, Mpipi output)
  ResidueBasedModel -- residue-bead CG
  Martini           -- Martini 2
  Martini3          -- Martini 3
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

log = logging.getLogger(__name__)

# Valid model types exported by cg2all
SUPPORTED_MODELS = (
    "CalphaBasedModel",
    "ResidueBasedModel",
    "Martini",
    "Martini3",
)


@dataclass
class BackmapResult:
    success: bool
    output_pdb: Optional[str] = None
    input_pdb: Optional[str] = None
    model_type: Optional[str] = None
    error: Optional[str] = None


class Backmapper:
    """
    Minimal wrapper around convert_cg2all.

    All the old source-detection / config-hunting / directory-walking
    logic has been removed.  The caller provides explicit paths.
    """

    def run(
        self,
        cg_pdb: str,
        output_dir: str,
        model_type: str = "CalphaBasedModel",
        device: str = "cpu",
        fix_atom: bool = True,
        topology_type: Optional[str] = None,
    ) -> BackmapResult:
        """
        Backmap a CG PDB to all-atom representation.

        Args:
            cg_pdb:        Path to the CG PDB file (CA-only or Martini beads).
            output_dir:    Directory to write the all-atom PDB into.
            model_type:    cg2all model name (see SUPPORTED_MODELS).
            device:        'cpu' or 'cuda' / 'cuda:0' etc.
            fix_atom:      Fix CA position during reconstruction (CalphaBasedModel only).
            topology_type: Pass 'slab' to re-centre the output along z axis.

        Returns:
            BackmapResult with output_pdb path on success.
        """
        cg_pdb = str(cg_pdb)
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        if not Path(cg_pdb).exists():
            return BackmapResult(success=False, input_pdb=cg_pdb,
                                 error=f"Input file not found: {cg_pdb}")

        if model_type not in SUPPORTED_MODELS:
            log.warning(
                f"Unknown model type '{model_type}'. "
                f"Supported: {SUPPORTED_MODELS}. Proceeding anyway."
            )

        output_pdb = str(out_dir / "backmapped.pdb")

        # For CalphaBasedModel, fix_atom=True requires CalphaBasedModel-FIX.ckpt.
        # Auto-degrade to fix_atom=False if that checkpoint is absent.
        effective_fix = fix_atom
        if fix_atom and model_type == "CalphaBasedModel":
            from .cg2all.lib.libconfig import MODEL_HOME
            fix_ckpt = MODEL_HOME / "CalphaBasedModel-FIX.ckpt"
            if not fix_ckpt.exists():
                log.info(
                    "CalphaBasedModel-FIX.ckpt not found; "
                    "using normal checkpoint (fix_atom=False)."
                )
                effective_fix = False

        try:
            from .cg2all import convert_cg2all
            convert_cg2all(
                in_pdb_fn=cg_pdb,
                out_fn=output_pdb,
                model_type=model_type,
                device=device,
                fix_atom=effective_fix,
            )
        except Exception as exc:
            log.exception("cg2all conversion failed")
            return BackmapResult(
                success=False,
                input_pdb=cg_pdb,
                model_type=model_type,
                error=str(exc),
            )

        if topology_type and topology_type.lower() == "slab":
            try:
                output_pdb = _center_slab_in_z(output_pdb)
            except Exception as exc:
                log.warning(f"Slab re-centering failed (non-fatal): {exc}")

        return BackmapResult(
            success=True,
            output_pdb=output_pdb,
            input_pdb=cg_pdb,
            model_type=model_type,
        )


# ---------------------------------------------------------------------------
# Slab geometry helper
# ---------------------------------------------------------------------------

def _center_slab_in_z(pdb_path: str) -> str:
    """Re-centre a slab PDB so the protein COM sits at box_z/2.

    Algorithm (matches old backmap.py):
      1. Read box_z from the CRYST1 record (Angstrom).
      2. Compute mass-weighted protein COM in z via MDAnalysis (preferred)
         or fall back to unweighted mean over ATOM records.
      3. Shift all atoms so that COM_z → box_z/2.
      4. After MDAnalysis writes the file, fix C-terminal atom names
         (OT1→O, OT2→OXT) and re-insert TER records.
    """
    pdb_path = str(pdb_path)

    # --- read box_z from CRYST1 -------------------------------------------------
    box_z_ang: Optional[float] = None
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith("CRYST1"):
            try:
                box_z_ang = float(line[24:33])
            except ValueError:
                pass
            break

    # --- try MDAnalysis (mass-weighted, respects PBC) ---------------------------
    try:
        import MDAnalysis as mda  # type: ignore

        u = mda.Universe(pdb_path)
        protein = u.select_atoms("protein")
        if len(protein) == 0:
            protein = u.atoms

        com = protein.center_of_mass()
        z_com = com[2]  # Angstrom

        if box_z_ang is None and u.dimensions is not None:
            box_z_ang = u.dimensions[2]
        if box_z_ang is None:
            box_z_ang = 2.0 * z_com  # last resort

        offset = box_z_ang / 2.0 - z_com
        log.info(
            f"  Slab z-centering: COM_z={z_com/10:.2f} nm  "
            f"box_z={box_z_ang/10:.2f} nm  offset={offset/10:+.2f} nm"
        )

        positions = u.atoms.positions.copy()
        positions[:, 2] += offset
        u.atoms.positions = positions
        u.atoms.write(pdb_path)

        _fix_c_terminus_atom_names(pdb_path)
        _insert_ter_after_oxt(pdb_path)
        return pdb_path

    except ImportError:
        log.warning("MDAnalysis not available; falling back to plain-text z-centering")

    # --- plain-text fallback (no MDAnalysis) ------------------------------------
    lines_in = Path(pdb_path).read_text().splitlines()
    atom_z = []
    for line in lines_in:
        if line.startswith("ATOM"):
            try:
                atom_z.append(float(line[46:54]))
            except ValueError:
                pass

    if not atom_z:
        return pdb_path

    z_com = float(np.mean(atom_z))
    if box_z_ang is None:
        box_z_ang = 2.0 * z_com
    offset = box_z_ang / 2.0 - z_com

    out_lines = []
    for line in lines_in:
        if line.startswith(("ATOM", "HETATM")):
            try:
                z = float(line[46:54]) + offset
                line = line[:46] + f"{z:8.3f}" + line[54:]
            except (ValueError, IndexError):
                pass
        out_lines.append(line)

    Path(pdb_path).write_text("\n".join(out_lines))
    return pdb_path


def _fix_c_terminus_atom_names(pdb_path: str) -> None:
    """Fix C-terminal oxygen names written by MDAnalysis (OT1/OT2 → O/OXT)."""
    rename_map = {"OT1": " O  ", "OT2": " OXT"}
    try:
        lines = Path(pdb_path).read_text().splitlines()
    except OSError:
        return
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16]
            if atom_name.strip() in rename_map:
                line = line[:12] + rename_map[atom_name.strip()] + line[16:]
        out.append(line)
    Path(pdb_path).write_text("\n".join(out))


def _insert_ter_after_oxt(pdb_path: str) -> None:
    """Insert TER records after OXT atoms if not already present."""
    try:
        lines = Path(pdb_path).read_text().splitlines()
    except OSError:
        return
    out = []
    for i, line in enumerate(lines):
        out.append(line)
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "OXT":
            next_line = lines[i + 1] if i + 1 < len(lines) else ""
            if not next_line.startswith("TER"):
                out.append("TER")
    Path(pdb_path).write_text("\n".join(out))
