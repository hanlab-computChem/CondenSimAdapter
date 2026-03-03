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

        try:
            from .cg2all import convert_cg2all
            convert_cg2all(
                in_pdb_fn=cg_pdb,
                out_fn=output_pdb,
                model_type=model_type,
                device=device,
                fix_atom=fix_atom,
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
    """
    Re-centre a slab PDB so the protein COM sits at z = 0.

    Reads the PDB, shifts all ATOM/HETATM coordinates, overwrites in place.
    Returns the same path.
    """
    lines_in = Path(pdb_path).read_text().splitlines()
    coords_z = []
    for line in lines_in:
        if line.startswith(("ATOM", "HETATM")):
            try:
                coords_z.append(float(line[46:54]))
            except ValueError:
                pass

    if not coords_z:
        return pdb_path

    z_shift = -np.mean(coords_z)
    out_lines = []
    for line in lines_in:
        if line.startswith(("ATOM", "HETATM")):
            try:
                z = float(line[46:54]) + z_shift
                line = line[:46] + f"{z:8.3f}" + line[54:]
            except (ValueError, IndexError):
                pass
        out_lines.append(line)

    Path(pdb_path).write_text("\n".join(out_lines))
    return pdb_path
