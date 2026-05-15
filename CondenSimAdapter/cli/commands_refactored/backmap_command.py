#!/usr/bin/env python3
"""Backmap command: CG -> all-atom."""

import sys
from pathlib import Path
from typing import Optional

import click


def _model_help() -> str:
    """Lazily resolve SUPPORTED_MODELS for help text."""
    try:
        from ...backmap.backmapper import SUPPORTED_MODELS

        return f"cg2all model type. Options: {', '.join(SUPPORTED_MODELS)}."
    except ImportError:
        return "cg2all model type (e.g. CalphaBasedModel)."


@click.command("backmap", context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--input",
    "-i",
    "input_path",
    type=click.Path(),
    default=None,
    help=(
        "Input CG PDB file or CG output directory "
        "(final.pdb used automatically).  "
        "If omitted, defaults to {system_name}_CG when -f is given."
    ),
)
@click.option(
    "--input-file",
    "-f",
    type=click.Path(exists=True),
    default=None,
    help="Simulation config YAML (for topology type and system name).",
)
@click.option(
    "--output",
    "-o",
    "output_dir",
    type=click.Path(),
    default=None,
    help="Output directory (default: {system_name}_backmap).",
)
@click.option(
    "--model-type",
    "-m",
    type=str,
    default="CalphaBasedModel",
    show_default=True,
    help=_model_help(),
)
@click.option(
    "--device",
    "-d",
    type=str,
    default="cpu",
    show_default=True,
    help="Compute device: cpu | cuda | cuda:0 etc.",
)
@click.option("--verbose", "-v", is_flag=True, default=False)
def backmap_command(
    input_path: Optional[str],
    input_file: Optional[str],
    output_dir: Optional[str],
    model_type: str,
    device: str,
    verbose: bool,
):
    """\b
    Convert a CG structure to all-atom via the cg2all neural network.

    \b
    Examples:
        adapter backmap -f FUS_LC.yaml              # auto-find FUS_LC_CG/
        adapter backmap -i FUS_LC_CG                # explicit directory
        adapter backmap -i FUS_LC_CG/final.pdb      # explicit file
        adapter backmap -i FUS_LC_CG -m Martini3
        adapter backmap -i FUS_LC_CG -f FUS_LC.yaml -d cuda
    """
    try:
        from ...backmap.backmapper import Backmapper
        from ...minimize.config_loader import load_config_from_yaml
    except ImportError as e:
        click.echo(f"Import error: {e}", err=True)
        sys.exit(1)

    # --- Resolve system_name and topology_type from config ----------------------
    topology_type = None
    system_name = None
    if input_file:
        try:
            import yaml

            name, _ = load_config_from_yaml(input_file)
            system_name = name
            with open(input_file) as fh:
                raw = yaml.safe_load(fh)
            topology_type = raw.get("topol") or raw.get("topology")
        except Exception as e:
            if verbose:
                click.echo(f"  Warning: could not parse config: {e}")

    # --- Auto-find input when -i is omitted -------------------------------------
    if input_path is None:
        if system_name is None:
            click.echo(
                "  Error: provide -i (input path) or -f (config YAML) so the "
                "CG directory can be located.",
                err=True,
            )
            sys.exit(1)
        input_path = f"{system_name}_CG"
        click.echo(f"  No -i given; trying default: {input_path}/")

    # --- Resolve directory → final.pdb -----------------------------------------
    # Track system_name from the directory name BEFORE resolving to final.pdb.
    p = Path(input_path).resolve()
    if p.is_dir():
        # Derive system_name from the directory name (strip trailing _CG suffix)
        dir_stem = p.name
        if system_name is None:
            system_name = dir_stem[:-3] if dir_stem.endswith("_CG") else dir_stem
        candidate = p / "final.pdb"
        if not candidate.exists():
            click.echo(
                f"  Error: directory '{input_path}' contains no final.pdb. "
                "Run 'adapter cg' first or specify the PDB file directly.",
                err=True,
            )
            sys.exit(1)
        cg_pdb = str(candidate)
    else:
        cg_pdb = str(p)
        if system_name is None:
            system_name = p.stem  # only for plain file input

    click.echo(f"\n{'=' * 60}\nBackmapping\n{'=' * 60}")
    click.echo(f"  Input:      {cg_pdb}")
    click.echo(f"  Model type: {model_type}")
    click.echo(f"  Device:     {device}")
    if topology_type:
        topo_note = {
            "slab": "z-centering enabled",
            "droplet": "droplet centering enabled",
            "cubic": "periodic boundary",
        }.get(topology_type, "")
        topo_str = f"  ({topo_note})" if topo_note else ""
        click.echo(f"  Topology:   {topology_type}{topo_str}")

    out = output_dir or f"{system_name}_backmap"

    backmapper = Backmapper()
    result = backmapper.run(
        cg_pdb=cg_pdb,
        output_dir=out,
        model_type=model_type,
        device=device,
        topology_type=topology_type,
    )

    if result.success:
        click.echo(f"\n  Completed. Output PDB: {result.output_pdb}")
    else:
        click.echo(f"\n  Backmapping failed: {result.error}", err=True)
        if verbose and result.error:
            click.echo(result.error, err=True)
        sys.exit(1)
