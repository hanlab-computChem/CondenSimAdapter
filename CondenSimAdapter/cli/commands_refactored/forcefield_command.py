#!/usr/bin/env python3
"""
Forcefield Command Group

Manage custom all-atom force fields for adapter minimize workflow.
"""

import click

from ...forcefield.registry import REGISTRY, BUILTIN_FORCE_FIELDS


@click.group("forcefield", context_settings={"help_option_names": ["-h", "--help"]})
def forcefield_command():
    """Manage all-atom force fields (add/remove/list custom entries)."""
    pass


@forcefield_command.command("list")
def list_forcefields():
    """List built-in and custom all-atom force fields."""
    click.echo("\nAll-atom force fields:")
    click.echo("-" * 60)

    for ff in BUILTIN_FORCE_FIELDS:
        click.echo(f"{ff.name:18} | builtin | {ff.pdb2gmx_name}")

    for ff_name in REGISTRY.list_custom_force_fields():
        ff = REGISTRY.get_force_field(ff_name)
        if ff:
            click.echo(f"{ff.name:18} | custom  | {ff.pdb2gmx_name}")

    click.echo("-" * 60)
    click.echo(f"Total: {REGISTRY.count()}")


@forcefield_command.command("add")
@click.option(
    "--ff-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=str),
    required=True,
    help="Path to custom force field directory (*.ff).",
)
@click.option(
    "--pdb2gmx-name",
    type=str,
    default=None,
    help="Name used for gmx pdb2gmx -ff (default: inferred from ff-dir).",
)
@click.option(
    "--family",
    type=click.Choice(["AMBER", "CHARMM"], case_sensitive=False),
    default="AMBER",
    show_default=True,
    help="Force field family.",
)
@click.option(
    "--water-model",
    type=str,
    default="tip3p",
    show_default=True,
    help="Water model for pdb2gmx.",
)
@click.option(
    "--solvate-cs",
    type=str,
    default="spc216",
    show_default=True,
    help="Water model used by gmx solvate -cs.",
)
@click.option(
    "--description",
    type=str,
    default="",
    help="Optional human-readable description.",
)
def add_forcefield(ff_dir, pdb2gmx_name, family, water_model, solvate_cs, description):
    """Register a custom all-atom force field and assign an aN id."""
    family_normalized = family.upper()
    gbsa_mapping = "AMBER99SB-ILDN" if family_normalized == "AMBER" else "CHARMM36"
    try:
        ff = REGISTRY.register_custom_force_field(
            ff_dir=ff_dir,
            pdb2gmx_name=pdb2gmx_name,
            family=family_normalized,
            water_model=water_model,
            solvate_cs=solvate_cs,
            description=description,
        )
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Registered custom force field: {ff.name}")
    click.echo(f"  pdb2gmx name: {ff.pdb2gmx_name}")
    click.echo(f"  family: {ff.family}")
    click.echo(f"  water model: {ff.water_model}")


@forcefield_command.command("remove")
@click.argument("force_field_id", type=str)
def remove_forcefield(force_field_id):
    """Remove a custom all-atom force field by aN id."""
    try:
        ff = REGISTRY.remove_custom_force_field(force_field_id)
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Removed custom force field: {ff.name} ({ff.pdb2gmx_name})")
