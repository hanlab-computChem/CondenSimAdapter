#!/usr/bin/env python3
"""
Adapter CLI

A workflow for CG and AA protein condensate simulation.

    Core Commands:
        cg               Run coarse-grained simulation
        backmap          Backmap CG structure to all-atom representation
        minimize         Energy minimization with AMBER/CHARMM force fields
        to_run           Generate production run scripts for minimize output

    Utility Commands:
        init             Initialize a new configuration template
        droplet-density  Estimate protein density in droplet geometry
        info             Display system and environment information
"""

import sys

import click

from .commands import init_command, cg_command, backmap_command, pace_opt_command, minimize_command, info_command, droplet_density_command, to_run_command, forcefield_command


class _LazyGroup(click.MultiCommand):
    """Lazily load a Click group to avoid heavy imports at CLI startup."""

    def __init__(self, import_path: str, **kwargs):
        super().__init__(**kwargs)
        self._import_path = import_path
        self._loaded = None

    def _load(self):
        if self._loaded is None:
            module_path, attr = self._import_path.rsplit(".", 1)
            import importlib
            mod = importlib.import_module(module_path)
            self._loaded = getattr(mod, attr)
        return self._loaded

    def list_commands(self, ctx):
        return self._load().list_commands(ctx)

    def get_command(self, ctx, name):
        return self._load().get_command(ctx, name)


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    """\b
    Adapter: a workflow for CG and AA protein condensate simulation.

    \b
    CORE COMMANDS:
        cg               Run coarse-grained simulation
        backmap          Backmap CG structure to all-atom representation
        minimize         Energy minimization with AMBER/CHARMM force fields
        to_run           Generate production run scripts for minimize output

    \b
    UTILITY COMMANDS:
        init             Initialize a new configuration template
        droplet-density  Estimate protein density in droplet geometry
        info             Display system and environment information
        forcefield       Manage custom all-atom force fields

    \b
    Available CG force fields:
        calvados, hps_urry, cocomo, mpipi_recharged

    \b
    Available all-atom force fields:
        1-a99SBdisp, 2-amber03wsc, 3-amber99sbws-stqp, 4-amber99sbws-stq,
        5-des-amber, 6-des-amber-sf1.0, 7-amber99sb-ildn, 8-amber14sb,
        9-charmm36m

    \b
    Typical workflow:
        1. adapter init my_project              # Create configuration template
        2. adapter cg -f config.yaml            # Run CG simulation
        3. adapter backmap -i output_CG -f config.yaml  # Backmap to all-atom
        4. adapter minimize -i output_backmap -f config.yaml  # Minimize
        5. adapter to_run -f config.yaml        # Generate production run scripts

    \b
    Additional examples:
        adapter init --topol droplet -c IIIMII  # Multi-component droplet
        adapter droplet-density -f config.yaml -r 15  # Estimate density
        adapter info                            # Check environment
    """
    pass


# Add commands (ordered: core commands first, then utility commands)
# Core commands
main.add_command(cg_command, 'cg')
main.add_command(backmap_command, 'backmap')
main.add_command(minimize_command, 'minimize')
main.add_command(to_run_command, 'to_run')
main.add_command(forcefield_command, 'forcefield')

# Utility commands
main.add_command(init_command, 'init')
main.add_command(droplet_density_command, 'droplet-density')
main.add_command(info_command, 'info')
_models_lazy = _LazyGroup(
    "CondenSimAdapter.cli.commands_refactored.models_command.models_command",
    name="models",
    help="Manage neural network models for backmapping.",
)
main.add_command(_models_lazy, 'models')

# Hidden/experimental commands
pace_opt_command.hidden = True
main.add_command(pace_opt_command, 'pace-opt')


def cli():
    """Entry point for CLI."""
    sys.exit(main())


if __name__ == '__main__':
    cli()

