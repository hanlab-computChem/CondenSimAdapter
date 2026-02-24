#!/usr/bin/env python3
"""
Info Command

Displays system and environment information.
"""

import sys

import click

from ..shared import CG_FORCE_FIELDS, get_minimize_force_fields


@click.command('info', context_settings={'help_option_names': ['-h', '--help']})
def info_command():
    """
    Display system and environment information.
    """
    click.echo(f"\n{'=' * 60}")
    click.echo(f"Adapter Environment")
    click.echo(f"{'=' * 60}")
    
    # Python version
    click.echo(f"\n  Python: {sys.version}")
    
    # Available CG force fields
    click.echo(f"\n  Available CG force fields:")
    for ff in CG_FORCE_FIELDS:
        click.echo(f"    - {ff}")
    
    # Available all-atom force fields
    click.echo(f"\n  Available all-atom force fields:")
    aa_ffs = get_minimize_force_fields()
    # Display in two columns
    for i in range(0, len(aa_ffs), 2):
        left = f"    {aa_ffs[i]}"
        if i + 1 < len(aa_ffs):
            right = f"{aa_ffs[i+1]}"
            click.echo(f"{left:30} {right}")
        else:
            click.echo(left)
    
    # Check GPU availability
    try:
        import torch
        cuda_available = torch.cuda.is_available()
        if cuda_available:
            click.echo(f"\n  CUDA: Available")
            click.echo(f"    Devices: {torch.cuda.device_count()}")
        else:
            click.echo(f"\n  CUDA: Not available")
    except ImportError:
        click.echo(f"\n  CUDA: Unknown (torch not installed)")
    
    click.echo()
