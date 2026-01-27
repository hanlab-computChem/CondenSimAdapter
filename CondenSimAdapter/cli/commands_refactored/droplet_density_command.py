#!/usr/bin/env python3
"""
Droplet Density Command

Estimates protein density in a droplet geometry.
"""

import math
import sys
from pathlib import Path

import click

from ...src import CGSimulationConfig


@click.command('droplet-density', context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--input-file', '-f',
    type=click.Path(exists=True),
    required=True,
    help='Configuration YAML file'
)
@click.option(
    '--radius', '-r',
    type=float,
    required=True,
    help='Droplet radius in nm'
)
@click.option(
    '--nmol', '-n',
    type=str,
    default=None,
    required=False,
    help='Number of molecules for each component (space-separated, e.g., -n "10 20")'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Show detailed calculation'
)
@click.argument('extra_nmol', nargs=-1, type=str, required=False)
def droplet_density_command(input_file: str, radius: float, nmol: str, verbose: bool, extra_nmol: tuple):
    """\b
    Estimate protein density in a droplet geometry.

    \b
    Calculates the protein concentration (mg/mL) based on:
        - Configuration YAML (components, sequences, residue counts)
        - Provided droplet radius
        - Optional: number of molecules per component (-n flag)

    \b
    Use -n to specify molecule counts and calculate achievable density.
    Example: adapter droplet-density -f config.yaml -r 15 -n 100 200

    \b
    Warnings are issued if density is below 300 mg/mL or above 800 mg/mL.

    \b
    Examples:
        adapter droplet-density -f config.yaml -r 15
        adapter droplet-density -f config.yaml -r 20 -n 10 20
        adapter droplet-density -f config.yaml -r 20 --verbose
    """
    # Combine -n value with extra positional arguments
    nmol_values = None
    if nmol:
        try:
            nmol_values = tuple(int(v) for v in nmol.split())
        except ValueError:
            click.echo(f"  ✗ Error: Invalid -n value: {nmol}")
            sys.exit(1)

    if extra_nmol:
        extra_values = tuple(int(v) for v in extra_nmol)
        if nmol_values:
            nmol_values = nmol_values + extra_values
        else:
            nmol_values = extra_values

    click.echo(f"\n{'=' * 60}")
    click.echo(f"Droplet Density Estimation")
    click.echo(f"{'=' * 60}\n")

    # Load configuration
    try:
        config = CGSimulationConfig.from_yaml(input_file)
    except Exception as e:
        click.echo(f"  ✗ Failed to load configuration: {e}")
        sys.exit(1)

    # Validate -n flag if provided
    if nmol_values:
        n_components = len(config.components)
        n_provided = len(nmol_values)

        if n_provided != n_components:
            click.echo(f"  ✗ Error: -n flag requires {n_components} values (one per component)")
            click.echo(f"    Found {n_components} component(s) in configuration:")
            for i, comp in enumerate(config.components):
                click.echo(f"      {i+1}. {comp.name} ({comp.type.value})")
            click.echo(f"\n    You provided {n_provided} value(s): {nmol_values}")
            click.echo(f"    Please provide exactly {n_components} integer(s).")
            sys.exit(1)

        # Check for negative values
        if any(n <= 0 for n in nmol_values):
            click.echo(f"  ✗ Error: All -n values must be positive integers")
            sys.exit(1)

    # Constants
    AVG_RESIDUE_MASS = 110.0  # Da (g/mol) per residue (fallback)
    AVOGADRO = 6.022e23  # molecules/mol
    WATER_MASS = 18.015  # Da (released during peptide bond formation)

    # Amino acid molecular weights (Da) - monoisotopic mass
    AA_WEIGHTS = {
        'A': 89.09, 'C': 121.15, 'D': 133.10, 'E': 147.13, 'F': 165.19,
        'G': 75.07, 'H': 155.16, 'I': 131.18, 'K': 146.19, 'L': 131.18,
        'M': 149.21, 'N': 132.12, 'P': 115.13, 'Q': 146.15, 'R': 174.20,
        'S': 105.09, 'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19,
        # Non-standard or ambiguous
        'U': 168.05,  # Selenocysteine
        'O': 255.31,  # Pyrrolysine
        'B': 132.61,  # Asx (average of D and N)
        'Z': 146.64,  # Glx (average of E and Q)
        'X': 110.0,   # Unknown (use average)
    }

    def calculate_protein_mass(sequence: str) -> float:
        """
        Calculate exact protein mass from sequence.

        Args:
            sequence: Amino acid sequence (single letter code)

        Returns:
            Molecular weight in Da
        """
        sequence = sequence.upper().strip()
        if not sequence:
            return 0.0

        # Sum residue weights
        total = sum(AA_WEIGHTS.get(aa, AVG_RESIDUE_MASS) for aa in sequence)

        # Subtract water molecules lost in peptide bond formation
        # For n residues, there are (n-1) peptide bonds
        peptide_bonds = len(sequence) - 1
        total -= peptide_bonds * WATER_MASS

        return total

    # Calculate droplet volume (nm³)
    volume_nm3 = (4.0 / 3.0) * math.pi * (radius ** 3)

    # Convert to liters (1 nm³ = 1e-24 L)
    volume_L = volume_nm3 * 1e-24

    # Calculate total mass and molecules
    total_mass_Da = 0.0
    total_molecules = 0
    component_details = []
    exact_mass_count = 0
    estimated_mass_count = 0
    nmol_source = []  # Track where nmol values come from

    for idx, comp in enumerate(config.components):
        # Use user-provided nmol if -n flag is set, otherwise use config value
        if nmol_values:
            current_nmol = nmol_values[idx]
            nmol_source.append('user')
        else:
            current_nmol = comp.nmol
            nmol_source.append('config')

        # Try to get actual sequence
        sequence = None
        nres = 0

        # Priority 1: comp.seq
        if comp.seq:
            sequence = comp.seq
            nres = len(sequence)

        # Priority 2: FASTA file
        if not sequence and comp.ffasta:
            try:
                fasta_path = Path(input_file).parent / comp.ffasta
                if fasta_path.exists():
                    with open(fasta_path, 'r') as f:
                        lines = f.readlines()
                        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
                        nres = len(sequence)
            except Exception:
                pass

        # Priority 3: PDB file
        if not sequence and comp.fpdb:
            try:
                pdb_path = Path(input_file).parent / comp.fpdb
                if pdb_path.exists():
                    from Bio.PDB import PDBParser
                    from Bio.SeqUtils import seq1
                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure('protein', str(pdb_path))
                    # Try to extract sequence from PDB
                    residues = list(structure.get_residues())
                    nres = len(residues)
                    try:
                        # Convert 3-letter codes to 1-letter
                        sequence = ''.join(seq1(res.get_resname()) for res in residues)
                    except Exception:
                        # If conversion fails, just count residues
                        sequence = None
            except Exception:
                pass

        # Priority 4: comp.nres (if set and we don't have sequence yet)
        if nres == 0 and comp.nres > 0:
            nres = comp.nres

        if nres == 0:
            click.echo(f"  ✗ Warning: Could not determine residue count for component '{comp.name}'")
            click.echo(f"    Please ensure sequence or structure files are accessible.")
            continue

        # Calculate mass for this component
        if sequence:
            # Use exact mass calculation from sequence
            mass_per_molecule = calculate_protein_mass(sequence)
            comp_mass_Da = mass_per_molecule * current_nmol
            exact_mass_count += 1
            mass_method = 'exact'
        else:
            # Use average residue mass
            mass_per_molecule = nres * AVG_RESIDUE_MASS
            comp_mass_Da = mass_per_molecule * current_nmol
            estimated_mass_count += 1
            mass_method = 'estimated'

        total_mass_Da += comp_mass_Da
        total_molecules += current_nmol

        component_details.append({
            'name': comp.name,
            'type': comp.type.value,
            'nmol': current_nmol,
            'nmol_source': nmol_source[idx],
            'nres': nres,
            'mass_per_mol': mass_per_molecule,
            'total_mass': comp_mass_Da,
            'mass_method': mass_method
        })

    if total_mass_Da == 0:
        click.echo(f"  ✗ Error: No valid components found or masses are zero")
        sys.exit(1)

    # Convert to grams
    total_mass_g = total_mass_Da / AVOGADRO

    # Calculate density (mg/mL)
    # Note: g/L = mg/mL (1 g/L = 1000 mg / 1000 mL = 1 mg/mL)
    density_mgmL = total_mass_g / volume_L  # g/L = mg/mL

    # Display results
    click.echo(f"  Input:")
    click.echo(f"    Configuration: {input_file}")
    click.echo(f"    Radius: {radius:.2f} nm")
    click.echo(f"    Volume: {volume_nm3:.2f} nm³ ({volume_L:.6e} L)")

    click.echo(f"\n  Composition:")
    click.echo(f"    Total components: {len(component_details)}")
    click.echo(f"    Total molecules: {total_molecules}")
    click.echo(f"    Total mass: {total_mass_Da:.2f} Da ({total_mass_g:.6e} g)")

    # Show nmol source
    if nmol_values:
        click.echo(f"    Molecule counts: User-provided via -n flag")
        for idx, comp_info in enumerate(component_details):
            click.echo(f"      {idx+1}. {comp_info['name']}: {comp_info['nmol']} molecules")
    else:
        click.echo(f"    Molecule counts: From configuration file")

    # Show mass calculation method
    if exact_mass_count > 0 and estimated_mass_count == 0:
        click.echo(f"    Mass calculation: Exact (from sequences)")
    elif exact_mass_count == 0 and estimated_mass_count > 0:
        click.echo(f"    Mass calculation: Estimated (avg {AVG_RESIDUE_MASS} Da/residue)")
    else:
        click.echo(f"    Mass calculation: Mixed ({exact_mass_count} exact, {estimated_mass_count} estimated)")

    if verbose:
        click.echo(f"\n  Component details:")
        for comp_info in component_details:
            nmol_str = "✓ user" if comp_info['nmol_source'] == 'user' else "≈ config"
            method_str = "✓ exact" if comp_info['mass_method'] == 'exact' else "≈ estimated"
            click.echo(f"    - {comp_info['name']} ({comp_info['type'].upper()}):")
            click.echo(f"        Molecules: {comp_info['nmol']} ({nmol_str})")
            click.echo(f"        Residues/molecule: {comp_info['nres']}")
            click.echo(f"        Mass/molecule: {comp_info['mass_per_mol']:.2f} Da ({method_str})")
            click.echo(f"        Total mass: {comp_info['total_mass']:.2f} Da")

    click.echo(f"\n  Density:")
    click.echo(f"    {density_mgmL:.1f} mg/mL (g/L)")

    # Check warnings
    if density_mgmL < 300:
        click.echo(f"\n  ⚠ WARNING: Density ({density_mgmL:.1f} mg/mL) is below 300 mg/mL")
        click.echo(f"    This may be too dilute for droplet formation.")
    elif density_mgmL > 800:
        click.echo(f"\n  ⚠ WARNING: Density ({density_mgmL:.1f} mg/mL) is above 800 mg/mL")
        click.echo(f"    This may be too concentrated and could cause simulation issues.")
    else:
        click.echo(f"    ✓ Density is within recommended range (300-800 mg/mL)")

    click.echo()
