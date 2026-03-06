#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import argparse
import sys
from pathlib import Path

try:
    import gromacs
except ImportError:
    print("Error: GromacsWrapper is not installed or not in the Python path.")
    print("Please install it, e.g., using 'pip install GromacsWrapper'")
    sys.exit(1)

def run_command(command, outfile=None, infile_content=None):
    """
    Helper function to run an external command robustly.

    Args:
        command (list): The command to execute as a list of strings.
        outfile (str, optional): Path to a file to redirect stdout to.
        infile_content (str, optional): String to be passed to stdin.
    """
    print(f"Running: {' '.join(command)}")
    try:
        if outfile:
            # If an output file is specified, open it and redirect stdout there.
            with open(outfile, 'w') as f_out:
                result = subprocess.run(
                    command,
                    check=True,
                    text=True,
                    stdout=f_out,
                    stderr=subprocess.PIPE,
                    input=infile_content
                )
        else:
            # If no output file, run the command and let its output go to the console.
            result = subprocess.run(
                command,
                check=True,
                text=True,
                stderr=subprocess.PIPE,
                input=infile_content
            )
        return result

    except FileNotFoundError:
        print(f"Error: Command '{command[0]}' not found. Please ensure it is in your PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error running {' '.join(command)}:")
        # Capture and print stdout/stderr on error for better debugging
        if e.stdout:
            print("--- STDOUT ---")
            print(e.stdout)
        if e.stderr:
            print("--- STDERR ---")
            print(e.stderr)
        sys.exit(1)

def main():
    """
    Main function to build a peptide structure and generate its GROMACS topology.
    """
    parser = argparse.ArgumentParser(
        description="Generate a peptide structure from a single-letter amino acid sequence "
                    "and prepare a GROMACS topology for it with uncapped (NH3+/COO-) termini."
    )
    parser.add_argument(
        "sequence",
        type=str,
        help="Amino acid sequence in single-letter code (e.g., 'GATTACA')."
    )
    parser.add_argument(
        "-ss", "--secondary-structure",
        type=str,
        default="l",
        help="Secondary structure for PCcli, 'l' for linear/loop (default)."
    )
    parser.add_argument(
        "--name",
        type=str,
        help="Protein name to use for output files (default: use sequence)."
    )
    args = parser.parse_args()

    sequence = args.sequence
    # Use protein name for files if provided, otherwise default to the sequence
    protein_name = args.name if args.name else sequence
    
    # Use Path objects for modern and robust file system interactions
    work_dir = Path(__file__).parent.resolve()
    os.chdir(work_dir)
    
    # Always use a local 'out' directory for intermediate files
    out_dir = work_dir / "out"

    print(f"--- Starting peptide preparation for sequence: {sequence} ---")
    print(f"Working directory: {work_dir}")

    # 1. Cleanup and Setup
    if out_dir.exists():
        print(f"Removing existing output directory: {out_dir}")
        shutil.rmtree(out_dir)
    out_dir.mkdir()

    # 2. Generate PDB from sequence using PCcli
    pdb_from_pccli = work_dir / f"{protein_name}.pdb"
    print(f"\nStep 1: Generating PDB structure using PCcli -> {pdb_from_pccli}")
    run_command(["PCcli", "-s", sequence, "-o", str(pdb_from_pccli), "-ss", args.secondary_structure])

    # 3. Pre-process the generated PDB file to remove unwanted lines
    pdb_processed = work_dir / f"{protein_name}_processed.pdb"
    print(f"Step 2: Pre-processing PDB file -> {pdb_processed}")
    with pdb_from_pccli.open('r') as infile, pdb_processed.open('w') as outfile:
        for line in infile:
            if "OXT" not in line and not line.startswith("HETATM") and 'USER' not in line:
                outfile.write(line)

    # 4. Run gmx pdb2gmx using GromacsWrapper for robust execution
    print("\nStep 3: Running gmx pdb2gmx to generate initial topology")
    try:
        gromacs.pdb2gmx(
            f=str(pdb_processed),
            o=f"{protein_name}-pace.pdb",
            p="draft.top",
            ff="pace-new",
            ter=True,
            ignh=True,
            input=('1', '0', '0')  # Amber ff, N-term: NH3+, C-term: COO-
        )
    except Exception as e:
        print(f"FATAL ERROR: gmx pdb2gmx failed. Please ensure GROMACS is correctly installed,")
        print(f"configured, and the 'pace-new' force field is available.")
        print(f"GromacsWrapper error: {e}")
        sys.exit(1)

    # 5. Process topology with a series of helper scripts
    print("\nStep 4: Processing topology with helper scripts...")
    
    # Run first set of scripts
    run_command([sys.executable, "change_resid.py", f"{protein_name}-pace.pdb"], outfile=f"{protein_name}-pace-resid.pdb")
    run_command([sys.executable, "genPair.py", "draft.top"], outfile=f"{protein_name}-pace.patch")

    # Read intermediate result needed for the next script
    with open("res.temp", "r") as f:
        residue_count = f.read().strip()
    
    # Run second set of scripts
    run_command([sys.executable, "insert_param.py", f"{protein_name}-pace.patch", "draft.top"], outfile=f"{protein_name}-pace.top")
    run_command([sys.executable, "change_molecule_name.py", f"{protein_name}-pace.top", protein_name])
    run_command([sys.executable, "insert_posre.py", f"{protein_name}-pace.top", "posre.itp"])
    run_command([sys.executable, "C-N-ter.py", f"{protein_name}-pace-resid.pdb", residue_count, f"{protein_name}-pace.top", "both"], outfile=f"{protein_name}-final.top")
    # Run ring-refine scripts
    run_command([sys.executable, "ring_refine.py", f"{protein_name}-final.top"], outfile=f"{protein_name}-pace-ring.top")
    # 6. Organize final output files into the output directory
    print("\nStep 5: Organizing final files...")
    shutil.move(f"{protein_name}-pace-resid.pdb", out_dir / "system.pdb")
    # Always output topol.top, the workflow generator will rename it
    shutil.move(f"{protein_name}-pace-ring.top", out_dir / "topol.top")

    # 7. Final cleanup of all intermediate and temporary files
    print("\nStep 6: Cleaning up temporary files...")
    temp_files = [
        f"{protein_name}.pdb", f"{protein_name}_processed.pdb", f"{protein_name}-pace.pdb",
        "draft.top", f"{protein_name}-pace.patch", "res.temp",
        f"{protein_name}-pace.top", "posre.itp", "mdout.mdp", f"{protein_name}-final.top", f"{protein_name}-pace-ring.top"
    ]
    for f_name in temp_files:
        f_path = Path(f_name)
        if f_path.exists():
            try:
                f_path.unlink()
            except OSError as e:
                print(f"Error cleaning up file {f_path}: {e}")

    print(f"\n--- Preparation complete! ---")
    print(f"Final structure and topology files are located in: '{out_dir}'")

if __name__ == "__main__":
    main()
