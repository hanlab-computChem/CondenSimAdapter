#!/usr/bin/env python3
"""
Prepare PACE-ASM topology from amino acid sequence.

This script generates a peptide structure and prepares its GROMACS topology
for PACE-ASM force field. Output format is compatible with prepare_peptide.py
(PACE-NEW), producing:
    - out/system.pdb
    - out/topol.top
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


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
        if e.stdout:
            print("--- STDOUT ---")
            print(e.stdout)
        if e.stderr:
            print("--- STDERR ---")
            print(e.stderr)
        sys.exit(1)


def main():
    """
    Main function to build a peptide structure and generate PACE-ASM topology.
    """
    parser = argparse.ArgumentParser(
        description="Generate a peptide structure from a single-letter amino acid sequence "
                    "and prepare a GROMACS topology for it with PACE-ASM force field."
    )
    parser.add_argument(
        "sequence",
        type=str,
        help="Amino acid sequence in single-letter code (e.g., 'GATTACA')."
    )
    parser.add_argument(
        "--name",
        type=str,
        help="Molecule name for output files (default: use sequence)."
    )
    args = parser.parse_args()

    sequence = args.sequence
    protein_name = args.name if args.name else sequence

    # Get script directory and change to it
    script_dir = Path(__file__).parent.resolve()
    os.chdir(script_dir)

    # Output directory
    out_dir = script_dir / "out"

    print(f"--- Starting PACE-ASM preparation for sequence: {sequence} ---")
    print(f"Working directory: {script_dir}")

    # 1. Cleanup and Setup
    if out_dir.exists():
        print(f"Removing existing output directory: {out_dir}")
        shutil.rmtree(out_dir)
    out_dir.mkdir()

    # Copy necessary files to working directory
    print("\nStep 0: Setting up environment...")

    # Helper to get script subdirectory
    scripts_dir = script_dir / "scripts"

    # Copy files from scripts directory
    shutil.copy(scripts_dir / "residuetypes.dat", ".")
    shutil.copy(scripts_dir / "comb_top.py", ".")

    # Copy MDP files
    mdp_dir = script_dir / "mdp"
    if mdp_dir.exists():
        for mdp_file in mdp_dir.glob("*"):
            shutil.copy(mdp_file, ".")

    # Note: Force field files (pace-asm.ff) are already accessible from the current directory
    # since we're running in pace_asm_top_builder. No need to copy.

    # 2. Generate PDB from sequence using PCcli
    pdb_from_pccli = script_dir / f"{protein_name}.pdb"
    print(f"\nStep 1: Generating PDB structure using PCcli -> {pdb_from_pccli}")
    run_command(["PCcli", "-s", sequence, "-o", str(pdb_from_pccli), "-ss", "l"])

    # 3. Pre-process the generated PDB file
    pdb_processed = script_dir / f"{protein_name}_processed.pdb"
    print(f"Step 2: Pre-processing PDB file -> {pdb_processed}")
    with pdb_from_pccli.open('r') as infile, pdb_processed.open('w') as outfile:
        for line in infile:
            if "OXT" not in line and not line.startswith("HETATM") and 'USER' not in line:
                outfile.write(line)

    # 4. Run gmx pdb2gmx with PACE-ASM force field (using charged termini: NH3+ / COO-)
    print("\nStep 3: Running gmx pdb2gmx with PACE-ASM force field")
    print("  Using charged termini (NH3+ / COO-)...")
    try:
        import gromacs
        gromacs.pdb2gmx(
            f=str(pdb_processed),
            o=f"{protein_name}-pace.pdb",
            p="draft.top",
            ff="pace-asm",
            ter=True,
            ignh=True,
            input=('1', '0', '0')  # Force field, N-term: NH3+, C-term: COO-
        )
    except ImportError:
        print("Error: GromacsWrapper is not installed.")
        print("Falling back to direct gmx pdb2gmx call...")
        # Fallback: run gmx pdb2gmx directly
        run_command([
            "gmx", "pdb2gmx", "-f", str(pdb_processed),
            "-o", f"{protein_name}-pace.pdb",
            "-p", "draft.top",
            "-ff", "pace-asm", "-ter", "-ignh"
        ], infile_content="1\n0\n0\n")
    except Exception as e:
        print(f"FATAL ERROR: gmx pdb2gmx failed: {e}")
        sys.exit(1)

    # 5. Process topology with helper scripts
    print("\nStep 4: Processing topology with helper scripts...")

    # Change residue numbering
    run_command([sys.executable, str(scripts_dir / "change_resid.py"),
                 f"{protein_name}-pace.pdb"], outfile=f"{protein_name}-pace-resid.pdb")

    # Get residue and atom counts
    print(f"  DEBUG: Reading {protein_name}-pace-resid.pdb to get counts...")
    with open(f"{protein_name}-pace-resid.pdb", 'r') as f:
        lines = f.readlines()
        for line in reversed(lines):
            if line.startswith("ATOM"):
                count_residue = line[22:26].strip()
                count_atom = line[6:11].strip()
                print(f"  DEBUG: Last atom line: {line[:60].strip()}")
                print(f"  DEBUG: count_residue = {count_residue}, count_atom = {count_atom}")
                break

    # Ensure genPairPACE is compiled
    genpair_script = scripts_dir / "genpair" / "genPairPACE"
    genpair_dir = scripts_dir / "genpair"
    genpair_c = genpair_dir / "genPairPACE.c"

    if not genpair_script.exists() or genpair_script.stat().st_mtime < genpair_c.stat().st_mtime:
        print("  genPairPACE not compiled or source modified. Compiling now...")
        try:
            run_command(["make", "-C", str(genpair_dir)])
            print("  genPairPACE compiled successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error: Failed to compile genPairPACE: {e}")
            sys.exit(1)

    if not genpair_script.exists():
        print(f"Error: genPairPACE not found at {genpair_script} after compilation attempt")
        sys.exit(1)

    # Run genPairPACE (use mode 1 for charged/uncapped termini)
    # Mode 1: forgiving mode (for charged/uncapped peptides)
    run_command([
        str(genpair_script), count_atom, count_residue,
        f"{protein_name}-pace-resid.pdb", "1"
    ], outfile=f"{protein_name}-pace.patch")

    # Insert parameters
    run_command([sys.executable, str(scripts_dir / "insert_param.py"),
                 f"{protein_name}-pace.patch", "draft.top"], outfile=f"{protein_name}-pace.top")

    # 6. Handle terminal modifications (using charged termini)
    print("  Processing terminal modifications for charged termini...")
    print(f"  DEBUG: count_residue = {count_residue}")
    run_command([sys.executable, str(scripts_dir / "C-N-ter.py"),
                 f"{protein_name}-pace-resid.pdb", count_residue,
                 f"{protein_name}-pace.top", "both"], outfile=f"{protein_name}-pace-final.top")

    # 7. Handle position restraint file
    if os.path.exists("posre.itp"):
        shutil.move("posre.itp", f"posre_{protein_name}.itp")
        run_command([sys.executable, str(scripts_dir / "rpl_posre.py"),
                     f"{protein_name}-pace-final.top", protein_name, f"{protein_name}-pace.top"])

    # 8. Organize final output files into the output directory
    print("\nStep 5: Organizing final files...")
    shutil.move(f"{protein_name}-pace-resid.pdb", out_dir / "system.pdb")
    shutil.move(f"{protein_name}-pace.top", out_dir / "topol.top")

    # 9. Final cleanup
    print("\nStep 6: Cleaning up temporary files...")
    temp_files = [
        f"{protein_name}.pdb",
        f"{protein_name}_processed.pdb",
        f"{protein_name}-pace.pdb",
        "draft.top",
        f"{protein_name}-pace.patch",
        f"{protein_name}-pace-resid.pdb",
        f"{protein_name}-pace.top",
        f"{protein_name}-pace-final.top",
        "posre.itp",
        "residuetypes.dat",
        "comb_top.py",
        "mdout.mdp",
    ]

    # Also clean up MDP files that were copied
    if mdp_dir.exists():
        for mdp_file in mdp_dir.glob("*"):
            temp_files.append(mdp_file.name)

    # Clean up force field directory
    if os.path.exists("pace-asm.ff"):
        import shutil as sh
        try:
            sh.rmtree("pace-asm.ff")
        except:
            pass

    for f_name in temp_files:
        f_path = Path(f_name)
        if f_path.exists():
            try:
                f_path.unlink()
            except OSError as e:
                print(f"Error cleaning up file {f_path}: {e}")

    print("\n--- PACE-ASM Preparation complete! ---")
    print(f"Final structure and topology files are located in: '{out_dir}'")


if __name__ == "__main__":
    main()
