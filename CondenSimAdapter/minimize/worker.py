#!/usr/bin/env python3
"""
Minimization Worker Script (standalone CLI wrapper)

Thin CLI entry-point around openmm_runner.run_minimization().
Used when the minimization needs to be run as a standalone script
(e.g. for debugging or manual invocation).

For normal use the minimizer.py calls openmm_runner directly in-process.
"""

import sys
from pathlib import Path

# Make the minimize/ package importable when run as a standalone script
sys.path.insert(0, str(Path(__file__).parent.parent))  # CondenSimAdapter/
sys.path.insert(0, str(Path(__file__).parent))          # minimize/

from CondenSimAdapter.minimize.openmm_runner import run_minimization


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Standalone minimization worker (thin wrapper around openmm_runner)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python worker.py -i conf.gro -t topol.top -o output_dir
    python worker.py -i conf.gro -t topol.top -o output_dir --device cuda
        """,
    )
    parser.add_argument("-i", "--input-gro",  required=True,  help="Input GRO file")
    parser.add_argument("-t", "--input-top",  required=True,  help="Input topology file")
    parser.add_argument("-o", "--output",     required=True,  help="Output directory")
    parser.add_argument("-d", "--device",     default="cuda",
                        choices=["cuda", "cpu", "opencl"],
                        help="Compute device (default: cuda)")
    parser.add_argument("-g", "--gpu-id",     type=int, default=0,
                        help="GPU device index (default: 0)")
    parser.add_argument("--iter",             type=int, default=5000,
                        help="Max iterations per stage (default: 5000)")
    parser.add_argument("--tolerance",        type=float, default=100.0,
                        help="Energy tolerance kJ/(mol·nm) (default: 100.0)")
    parser.add_argument("-l", "--level",      default="medium",
                        choices=["high", "medium", "low"],
                        help="Optimization level (default: medium)")
    parser.add_argument("--ff-type",          default="amber",
                        choices=["amber", "charmm"],
                        help="Force field family (default: amber)")
    parser.add_argument("--ff-name",          default="amber99sb-ildn",
                        help="pdb2gmx force-field name (informational)")
    parser.add_argument("--gb-model",         default="GBn2",
                        choices=["GBn2", "OBC2"],
                        help="Implicit solvent model (default: GBn2)")
    parser.add_argument("--salt-conc",        type=float, default=0.15,
                        help="Salt concentration in M (default: 0.15)")
    parser.add_argument("--cutoff",           type=float, default=2.0,
                        help="Nonbonded cutoff in nm (default: 2.0)")

    args = parser.parse_args()

    run_minimization(
        input_gro=args.input_gro,
        input_top=args.input_top,
        output_dir=args.output,
        device=args.device,
        gpu_id=args.gpu_id,
        max_iterations=args.iter,
        optimization_level=args.level,
        ff_type=args.ff_type,
        ff_name=args.ff_name,
        gb_model=args.gb_model,
        salt_conc=args.salt_conc,
        cutoff=args.cutoff,
        tolerance=args.tolerance,
    )
