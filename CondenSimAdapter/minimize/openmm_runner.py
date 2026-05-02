"""
OpenMM in-process energy minimization runner.

This module provides run_minimization() for use directly from the main process
(imported by minimizer.py).  Running in-process avoids subprocess environment
issues that can prevent CUDA/GPU discovery.

worker.py is kept as a thin CLI wrapper around this function for standalone use.
"""

from __future__ import annotations

import gc
import os
import random
import shutil
from pathlib import Path

try:
    import openmm.unit as unit
    import openmm as mm
    from openmm.app import GromacsGroFile, Simulation, PDBFile, forcefield as ff
    from openmm import LangevinIntegrator, Platform
except ImportError:
    import simtk.unit as unit
    import simtk.openmm as mm
    from simtk.openmm.app import GromacsGroFile, Simulation, PDBFile, forcefield as ff
    from simtk.openmm import LangevinIntegrator, Platform

from .softcore import (
    GromacsTopFileWithSoftcore,
    NONBONDED_GAUSSIAN,
    NONBONDED_SOFTCORE,
    NONBONDED_STANDARD,
    IMPLICIT_GBN2,
    IMPLICIT_OBC2,
)
from ..src.device_utils import resolve_openmm_platform


def _get_lambda_values(level: str):
    """Return lambda schedule for a given optimization level."""
    configs = {
        "high":   [0.65, 0.75, 0.85, 0.95],
        "medium": [0.75, 0.85, 0.95],
        "low":    [0.85, 0.95],
    }
    level = level.lower()
    if level not in configs:
        raise ValueError(f"Unknown optimization level: {level!r}. Valid: {list(configs)}")
    return configs[level]


def _select_platform(device: str, gpu_id: int):
    """
    Resolve and probe the best available OpenMM platform.

    Fallback chain:  CUDA → OpenCL → CPU

    The probe context is created WITHOUT per-device properties to avoid
    CUDA_ERROR_OPERATING_SYSTEM false-negatives on some Linux configurations.
    The real simulation contexts are then created WITH the full properties dict.

    Returns
    -------
    platform : Platform
    properties : dict
    platform_name : str
    """
    platform_name, properties = resolve_openmm_platform(device, gpu_id, precision="mixed")

    if platform_name == "CUDA":
        fallback_chain = [
            ("CUDA",   properties),
            ("OpenCL", {"Precision": "mixed", "DeviceIndex": str(gpu_id)}),
            ("CPU",    {}),
        ]
    elif platform_name == "OPENCL":
        fallback_chain = [
            ("OpenCL", properties),
            ("CPU",    {}),
        ]
    else:
        fallback_chain = [(platform_name, properties)]

    def _probe(name):
        p = Platform.getPlatformByName(name)
        if name not in ("CPU", "Reference"):
            _s = mm.System(); _s.addParticle(1.0)
            _i = mm.VerletIntegrator(0.001)
            _c = mm.Context(_s, _i, p)   # no properties — avoids false failures
            del _c, _i, _s
        return p

    for try_name, try_props in fallback_chain:
        try:
            p = _probe(try_name)
            print(f"  [platform] Using {try_name}", flush=True)
            return p, try_props, try_name
        except Exception as exc:
            print(f"  [warn] {try_name} unavailable ({exc}), trying next ...", flush=True)

    raise RuntimeError("No usable OpenMM platform found (tried CUDA, OpenCL, CPU).")


def run_minimization(
    input_gro: str,
    input_top: str,
    output_dir: str,
    device: str = "cuda",
    gpu_id: int = 0,
    max_iterations: int = 5000,
    optimization_level: str = "medium",
    ff_type: str = "amber",
    ff_name: str = "amber99sb-ildn",
    gb_model: str = "OBC2",
    salt_conc: float = 0.15,
    cutoff: float = 2.0,
    tolerance: float = 100.0,
) -> dict:
    """
    Run multi-stage softcore energy minimization in the current process.

    Stages
    ------
    1. Gaussian repulsion  — resolves severe clashes
    2. Softcore (N steps)  — progressive introduction of full LJ/Coulomb
    3. Standard            — final minimization with implicit solvent (AMBER)
                             or high-lambda softcore (CHARMM)

    Parameters
    ----------
    input_gro        : path to GROMACS .gro structure
    input_top        : path to GROMACS .top topology
    output_dir       : directory for output PDB files
    device           : "cuda" | "opencl" | "cpu"
    gpu_id           : GPU device index
    max_iterations   : max minimisation iterations per stage
    optimization_level: "high" | "medium" | "low"
    ff_type          : "amber" | "charmm"
    ff_name          : GROMACS pdb2gmx force-field name (informational)
    gb_model         : "GBn2" | "OBC2" (implicit solvent, AMBER only)
    salt_conc        : salt concentration in M
    cutoff           : nonbonded cutoff in nm
    tolerance        : minimisation energy tolerance in kJ/(mol·nm)

    Returns
    -------
    dict with keys "step_info" and "total_steps"
    """
    input_gro_path  = Path(input_gro).resolve()
    input_top_path  = Path(input_top).resolve()
    output_path     = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    # Copy inputs into output dir so topology relative #includes resolve
    conf_gro     = output_path / "conf.gro"
    minimize_top = output_path / "topol.top"
    if str(input_gro_path) != str(conf_gro):
        shutil.copy2(str(input_gro_path), str(conf_gro))
    if str(input_top_path) != str(minimize_top):
        shutil.copy2(str(input_top_path), str(minimize_top))

    original_cwd = os.getcwd()
    os.chdir(output_path)

    try:
        # ── Load structure ────────────────────────────────────────────────────
        conf = GromacsGroFile("conf.gro")
        gro_positions = conf.getPositions()
        box_vectors   = conf.getPeriodicBoxVectors()

        # includeDir = parent of minimize/ (where the FF folder lives, e.g. a99SBdisp.ff/)
        top = GromacsTopFileWithSoftcore(
            "topol.top",
            periodicBoxVectors=box_vectors,
            includeDir=str(output_path.parent),
            forcefield_type=ff_type.upper(),
        )

        # ── Platform selection ───────────────────────────────────────────────
        platform, properties, platform_name = _select_platform(device, gpu_id)

        current_positions = gro_positions
        gb_constant = IMPLICIT_GBN2 if gb_model.upper() == "GBN2" else IMPLICIT_OBC2

        # ── Stage 1: Gaussian repulsion ───────────────────────────────────────
        system_gauss = top.createSystem(
            nonbondedCutoff=cutoff * unit.nanometer,
            nonbondedMethod=ff.CutoffPeriodic,
            nonbonded_type=NONBONDED_GAUSSIAN,
            add_implicit_solvent=False,
        )
        integ = LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picosecond)
        integ.setRandomNumberSeed(random.randint(0, 2**31 - 1))
        sim = Simulation(top.topology, system_gauss, integ, platform, properties)
        sim.context.setPositions(current_positions)
        sim.context.setPeriodicBoxVectors(*box_vectors)
        sim.minimizeEnergy(maxIterations=max_iterations, tolerance=50.0)

        state = sim.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
        PDBFile.writeFile(top.topology, state.getPositions(), open("step1_gaussian.pdb", "w"))
        print("  Generated step1_gaussian.pdb", flush=True)
        current_positions = state.getPositions()

        del system_gauss, sim, state, integ
        gc.collect()

        # ── Stage 2: Softcore (progressive lambda) ────────────────────────────
        lambda_values = _get_lambda_values(optimization_level)
        add_implicit   = ff_type.upper() == "AMBER"

        for step_num, lam in enumerate(lambda_values, 1):
            system_sc = top.createSystem(
                nonbondedCutoff=cutoff * unit.nanometer,
                nonbondedMethod=ff.CutoffPeriodic,
                nonbonded_type=NONBONDED_SOFTCORE,
                add_implicit_solvent=add_implicit,
                gb_model=gb_constant if add_implicit else None,
                salt_conc=salt_conc,
                soft_lambda=lam,
            )
            integ = LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picosecond)
            integ.setRandomNumberSeed(random.randint(0, 2**31 - 1))
            sim = Simulation(top.topology, system_sc, integ, platform, properties)
            sim.context.setPositions(current_positions)
            sim.context.setPeriodicBoxVectors(*box_vectors)
            sim.minimizeEnergy(maxIterations=max_iterations, tolerance=tolerance)

            state = sim.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
            current_positions = state.getPositions()
            out_pdb = f"step2_softcore_{step_num}.pdb"
            PDBFile.writeFile(top.topology, current_positions, open(out_pdb, "w"))
            print(f"  Generated {out_pdb} (lambda={lam})", flush=True)

            del system_sc, sim, state, integ
            gc.collect()

        # ── Stage 3: Final minimization ───────────────────────────────────────
        if ff_type.upper() == "CHARMM":
            # CHARMM: GB forces are unreliable — use very-high-lambda softcore
            final_lambda = 0.99999
            print(f"  Using softcore final minimize (lambda={final_lambda}) for CHARMM", flush=True)
            system_final = top.createSystem(
                nonbondedCutoff=cutoff * unit.nanometer,
                nonbondedMethod=ff.CutoffPeriodic,
                nonbonded_type=NONBONDED_SOFTCORE,
                add_implicit_solvent=False,
                soft_lambda=final_lambda,
            )
        else:
            # AMBER: standard potential with implicit solvent
            system_final = top.createSystem(
                nonbondedCutoff=cutoff * unit.nanometer,
                nonbondedMethod=ff.CutoffPeriodic,
                nonbonded_type=NONBONDED_STANDARD,
                add_implicit_solvent=True,
                gb_model=gb_constant,
                salt_conc=salt_conc,
            )

        integ = LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picosecond)
        integ.setRandomNumberSeed(random.randint(0, 2**31 - 1))
        sim = Simulation(top.topology, system_final, integ, platform, properties)
        sim.context.setPositions(current_positions)
        sim.context.setPeriodicBoxVectors(*box_vectors)
        sim.minimizeEnergy(maxIterations=max_iterations, tolerance=tolerance)

        state = sim.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
        PDBFile.writeFile(top.topology, state.getPositions(), open("final.pdb", "w"))
        print("  Generated final.pdb", flush=True)

        del system_final, sim, state, integ
        gc.collect()

        total_steps = 1 + len(lambda_values) + 1
        step_info   = f"{optimization_level.capitalize()} ({total_steps} steps)"
        print(f"  Optimization: {step_info}", flush=True)
        print("\n  Done!", flush=True)

        return {"step_info": step_info, "total_steps": total_steps}

    finally:
        os.chdir(original_cwd)
