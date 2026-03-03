"""
Unified CG simulation driver.

Replaces the five separate run_*() methods that previously lived in
the 3358-line cg.py.  All force fields go through the same pipeline:

  build coords -> build topology -> assemble forces -> run MD
"""

from __future__ import annotations

import logging
import os
import shutil
import time
from pathlib import Path
from typing import List, Optional

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from mdtraj.reporters import XTCReporter
from tqdm import tqdm

from .config import CGConfig, SimulationResult, TopologyType
from .molecule import build_all_chains
from .topology import build_topology, get_masses
from .forcefield import create_forcefield
from .forcefield.base import CGForceField
from .forcefield.cocomo import CocomoFF

log = logging.getLogger(__name__)


class CGSimulation:
    """
    Unified coarse-grained simulation engine.

    All CG force fields (CALVADOS, HPS, COCOMO, Mpipi) are handled
    identically -- the only difference is the force field object used to
    populate the OpenMM System.
    """

    def __init__(self, config: CGConfig):
        self.config = config
        self.ff: CGForceField = create_forcefield(config.force_field)
        self._result: Optional[SimulationResult] = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        output_dir: str,
        gpu_id: Optional[int] = None,
        overwrite: bool = False,
    ) -> SimulationResult:
        """
        Build and run the CG simulation.

        Args:
            output_dir: Directory to write output files into.
            gpu_id:     GPU device index; None uses config default.
            overwrite:  If False, raises an error if output_dir exists.

        Returns:
            SimulationResult with paths to trajectory, final PDB, and log.
        """
        out = Path(output_dir)
        if out.exists() and not overwrite:
            raise FileExistsError(
                f"Output directory '{out}' already exists. "
                "Pass overwrite=True to proceed."
            )
        out.mkdir(parents=True, exist_ok=True)

        gpu_id = gpu_id if gpu_id is not None else self.config.simulation.gpu_id
        t0 = time.perf_counter()

        try:
            result = self._run_pipeline(out, gpu_id)
        except Exception as exc:
            log.exception("CG simulation failed")
            result = SimulationResult(
                success=False,
                output_dir=str(out),
                error=str(exc),
                elapsed_seconds=time.perf_counter() - t0,
            )

        result.elapsed_seconds = time.perf_counter() - t0
        self._result = result
        return result

    # ------------------------------------------------------------------
    # Internal pipeline
    # ------------------------------------------------------------------

    def _run_pipeline(self, out: Path, gpu_id: int) -> SimulationResult:
        cfg = self.config

        # 1. Build CG coordinates + chain metadata
        log.info("Building initial coordinates ...")
        positions, chain_meta = build_all_chains(cfg)

        # 2. Build OpenMM Topology
        topology, pos_qty = build_topology(chain_meta, positions, cfg.box)

        # 3. Create OpenMM System
        log.info("Assembling force field ...")
        system = self._create_system(topology, chain_meta, positions, cfg)

        # 4. Write initial structure
        top_pdb = str(out / "top.pdb")
        _save_pdb(topology, pos_qty, top_pdb, cfg.box)

        # 5. Create Simulation object
        sim_params = cfg.simulation
        integrator = mm.LangevinMiddleIntegrator(
            cfg.temperature * unit.kelvin,
            sim_params.friction / unit.picosecond,
            sim_params.dt * unit.picosecond,
        )

        openmm_platform, platform_props = _resolve_platform(
            sim_params.platform, gpu_id
        )
        simulation = app.Simulation(
            topology, system, integrator,
            openmm_platform, platform_props
        )
        simulation.context.setPositions(pos_qty)

        # 6. Energy minimise before production
        log.info("Energy minimisation ...")
        simulation.minimizeEnergy(maxIterations=2000)

        # 7. Production MD
        xtc_file = str(out / f"{cfg.system_name}.xtc")
        log_file = str(out / f"{cfg.system_name}.log")
        chk_file = str(out / "checkpoint.chk")

        simulation.reporters.append(XTCReporter(xtc_file, sim_params.wfreq))
        simulation.reporters.append(
            app.StateDataReporter(
                log_file,
                sim_params.log_freq,
                step=True,
                potentialEnergy=True,
                temperature=True,
                remainingTime=True,
                totalSteps=sim_params.steps,
                speed=True,
            )
        )

        log.info(f"Starting MD: {sim_params.steps:,} steps ...")
        batch = max(sim_params.wfreq, 1000)
        n_batches = sim_params.steps // batch

        with tqdm(total=sim_params.steps, unit="steps", smoothing=0.05) as pbar:
            for _ in range(n_batches):
                simulation.step(batch)
                simulation.saveCheckpoint(chk_file)
                pbar.update(batch)

        # 8. Save final structure
        final_pdb = str(out / "final.pdb")
        _save_final_pdb(simulation, topology, cfg.box, final_pdb)
        log.info(f"Simulation complete. Final PDB: {final_pdb}")

        return SimulationResult(
            success=True,
            output_dir=str(out),
            final_pdb=final_pdb,
            trajectory=xtc_file,
            log_file=log_file,
        )

    # ------------------------------------------------------------------
    # System assembly
    # ------------------------------------------------------------------

    def _create_system(
        self,
        topology: app.Topology,
        chain_meta: List[dict],
        positions: np.ndarray,
        config: CGConfig,
    ) -> mm.System:
        system = mm.System()

        # 1. Add particles (masses)
        self.ff.add_masses(system, chain_meta)

        # 2. Periodic box
        bx, by, bz = config.box
        system.setDefaultPeriodicBoxVectors(
            mm.Vec3(bx, 0, 0) * unit.nanometer,
            mm.Vec3(0, by, 0) * unit.nanometer,
            mm.Vec3(0, 0, bz) * unit.nanometer,
        )

        # 3. Backbone harmonic bonds
        system.addForce(self.ff.build_harmonic_bonds(topology))

        # 4. COCOMO-specific angle force
        if isinstance(self.ff, CocomoFF):
            system.addForce(self.ff.build_angle_force(topology))

        # 5. ENM for folded domains
        has_domains = any(meta["folded_domains"] for meta in chain_meta)
        if has_domains:
            enm = self.ff.build_enm_bonds(positions, chain_meta,
                                           restraint_type="harmonic")
            if enm is not None:
                system.addForce(enm)

        # 6. Non-bonded forces (force field specific)
        if isinstance(self.ff, CocomoFF):
            nb_forces = self.ff.create_nonbonded_forces(
                topology, chain_meta,
                config.temperature, config.ionic_strength,
                positions=positions,
            )
        else:
            nb_forces = self.ff.create_nonbonded_forces(
                topology, chain_meta,
                config.temperature, config.ionic_strength,
            )
        for f in nb_forces:
            system.addForce(f)

        # 7. Droplet confinement
        if config.topology == TopologyType.DROPLET:
            r = config.droplet_radius or (min(config.box) * 0.4)
            drop_force = self.ff.build_droplet_force(topology, positions, r, config.droplet_k)
            system.addForce(drop_force)

        # 8. CMMotionRemover
        system.addForce(mm.CMMotionRemover(1000))

        return system


# ---------------------------------------------------------------------------
# Platform selection
# ---------------------------------------------------------------------------

def _resolve_platform(
    platform_name: str,
    gpu_id: int,
):
    """
    Select the best available OpenMM platform.

    Falls back to CPU if the requested platform is unavailable.
    """
    preferred = platform_name.upper()

    if gpu_id >= 0 and preferred in ("CUDA", "OPENCL"):
        try:
            plat = mm.Platform.getPlatformByName(preferred.capitalize()
                                                 if preferred == "Opencl"
                                                 else preferred)
            props = {"DeviceIndex": str(gpu_id), "Precision": "mixed"}
            return plat, props
        except Exception:
            log.warning(f"{preferred} not available -- falling back to CPU.")

    try:
        plat = mm.Platform.getPlatformByName("CPU")
    except Exception:
        plat = mm.Platform.getPlatformByName("Reference")
    return plat, {}


# ---------------------------------------------------------------------------
# PDB I/O helpers
# ---------------------------------------------------------------------------

def _save_pdb(
    topology: app.Topology,
    positions,
    path: str,
    box: List[float],
) -> None:
    with open(path, "w") as f:
        app.PDBFile.writeFile(topology, positions, f)


def _save_final_pdb(
    simulation: app.Simulation,
    topology: app.Topology,
    box: List[float],
    path: str,
) -> None:
    state = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=True
    )
    bx, by, bz = box
    state.getPositions(asNumpy=True)
    simulation.context.setPeriodicBoxVectors(
        mm.Vec3(bx, 0, 0) * unit.nanometer,
        mm.Vec3(0, by, 0) * unit.nanometer,
        mm.Vec3(0, 0, bz) * unit.nanometer,
    )
    state2 = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=True
    )
    with open(path, "w") as f:
        app.PDBFile.writeFile(topology, state2.getPositions(), f)
