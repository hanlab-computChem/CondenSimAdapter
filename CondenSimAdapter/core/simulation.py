"""
Unified CG simulation driver.

Replaces the five separate run_*() methods that previously lived in
the 3358-line cg.py.  All force fields go through the same pipeline:

  build coords -> build topology -> assemble forces -> run MD
"""

from __future__ import annotations

import logging
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
from .entanglement import EntanglementAnalyzer
from .forcefield import create_forcefield
from .forcefield.base import CGForceField
from .forcefield.cocomo import CocomoFF
from .molecule import build_all_chains
from .topology import build_topology
from .z1plus import Z1PlusWrapper

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
        # Use resolved_force_field so that bare 'calvados' is auto-versioned
        # (CALVADOS2 for all-IDP systems, CALVADOS3 when MDP components are present)
        self.ff: CGForceField = create_forcefield(config.resolved_force_field)
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
                f"Output directory '{out}' already exists. Pass overwrite=True to proceed."
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

        openmm_platform, platform_props, actual_platform = _resolve_platform(
            sim_params.platform, gpu_id
        )
        print(f"  Platform:    {actual_platform}", flush=True)
        simulation = app.Simulation(topology, system, integrator, openmm_platform, platform_props)
        simulation.context.setPositions(pos_qty)

        # 6. Energy minimise before production
        # Use 10 000 iterations to handle large initial bond deviations that arise
        # when MDP chains have compact IDR stubs displaced from the folded domain.
        log.info("Energy minimisation ...")
        simulation.minimizeEnergy(maxIterations=10000)

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
        remainder = sim_params.steps % batch

        with tqdm(total=sim_params.steps, unit="steps", smoothing=0.05) as pbar:
            for _ in range(n_batches):
                simulation.step(batch)
                simulation.saveCheckpoint(chk_file)
                pbar.update(batch)

            if remainder > 0:
                simulation.step(remainder)
                simulation.saveCheckpoint(chk_file)
                pbar.update(remainder)

        # 8. Save final structure
        final_pdb = str(out / "final.pdb")
        _save_final_pdb(simulation, topology, cfg.box, final_pdb)
        log.info(f"Simulation complete. Final PDB: {final_pdb}")

        # 9. Entanglement check (warn-only; never blocks the run)
        ent_stats: dict = {}
        if cfg.check_entanglement:
            try:
                ent_stats = _check_entanglement(final_pdb, cfg, chain_meta)
            except Exception as exc:  # noqa: BLE001
                log.warning("Entanglement check failed (non-fatal): %s", exc)

        return SimulationResult(
            success=True,
            output_dir=str(out),
            final_pdb=final_pdb,
            trajectory=xtc_file,
            log_file=log_file,
            entanglement_mean_z=ent_stats.get("mean_z"),
            entanglement_max_z=ent_stats.get("max_z"),
            entanglement_method=ent_stats.get("method"),
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
        # Collect ENM pairs for exclusion from nonbonded forces, matching the
        # original OpenMpipi behavior (topology.addBond for ENM → createExclusionsFromBonds)
        # and CALVADOS behavior (add_exclusions called for every restrained pair).
        has_domains = any(meta["folded_domains"] for meta in chain_meta)
        enm_pairs: list = []
        if has_domains:
            enm = self.ff.build_enm_bonds(positions, chain_meta, restraint_type="harmonic")
            if enm is not None:
                system.addForce(enm)
                enm_pairs = _extract_bond_pairs(enm)

        # 6. Non-bonded forces (force field specific)
        if isinstance(self.ff, CocomoFF):
            nb_forces = self.ff.create_nonbonded_forces(
                topology,
                chain_meta,
                config.temperature,
                config.ionic_strength,
                positions=positions,
            )
        else:
            nb_forces = self.ff.create_nonbonded_forces(
                topology,
                chain_meta,
                config.temperature,
                config.ionic_strength,
            )
        for f in nb_forces:
            # Conditionally add ENM pair exclusions to every CustomNonbondedForce.
            # - mpipi / calvados / hps: exclude ENM pairs (matching OpenMpipi which
            #   adds all bonds to topology before createExclusionsFromBonds, and
            #   CALVADOS which calls add_exclusions for every restrained pair).
            # - cocomo: do NOT exclude ENM pairs (original COCOMO intentionally lets
            #   the 10-5 LJ act on native contacts, adding a small attractive
            #   contribution that deepens the native-contact potential well).
            if (
                enm_pairs
                and self.ff._exclude_enm_from_nonbonded
                and isinstance(f, mm.CustomNonbondedForce)
            ):
                for i, j in enm_pairs:
                    try:
                        f.addExclusion(i, j)
                    except Exception:
                        pass  # pair already in exclusion list
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
    """Select the best available OpenMM platform with a real GPU probe.

    Unlike getPlatformByName() alone, this actually tries to create a tiny
    mm.Context to confirm the GPU is usable before returning it.  Falls back
    to CPU and prints a warning if the GPU fails.

    Returns (platform, properties, actual_platform_name).
    """
    preferred = platform_name.upper()
    plat_name_map = {"CUDA": "CUDA", "OPENCL": "OpenCL"}

    if preferred in plat_name_map:
        try:
            plat = mm.Platform.getPlatformByName(plat_name_map[preferred])
            props = {"DeviceIndex": str(gpu_id), "Precision": "mixed"}
            # Real probe: create a tiny context to verify the GPU is accessible
            _sys = mm.System()
            _sys.addParticle(1.0)
            _ctx = mm.Context(_sys, mm.VerletIntegrator(0.001), plat, props)
            del _ctx, _sys
            return plat, props, f"{plat_name_map[preferred]}:{gpu_id}"
        except Exception as e:
            print(
                f"  [warn] {preferred} GPU {gpu_id} unavailable ({e}); falling back to CPU.",
                flush=True,
            )

    try:
        plat = mm.Platform.getPlatformByName("CPU")
    except Exception:
        plat = mm.Platform.getPlatformByName("Reference")
    return plat, {}, plat.getName()


# ---------------------------------------------------------------------------
# PDB I/O helpers
# ---------------------------------------------------------------------------


def _extract_bond_pairs(force: mm.Force) -> list:
    """
    Extract (i, j) index pairs from a HarmonicBondForce or CustomBondForce.

    Used to collect ENM bond pairs so they can be added as exclusions to
    CustomNonbondedForce objects (WF, Yukawa) — matching both the original
    OpenMpipi approach (topology.addBond + createExclusionsFromBonds) and
    CALVADOS (explicit add_exclusions for every restrained pair).
    """
    pairs = []
    if isinstance(force, mm.HarmonicBondForce):
        for k in range(force.getNumBonds()):
            p1, p2, *_ = force.getBondParameters(k)
            pairs.append((int(p1), int(p2)))
    elif isinstance(force, mm.CustomBondForce):
        for k in range(force.getNumBonds()):
            p1, p2, *_ = force.getBondParameters(k)
            pairs.append((int(p1), int(p2)))
    return pairs


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
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    bx, by, bz = box
    state.getPositions(asNumpy=True)
    simulation.context.setPeriodicBoxVectors(
        mm.Vec3(bx, 0, 0) * unit.nanometer,
        mm.Vec3(0, by, 0) * unit.nanometer,
        mm.Vec3(0, 0, bz) * unit.nanometer,
    )
    state2 = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    with open(path, "w") as f:
        app.PDBFile.writeFile(topology, state2.getPositions(), f)


# ---------------------------------------------------------------------------
# Entanglement check helpers
# ---------------------------------------------------------------------------

# Unified per-system mean-Z thresholds (apply to both built-in and Z1+).
# Based on benchmark of 32 equilibrated condensate systems:
#   built-in range: 1.4–6.9 ; Z1+ range: 1.1–6.3
# Thresholds are calibrated on the built-in algorithm; when Z1+ is used the
# values are ~1.4x lower (see plan for derivation).
_ENT_THRESH_ELEVATED = 4.0  # exceeds most IDP baselines
_ENT_THRESH_HIGH = 7.0  # exceeds all benchmark systems → likely artifact


def _load_pdb_for_entanglement(
    pdb_path: str,
    chain_meta: List[dict],
) -> tuple:
    """
    Load a PDB written by the CG pipeline and return (positions, boundaries, box_angstrom).

    Uses chain_meta to map from the segment order to per-chain boundaries so
    we do not depend on MDAnalysis segment naming.  Positions are returned
    in Angstroms (as in the PDB) because EntanglementAnalyzer and Z1PlusWrapper
    both accept Angstroms.
    """
    import warnings

    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise ImportError("MDAnalysis is required for the entanglement check.") from exc

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(pdb_path)

    box_ang = u.dimensions[:3].astype(np.float64)  # Angstroms

    all_pos: list = []
    boundaries: list = []
    offset = 0
    # Rely on segment order (same order as chain_meta) rather than segids,
    # because the CG PDB reuses single-letter segids for many chains.
    for seg in u.segments:
        ca = seg.atoms.select_atoms("name CA")
        if len(ca) == 0:
            continue
        all_pos.append(ca.positions.astype(np.float64))
        boundaries.append((offset, offset + len(ca)))
        offset += len(ca)

    if not all_pos:
        raise ValueError("No CA atoms found in %s" % pdb_path)

    positions = np.vstack(all_pos)
    return positions, boundaries, box_ang


def _check_entanglement(
    final_pdb: str,
    cfg: CGConfig,
    chain_meta: List[dict],
) -> dict:
    """
    Run the entanglement check on *final_pdb* and return a stats dictionary.

    Priority:
      1. Z1+ (if cfg.z1plus_executable is set and points to an executable)
      2. Built-in Z-code PPA (EntanglementAnalyzer)

    Prints a formatted summary to stdout and logs warnings at the appropriate
    level.  Never raises — callers should catch and log any exception.

    Returns
    -------
    dict with keys: mean_z, max_z, method, z_values (np.ndarray), n_chains
    """
    positions, boundaries, box_ang = _load_pdb_for_entanglement(final_pdb, chain_meta)
    n_chains = len(boundaries)

    # Choose analysis method
    method = "builtin"
    use_z1plus = False
    if cfg.z1plus_executable:
        wrapper = Z1PlusWrapper(executable=cfg.z1plus_executable)
        if wrapper.is_available():
            use_z1plus = True
            method = "z1plus"
        else:
            log.warning(
                "z1plus_executable '%s' is not executable; falling back to built-in Z-code PPA.",
                cfg.z1plus_executable,
            )

    if use_z1plus:
        result = wrapper.run(positions, boundaries, box_ang.tolist())
        z_values = result["z_values"]
        mean_z = float(result["mean_z"])
        max_z = float(result["max_z"])
        method_label = "Z1+"
    else:
        analyzer = EntanglementAnalyzer(positions, boundaries, box_ang, use_pbc=True)
        report = analyzer.run(max_iter=100)
        z_values = report.z_values
        mean_z = float(report.mean_z)
        max_z = float(report.max_z)
        method_label = "built-in Z-code PPA"

    n_entangled = int(np.sum(z_values > 0))
    frac_entangled = 100.0 * n_entangled / max(n_chains, 1)

    stats = {
        "mean_z": mean_z,
        "max_z": max_z,
        "method": method,
        "z_values": z_values,
        "n_chains": n_chains,
    }

    # --- Print summary (always) ---
    print(
        f"\n[Entanglement check]  Method: {method_label}"
        f"{'  (Z1+ not found; using built-in)' if method == 'builtin' and cfg.z1plus_executable else ''}"
    )
    print(f"  Chains analysed : {n_chains}")
    print(f"  Mean Z          : {mean_z:.2f}")
    print(f"  Max Z           : {int(max_z)}")
    print(f"  Fraction Z > 0  : {frac_entangled:.1f}%")

    _log_entanglement_verdict(stats)
    return stats


def _log_entanglement_verdict(stats: dict) -> None:
    """
    Print a tiered verdict line based on unified mean-Z thresholds.

    Thresholds (calibrated on built-in; Z1+ values ~1.4x lower):
      OK       : mean Z ≤ 0.5
      NORMAL   : 0.5 < mean Z ≤ 4.0
      ELEVATED : 4.0 < mean Z ≤ 7.0
      HIGH     : mean Z > 7.0
    """
    mean_z = stats["mean_z"]

    if mean_z <= 0.5:
        verdict = "OK"
        detail = "No significant entanglement detected."
        log.info("[Entanglement] %s  mean Z = %.2f", verdict, mean_z)
        print(f"  Verdict         : {verdict} — {detail}")

    elif mean_z <= _ENT_THRESH_ELEVATED:
        verdict = "NORMAL"
        detail = "Within the expected range for condensate simulations."
        log.info("[Entanglement] %s  mean Z = %.2f", verdict, mean_z)
        print(f"  Verdict         : {verdict} — {detail}")

    elif mean_z <= _ENT_THRESH_HIGH:
        verdict = "ELEVATED"
        detail = (
            f"Mean Z ({mean_z:.2f}) is moderately high.  For IDP-only systems this "
            "may indicate topology artifacts; for MDP systems with long IDRs it can "
            "be physically expected.  Consider inspecting the structure if unexpected."
        )
        msg = f"[Entanglement] {verdict} — {detail}"
        log.warning(msg)
        print(f"  Verdict         : *** {verdict} ***\n                    {detail}")

    else:
        verdict = "HIGH"
        detail = (
            f"Mean Z ({mean_z:.2f}) exceeds the range observed in all 32 benchmark "
            "condensate systems (max: 6.9 built-in / 6.3 Z1+).  Possible causes:\n"
            "                      - Chains threading through folded domains during IDR placement\n"
            "                      - Loop/ring peptide topological linking\n"
            "                    Recommendation: inspect the final structure visually "
            "and consider re-running with a lower initial density or a different "
            "placement strategy.\n"
            "                    To disable this check: set  check_entanglement: false  "
            "in your config file."
        )
        msg = f"[Entanglement] {verdict} — mean Z = {mean_z:.2f}"
        log.warning(msg)
        print(
            f"  Verdict         : *** WARNING: {verdict} ENTANGLEMENT ***\n"
            f"                    {detail}"
        )
