# Tutorial 3: Production Simulation with OpenMM or AMBER

This tutorial explains how to use the `adapter to_run` command to generate production molecular dynamics (MD) scripts for OpenMM or AMBER users. This tool is optional and provides a convenient bridge between the CondenSimAdapter's `minimize` output and common MD engines.

## Overview

The `to_run` command processes the output from `adapter minimize` to create input files for either OpenMM or AMBER.

- **OpenMM**: Utilizes its native capability to process GROMACS topology and coordinate files directly.
- **AMBER**: Uses [ParmEd](https://github.com/pearlDingzhen/ParmEd.git) to convert GROMACS topology and coordinates into AMBER's `prmtop` and `inpcrd` formats.

## Prerequisites

- Output from `adapter minimize` (e.g., `minimize_final_solvated.gro`, `topol.top`, and the `.ff` directory).
- A corresponding YAML configuration file.

Example data is provided in the `example/` directory.

---

## 3.1 OpenMM Production Scripts

OpenMM has built-in support for GROMACS files. The `adapter to_run` command generates a Python script that uses OpenMM's `GromacsTopFile` and `GromacsGroFile` classes to load your system.

### Generating the Script

```bash
cd example
adapter to_run -f H1_Prota.yaml -e openmm
```

### Output Files

In the `{system_name}_production/` directory, you will find:
- `openmm_run.py`: A complete Python script for running the simulation.
- `minimize_final_solvated.gro`: Coordinate file.
- `topol.top`: Topology file.
- Force field directory (`.ff`).

### Customizing Parameters

The generated `openmm_run.py` includes default parameters (100 ns, 2 fs timestep, 300 K, NPT). You can easily modify these at the top of the script:

```python
# Simulation parameters
TOTAL_TIME_NS = 100.0
TIMESTEP_FS = 2.0
TEMPERATURE_K = 300.0
RUN_NPT = True
```

---

## 3.2 AMBER Production Files

For AMBER users, `adapter to_run` performs a format conversion. Since GROMACS and AMBER use different topology structures, we use ParmEd to handle the translation.

### Recommendation: Custom ParmEd Installation

In certain cases, standard ParmEd may not correctly handle 4-point water models (such as TIP4P) that include virtual sites (dummy atoms) during GROMACS-to-AMBER conversions. To ensure accurate topology generation for these systems, we recommend using a patched version of ParmEd:

```bash
# Uninstall existing parmed if necessary
pip uninstall parmed
# Install the patched version
pip install git+https://github.com/pearlDingzhen/ParmEd.git
```

### Generating AMBER Files

```bash
cd example
adapter to_run -f H1_Prota.yaml -e amber
```

### Output Files

The command generates:
- `system.prmtop`: AMBER topology file.
- `system.inpcrd`: AMBER coordinate file.
- `minimize_final_solvated.gro`: Original GROMACS coordinates.
- `topol.top`: Original GROMACS topology (for reference).

### Running AMBER

You will need to provide your own `.mdin` file, which typically includes the following stages: energy minimization (em), heating and equilibration (heat/eq), and production MD (nvt or npt).

Example execution:
```bash
pmemd -O -i your_input.mdin -c system.inpcrd -p system.prmtop -r restart.ncrst -x trajectory.nc
```
