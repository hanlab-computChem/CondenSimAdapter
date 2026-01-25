# CondenSimAdapter Configuration Guide

## Overview

CondenSimAdapter is a multi-stage workflow for simulating protein condensates and membrane-disrupting peptides using coarse-grained molecular dynamics simulations. The package supports two main task types:

- **IDP (Intrinsically Disordered Proteins)**: For simulating intrinsically disordered proteins and their phase separation behavior
- **MDP (Membrane-Disrupting Peptides)**: For simulating membrane-disrupting peptides and their interactions with lipid bilayers

## Configuration File Structure

The configuration is defined in a YAML file (`config.yaml`) with the following structure:

```yaml
# Input files configuration
input_files:
  structure_pdb: "protein.pdb"    # For MDP tasks
  sequence_fasta: "protein.fasta" # For IDP tasks

# Protein configuration
protein:
  name: "protein_name"
  nmol: 20

# Stage 1: Coarse-grained CALVADOS simulation
cg_calvados:
  task_type: "IDP"  # or "MDP"
  # ... CALVADOS parameters

# Future stages (to be implemented)
backmapping:
  pass

openmm_refinement:
  pass

gromacs_aa:
  pass

aa_simulation:
  pass
```

## Input Files Configuration

### `input_files` Section

This section defines the input files required for the simulation.

#### For IDP Tasks:
- **`sequence_fasta`** (string, required): Path to the FASTA file containing protein sequences. The file should contain one or more protein sequences with proper headers.

#### For MDP Tasks:
- **`structure_pdb`** (string, required): Path to the PDB file containing the protein structure. This should be a high-quality structure file suitable for coarse-grained simulation.

## Protein Configuration

### `protein` Section

This section defines the protein system parameters.

- **`name`** (string, required): The name of the protein. For IDP tasks, this should match one of the sequence names in the FASTA file. For MDP tasks, this is typically the PDB filename without the `.pdb` extension.

- **`nmol`** (integer, required): Number of protein molecules to include in the simulation. Typical values range from 10-50 depending on the system size and computational resources.

## CALVADOS Configuration

### `cg_calvados` Section

This section contains all parameters for the coarse-grained CALVADOS simulation.

#### Task Type
- **`task_type`** (string, required): Specifies the type of simulation task.
  - `"IDP"`: For intrinsically disordered protein simulations
  - `"MDP"`: For membrane-disrupting peptide simulations

#### General Parameters
- **`platform`** (string, default: "CUDA"): The computing platform to use for the simulation.
  - `"CUDA"`: Use GPU acceleration (recommended)
  - `"CPU"`: Use CPU-only computation

- **`temp`** (float, default: 310.0): Temperature of the simulation in Kelvin. Typical values range from 300-320 K.

- **`box`** (list of floats, required): Simulation box dimensions in nanometers [x, y, z]. The box should be large enough to accommodate all proteins and allow for proper dynamics.

- **`ionic`** (float, default: 0.15): Ionic strength in molar units. Typical values range from 0.1-0.2 M.

- **`pH`** (float, default: 7.0): pH of the solution. Typical values range from 6.5-7.5.

#### Simulation Parameters
- **`steps`** (integer, required): Total number of simulation steps. Each step corresponds to 10 femtoseconds.
  - For IDP tasks: Typically 200,000 - 1,000,000 steps
  - For MDP tasks: Typically 500,000 - 5,000,000 steps

- **`wfreq`** (integer, default: 1000): Frequency of trajectory writing (every N steps). Smaller values provide more frequent output but larger file sizes.

- **`logfreq`** (integer, default: 1000): Frequency of log output (every N steps).

- **`friction`** (float, default: 0.01): Friction coefficient for Langevin dynamics in ps⁻¹.

#### Topology Parameters
- **`topol`** (string, default: "random"): Initial topology for protein placement.
  - `"random"`: Random placement of proteins
  - `"slab"`: Slab topology with proteins placed in a high-density central region
  - `"grid"`: Grid topology with proteins placed in a regular grid pattern

- **`slab_width`** (float, required for slab topology): Width of the central high-density region in nanometers.
  - For IDP tasks: Default 20 nm
  - For MDP tasks: Default 40 nm (thicker for membrane interactions)

#### Checkpoint Parameters
- **`frestart`** (string, default: "restart.chk"): Name of the restart checkpoint file. The simulation automatically uses checkpoint restart mode.

#### MDP-Specific Parameters
- **`fdomains`** (string, required for MDP tasks): Path to the domains definition file (YAML format). This file defines the structural domains of the protein for harmonic restraints.

- **`restraint`** (boolean, default: true for MDP, false for IDP): Whether to apply harmonic restraints to maintain protein structure. This parameter has different default values for different task types.

#### Advanced Parameters
- **`verbose`** (boolean, default: true): Enable verbose output during simulation.

- **`slab_eq`** (boolean, default: false): Whether to perform slab equilibration.

- **`steps_eq`** (integer, default: 1000): Number of equilibration steps (fixed value, not configurable).

- **`friction`** (float, default: 0.01): Friction coefficient for Langevin dynamics in ps⁻¹.

## Example Configurations

### IDP Task Example
```yaml
input_files:
  sequence_fasta: "protein.fasta"

protein:
  name: "FUS_LC"
  nmol: 20

cg_calvados:
  task_type: "IDP"
  platform: "CUDA"
  temp: 310.0
  box: [25.0, 25.0, 50.0]
  steps: 200000
  wfreq: 5000
  topol: "slab"
  slab_width: 20
  ionic: 0.15
  pH: 7.0
  restraint: false
```

### MDP Task Example
```yaml
input_files:
  structure_pdb: "TDP43.pdb"

protein:
  name: "TDP43"
  nmol: 15

cg_calvados:
  task_type: "MDP"
  platform: "CUDA"
  temp: 310.0
  box: [10.0, 10.0, 40.0]
  steps: 500000
  wfreq: 5000
  topol: "slab"
  slab_width: 40
  fdomains: "domains.yaml"
  restraint: true
  ionic: 0.15
  pH: 7.5
```

## File Requirements

### For IDP Tasks
- FASTA file with protein sequences
- `residues_CALVADOS2.csv` (automatically provided)

### For MDP Tasks
- PDB file with protein structure
- `domains.yaml` file defining protein domains
- `residues_CALVADOS3.csv` (automatically provided)

## Output Files

The simulation generates several output files:
- **Trajectory files** (`.dcd`): Molecular dynamics trajectory
- **Checkpoint files** (`.chk`): Simulation restart points
- **Configuration files** (`.yaml`): CALVADOS configuration
- **Log files** (`.log`): Simulation progress and statistics
- **Analysis files**: Density profiles, contact maps, and other analysis results

## Parameter Usage Notes

### Actually Used Parameters
The following parameters are actually passed to CALVADOS simulations:
- `task_type`, `platform`, `temp`, `box`, `ionic`, `pH`, `topol`, `slab_width`
- `steps`, `wfreq`, `friction`, `verbose`, `slab_eq`, `frestart`
- `restraint`, `fdomains` (MDP only)

**Note**: The `friction` parameter is used as `friction` in IDP tasks and `friction_coeff` in MDP tasks.

**Important**: The package uses the CalvadosGenerator approach for parameter handling:
- **CalvadosGenerator**: Generates CALVADOS prepare.py scripts with embedded parameters

### Fixed Parameters
Some parameters have fixed values and cannot be changed:
- `restart`: Always set to 'checkpoint'
- `steps_eq`: Always set to 1000
- `runtime`: Always set to 0

### Unused Parameters
The following parameters are not currently used by the simulation:
- `eps_lj`, `logfreq`

### Auto-Generated Parameters
The following parameters are automatically generated and cannot be configured:
- `sysname`: Automatically generated as `{protein_name}_{replica_number}`

## Best Practices

1. **Box Size**: Ensure the simulation box is large enough to prevent artifacts from periodic boundary conditions.

2. **Protein Concentration**: Choose appropriate `nmol` values to achieve desired protein concentrations.

3. **Simulation Length**: Use sufficient `steps` to observe the phenomena of interest (phase separation for IDP, membrane disruption for MDP).

4. **Output Frequency**: Balance `wfreq` between file size and analysis resolution.

5. **Temperature**: Use physiological temperatures (300-320 K) for biological relevance.

6. **Ionic Strength**: Use physiological salt concentrations (0.1-0.2 M) unless studying specific salt effects.

7. **Topology**: Use 'slab' topology for most simulations, 'grid' for systematic studies, 'random' for initial exploration.
