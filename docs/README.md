# CondenSimAdapter

A multi-stage workflow for simulating protein condensates and membrane-disrupting peptides using coarse-grained molecular dynamics simulations.

## Overview

CondenSimAdapter is a Python package that provides a streamlined workflow for running coarse-grained molecular dynamics simulations using CALVADOS. The package supports two main types of simulations:

- **IDP (Intrinsically Disordered Proteins)**: Simulations of intrinsically disordered proteins to study phase separation and condensate formation
- **MDP (Membrane-Disrupting Peptides)**: Simulations of membrane-disrupting peptides to study their interactions with lipid bilayers

## Features

- **Automated workflow**: Generate CALVADOS scripts from simple YAML configuration files
- **Two task types**: Support for both IDP and MDP simulations
- **Flexible configuration**: Comprehensive parameter control through YAML files
- **GPU acceleration**: Native support for CUDA-accelerated simulations
- **Analysis integration**: Built-in analysis tools for density profiles, contact maps, and more

## Installation

### Prerequisites

1. **CALVADOS**: Install CALVADOS in your conda environment
   ```bash
   conda install -c conda-forge calvados
   ```

2. **CUDA**: For GPU acceleration, ensure CUDA is properly installed

### Install CondenSimAdapter

```bash
# Clone the repository
git clone <repository-url>
cd CondenSimAdapter

# Install in development mode
pip install -e .
```

## Quick Start

### 1. Create a Configuration File

Create a `config.yaml` file for your simulation:

```yaml
# For IDP simulation
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
  restraint: false
```

### 2. Run the Simulation

```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator

# Initialize generator
generator = CalvadosGenerator("config.yaml")

# Generate and run simulation
cmd = generator.generate_and_run("output_dir", protein_name="FUS_LC", gpu_id=0, replica=1)
```

### 3. Example Scripts

Use the provided example scripts for quick testing:

```bash
# Run IDP example
cd example/IDP
python run_example.py

# Run MDP example
cd example/MDP
python run_example.py
```

## Configuration

The package uses YAML configuration files to define simulation parameters. See the [Configuration Guide](configuration.md) for detailed parameter descriptions.

### Key Configuration Sections

- **`input_files`**: Define input files (FASTA for IDP, PDB for MDP)
- **`protein`**: Protein system parameters (name, number of molecules)
- **`cg_calvados`**: CALVADOS simulation parameters (temperature, box size, steps, etc.)

## Task Types

### IDP (Intrinsically Disordered Proteins)

For simulating intrinsically disordered proteins and studying their phase separation behavior:

- Uses CALVADOS2 residue definitions
- Requires FASTA file with protein sequences
- Optimized for studying protein condensates
- Typical simulation length: 200,000 - 1,000,000 steps

### MDP (Membrane-Disrupting Peptides)

For simulating membrane-disrupting peptides and their interactions with lipid bilayers:

- Uses CALVADOS3 residue definitions
- Requires PDB file with protein structure
- Requires domains.yaml file for structural restraints
- Optimized for studying membrane interactions
- Typical simulation length: 500,000 - 5,000,000 steps

## Output Files

The simulation generates several output files:

- **Trajectory files** (`.dcd`): Molecular dynamics trajectory
- **Checkpoint files** (`.chk`): Simulation restart points
- **Configuration files** (`.yaml`): CALVADOS configuration
- **Log files** (`.log`): Simulation progress and statistics
- **Analysis files**: Density profiles, contact maps, and other analysis results

## Examples

### IDP Example

```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator

# Initialize for IDP simulation
generator = CalvadosGenerator("config_idp.yaml")

# Generate and run
cmd = generator.generate_and_run("output_idp", protein_name="FUS_LC", gpu_id=0, replica=1)
```

### MDP Example

```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator

# Initialize for MDP simulation
generator = CalvadosGenerator("config_mdp.yaml")

# Generate and run
cmd = generator.generate_and_run("output_mdp", protein_name="TDP43", gpu_id=0, replica=1)
```

## Best Practices

1. **Environment Setup**: Always use a dedicated conda environment with CALVADOS installed
2. **Box Size**: Ensure simulation box is large enough to prevent periodic boundary artifacts
3. **Protein Concentration**: Choose appropriate molecule counts for your system
4. **Simulation Length**: Use sufficient steps to observe phenomena of interest
5. **Output Frequency**: Balance between file size and analysis resolution
6. **Temperature**: Use physiological temperatures (300-320 K) for biological relevance

## Troubleshooting

### Common Issues

1. **CALVADOS not found**: Ensure CALVADOS is installed in your conda environment
2. **CUDA errors**: Check CUDA installation and GPU availability
3. **File not found**: Verify input file paths in configuration
4. **Memory issues**: Reduce box size or number of molecules

### Getting Help

- Check the [Configuration Guide](configuration.md) for parameter details
- Review example configurations in the `example/` directory
- Ensure all required files are present and properly formatted

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## License

[Add your license information here]

## Citation

If you use CondenSimAdapter in your research, please cite:

[Add citation information here]
