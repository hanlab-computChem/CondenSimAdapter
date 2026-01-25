# CondenSimAdapter Documentation

Welcome to the CondenSimAdapter documentation! This package provides a streamlined workflow for running coarse-grained molecular dynamics simulations using CALVADOS.

## Quick Navigation

- **[README](README.md)** - Overview, installation, and quick start guide
- **[Configuration Guide](configuration.md)** - Detailed parameter descriptions for config.yaml
- **[API Documentation](api.md)** - Complete API reference and usage examples

## What is CondenSimAdapter?

CondenSimAdapter is a Python package designed to simplify the process of running coarse-grained molecular dynamics simulations for protein condensates and membrane-disrupting peptides. It provides:

- **Automated workflow** for CALVADOS simulations
- **Two task types**: IDP (Intrinsically Disordered Proteins) and MDP (Membrane-Disrupting Peptides)
- **YAML-based configuration** for easy parameter management
- **GPU acceleration** support
- **Built-in analysis tools**

## Getting Started

1. **Install the package**: See [README](README.md#installation)
2. **Create a configuration**: See [Configuration Guide](configuration.md)
3. **Run your first simulation**: See [API Documentation](api.md#usage-examples)

## Task Types

### IDP (Intrinsically Disordered Proteins)
- Simulate intrinsically disordered proteins
- Study phase separation and condensate formation
- Uses CALVADOS2 residue definitions
- Requires FASTA file with protein sequences

### MDP (Membrane-Disrupting Peptides)
- Simulate membrane-disrupting peptides
- Study interactions with lipid bilayers
- Uses CALVADOS3 residue definitions
- Requires PDB file with protein structure

## Example Configurations

### IDP Example
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
  restraint: false
```

### MDP Example
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
```

## Quick Code Example

```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator
import subprocess

# Initialize generator
generator = CalvadosGenerator("config.yaml")

# Generate and run simulation
cmd = generator.generate_and_run("output_dir", protein_name="FUS_LC", gpu_id=0, replica=1)

# Execute the simulation
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

if result.returncode == 0:
    print("Simulation completed successfully!")
else:
    print("Simulation failed:", result.stderr)
```

## Documentation Structure

### For New Users
1. Start with the **[README](README.md)** for installation and basic usage
2. Read the **[Configuration Guide](configuration.md)** to understand all available parameters
3. Use the **[API Documentation](api.md)** for detailed programming examples

### For Advanced Users
- **[API Documentation](api.md)** - Complete reference for all classes and methods
- **[Configuration Guide](configuration.md)** - Detailed parameter descriptions and best practices

### For Developers
- **[API Documentation](api.md)** - Internal API structure and extension points
- Source code in the `CondenSimAdapter/` directory

## Support and Contributing

- **Issues**: Report bugs and request features through GitHub issues
- **Contributing**: Submit pull requests for improvements
- **Documentation**: Help improve these docs by submitting corrections or additions

## Related Resources

- **CALVADOS**: The underlying simulation engine - [CALVADOS Documentation](https://github.com/julianegres/calvados)
- **Molecular Dynamics**: General concepts and theory
- **Protein Condensates**: Background on phase separation phenomena
- **Membrane Proteins**: Information on membrane-disrupting peptides

---

*This documentation is maintained as part of the CondenSimAdapter package. For the latest version, please refer to the GitHub repository.*
