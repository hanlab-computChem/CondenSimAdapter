# CondenSimAdapter

CondenSimAdapter is an automated workflow for protein condensate simulations, covering the main stages from coarse-grained (CG) to all-atom (AA).

## Installation

### Quick Start

```bash
# 1. Create a conda environment
conda create -n conden python=3.11 -y
conda activate conden

# 2. Install from PyPI
pip install CondenSimAdapter

# 3. Verify installation
adapter --version
```

### Standard Installation

```bash
pip install CondenSimAdapter
```

This installs everything including:

- Core simulation tools (OpenMM, MDAnalysis, mdtraj)
- Neural network backmapping (PyTorch, DGL, e3nn)

Requirements: CUDA-capable GPU with CUDA >= 12.1

### Development Installation

If you want to modify the source code:

```bash
git clone https://github.com/hanlab-computChem/CondenSimAdapter.git
cd CondenSimAdapter
pip install -e ".[dev]"
```

## Testing

```bash
# Run tests
pytest tests/

# Run with coverage
pytest tests/ --cov=CondenSimAdapter
```

## Usage

### Command Line Interface

```bash
# Show help
adapter --help

# Run backmapping
adapter backmap -c cg_structure.pdb -o aa_structure.pdb

# Check model status
adapter models status
```

### Python API

```python
from CondenSimAdapter import backmap

# Backmap CG structure to AA
backmap.convert("cg_input.pdb", "aa_output.pdb")
```

## Requirements

- Python >= 3.10, < 3.12 (3.11 recommended)
- CUDA >= 12.1 (for ML features)
- GROMACS >= 2023 (install separately)

## Links

- PyPI: https://pypi.org/project/CondenSimAdapter/
- Documentation: https://github.com/hanlab-computChem/CondenSimAdapter#readme
- Issues: https://github.com/hanlab-computChem/CondenSimAdapter/issues

## License

GPL-3.0
