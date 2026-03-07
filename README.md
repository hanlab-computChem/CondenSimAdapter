# CondenSimAdapter

CondenSimAdapter is an automated workflow for protein condensate simulations, covering the main stages from coarse-grained (CG) to all-atom (AA).

## Installation

### Quick Start (Recommended)

**Step 1: Create conda environment and install heavy dependencies**

```bash
# Create environment with Python 3.11
conda create -n condensim python=3.11 -y
conda activate condensim

# Install OpenMM (from conda-forge)
conda install -c conda-forge openmm=8.2.0

# Install PyTorch (adjust CUDA version as needed, here using CUDA 12.1)
conda install pytorch=2.4.1 pytorch-cuda=12.1 -c pytorch -c nvidia

# Install DGL for CUDA 12.1
conda install -c dglteam/label/cu121 dgl=1.1.3
```

For different CUDA versions, adjust the PyTorch and DGL installation:
- **CUDA 11.8**: `pytorch-cuda=11.8` and `conda install -c dglteam/label/cu118 dgl=1.1.3`
- **CPU only**: `conda install pytorch=2.4.1 cpuonly -c pytorch` and `conda install -c dglteam dgl=1.1.3`

**Step 2: Install CondenSimAdapter from PyPI**

```bash
pip install CondenSimAdapter
```

This installs the core package with ML backmapping support (requires Step 1 conda dependencies to be pre-installed).

**Step 3: Verify installation**

```bash
adapter --version
```

### Alternative: Pure pip Installation (Not Recommended)

If you cannot use conda, you can install ML dependencies via pip:

```bash
pip install CondenSimAdapter[ml]
```

⚠️ **Warning**: Installing PyTorch and DGL via pip may fail or cause CUDA compatibility issues. Use conda installation (above) for best results.

### Minimal Installation (No ML Backmapping)

If you don't need backmapping functionality:

```bash
pip install CondenSimAdapter[minimal]
```

### Development Installation

```bash
# 1. Follow Step 1 above to install conda dependencies

# 2. Clone and install in editable mode
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
- CUDA >= 12.1 for cg2all backmapping (adjust conda packages for your CUDA version)
- GROMACS >= 2023 (install separately)

## Links

- PyPI: https://pypi.org/project/CondenSimAdapter/
- Documentation: https://github.com/hanlab-computChem/CondenSimAdapter#readme
- Issues: https://github.com/hanlab-computChem/CondenSimAdapter/issues

## License

GPL-3.0
