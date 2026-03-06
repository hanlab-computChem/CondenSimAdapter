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

### Installation Options

#### Core Package Only
```bash
pip install CondenSimAdapter
```

#### With Simulation Dependencies (recommended)
```bash
pip install "CondenSimAdapter[sim]"
```
This includes OpenMM, MDAnalysis, mdtraj, and other simulation tools.

#### With Machine Learning Dependencies (for neural network backmapping)
```bash
pip install "CondenSimAdapter[ml]"
```
This includes PyTorch, DGL, and e3nn for CG-to-AA backmapping.

#### Full Installation (everything)
```bash
pip install "CondenSimAdapter[all]"
```

## Development Installation

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

- **PyPI**: https://pypi.org/project/CondenSimAdapter/
- **Documentation**: https://github.com/hanlab-computChem/CondenSimAdapter#readme
- **Issues**: https://github.com/hanlab-computChem/CondenSimAdapter/issues

## License

GPL-3.0
