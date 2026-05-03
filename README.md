# CondenSimAdapter

An automated workflow for protein condensate simulations, covering the main stages from coarse-grained (CG) to all-atom (AA).

## Installation

### Quick Start (Recommended)

**Step 1: Create environment and install heavy dependencies**

```bash
conda create -n condensim python=3.11 -y
conda activate condensim

# PyTorch + DGL (CUDA 12.x)
pip install "torch>=2.4,<2.5" --index-url https://download.pytorch.org/whl/cu124
pip install dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html

# OpenMM with CUDA 12
pip install "openmm[cuda12]>=8.2"
```

For other setups:
- **CUDA 13.x**: replace `cu124` with `cu13` and use `openmm[cuda13]`
- **AMD GPU**: `pip install "openmm[hip6]>=8.2"` or `openmm[hip7]`
- **CPU only**: use `--index-url https://download.pytorch.org/whl/cpu` and plain `openmm>=8.2`

**Step 2: Install CondenSimAdapter**

```bash
pip install CondenSimAdapter
```

**Step 3: Download neural network models (~180 MB)**

```bash
adapter models download
```

Models are hosted on GitHub Releases and cached locally.

**Step 4: Verify**

```bash
adapter --help
adapter info
adapter models status
```

### Development Installation

```bash
git clone https://github.com/hanlab-computChem/CondenSimAdapter.git
cd CondenSimAdapter

# Install heavy deps (Step 1 above), then:
pip install -e ".[ml,openmm,dev]"

# Download models
adapter models download

# Run tests
pytest tests/unit -v
```

## Usage

### Typical Workflow

```bash
# 1. Create project configuration
adapter init my_project --topol cubic -c FUS:10

# 2. Run CG simulation
adapter cg -f config.yaml

# 3. Backmap to all-atom
adapter backmap -i output_CG -f config.yaml

# 4. Energy minimization
adapter minimize -i output_backmap -f config.yaml

# 5. Generate production run scripts
adapter to_run -f config.yaml
```

### CLI Commands

```
CORE COMMANDS:
    cg               Run coarse-grained simulation
    backmap          Backmap CG structure to all-atom representation
    minimize         Energy minimization with AMBER/CHARMM force fields
    to_run           Generate production run scripts for minimize output

UTILITY COMMANDS:
    init             Initialize a new configuration template
    droplet-density  Estimate protein density in droplet geometry
    info             Display system and environment information
    forcefield       Manage custom all-atom force fields
    models           Manage neural network models for backmapping
```

## Requirements

- Python >= 3.10, < 3.12 (3.11 recommended)
- Linux x86_64
- CUDA >= 12.x capable GPU (optional, CPU mode works too)
- GROMACS >= 2023 (for topology preparation)

## Package Structure

| Component | Size | Notes |
|-----------|------|-------|
| Source code + data | ~3 MB | Installed via pip |
| Neural network models | ~235 MB | Downloaded via `adapter models download` |
| PyTorch + DGL + OpenMM | ~2 GB | Installed via pip (Step 1) |

## Links

- PyPI: https://pypi.org/project/CondenSimAdapter/
- GitHub: https://github.com/hanlab-computChem/CondenSimAdapter
- Issues: https://github.com/hanlab-computChem/CondenSimAdapter/issues

## License

GPL-3.0
