# CondenSimAdapter Installation Guide

## 🚀 Quick Start

### Option 1: pip (Simplest - Recommended)

```bash
# Install PyTorch + dependencies first (one-time setup)
pip install torch==2.1.2 --index-url https://download.pytorch.org/whl/cu121  # or /cpu for CPU
pip install dgl -f https://data.dgl.ai/wheels/torch-2.1/cu121/repo.html

# Install CondenSimAdapter (includes ~200MB models)
pip install CondenSimAdapter

# Done!
adapter --help
```

Or use the install script:
```bash
bash install-pip.sh gpu   # For GPU
bash install-pip.sh cpu   # For CPU only
```

### Option 2: conda (Coming Soon)

```bash
conda create -n csa -c conda-forge conden-sim-adapter
conda activate csa
```

## 📦 What's Included

The pip package (~200 MB) contains everything:

| Component | Size | Included |
|-----------|------|----------|
| Source code + data | ~10 MB | ✅ |
| Neural network models | ~189 MB | ✅ |
| **Total** | **~200 MB** | **✅ One download** |

Models included:
- `CalphaBasedModel.ckpt` - Standard CA-based CG backmapping
- `CalphaBasedModel-FIX.ckpt` - With CA position fixing
- `ResidueBasedModel.ckpt` - Residue-level CG
- `Martini.ckpt` - Martini 2 backmapping
- `Martini3.ckpt` - Martini 3 backmapping

## 🔧 Requirements

- Python 3.10 or 3.11
- Linux x86_64 (primary support)
- 500 MB disk space
- CUDA 12.1 capable GPU (optional, CPU mode works too)

## 📋 Detailed Steps

### Step 1: Prepare Environment

```bash
# Create virtual environment (recommended)
python3.11 -m venv csa-env
source csa-env/bin/activate

# Or use conda
conda create -n csa python=3.11
conda activate csa
```

### Step 2: Install PyTorch Stack

**For GPU (CUDA 12.1):**
```bash
pip install torch==2.1.2 --index-url https://download.pytorch.org/whl/cu121
pip install dgl -f https://data.dgl.ai/wheels/torch-2.1/cu121/repo.html
```

**For CPU only:**
```bash
pip install torch==2.1.2 --index-url https://download.pytorch.org/whl/cpu
pip install dgl -f https://data.dgl.ai/wheels/torch-2.1/repo.html
```

### Step 3: Install CondenSimAdapter

```bash
# Install from PyPI
pip install CondenSimAdapter

# Or install with all simulation dependencies
pip install "CondenSimAdapter[all]"
```

This downloads ~200MB including all neural network models.

### Step 4: Verify Installation

```bash
# Check CLI
adapter --version
adapter --help

# Check models
adapter models status

# Run tests (optional)
pip install pytest
pytest tests/unit -v --tb=short
```

## 🎯 First Workflow

```bash
# 1. Create project
adapter init my_project --topol cubic -c FUS:10
cd my_project

# 2. Edit config.yaml if needed

# 3. Run CG simulation
adapter cg -f config.yaml

# 4. Backmap to all-atom
adapter backmap -i output_CG -f config.yaml

# 5. Minimize
adapter minimize -i output_backmap -f config.yaml
```

## 🔍 Troubleshooting

### ImportError: No module named 'torch'

Install PyTorch first (see Step 2 above).

### CUDA errors

Check CUDA availability:
```bash
python -c "import torch; print('CUDA:', torch.cuda.is_available())"
```

If False, reinstall with CPU version or check GPU drivers.

### Model files not found

If models are missing (shouldn't happen with pip install):
```bash
adapter models download
```

### Permission denied

Use `--user` flag:
```bash
pip install --user CondenSimAdapter
```

Or use virtual environment (recommended).

## 📦 Package Size Note

The package is ~200 MB because it includes neural network models. This is larger than typical Python packages but:

- **One-time download** - Models included, no separate downloads
- **Academic standard** - Similar to PyTorch, TensorFlow model packages
- **Fast startup** - No waiting for model downloads on first use

## 🆚 pip vs conda

| Feature | pip | conda |
|---------|-----|-------|
| **Install speed** | ⭐⭐⭐ Fast | ⭐⭐ Slower |
| **CUDA handling** | Manual | Automatic |
| **Package size** | ~200 MB | ~200 MB |
| **Dependencies** | Need PyTorch first | All included |
| **Best for** | Users with PyTorch already | Fresh installs |

## 📝 Development Install

For contributing to the project:

```bash
git clone https://github.com/hanlab-computChem/CondenSimAdapter.git
cd CondenSimAdapter

pip install torch  # or your preferred PyTorch version
pip install -e ".[sim,ml,dev]"

# Run tests
pytest tests/ -v
```

## 🔗 Links

- PyPI: https://pypi.org/project/CondenSimAdapter/
- GitHub: https://github.com/hanlab-computChem/CondenSimAdapter
- Issues: https://github.com/hanlab-computChem/CondenSimAdapter/issues
- Documentation: See `Tutorials/` directory

## 💡 Pro Tips

1. **Use virtual environments** - Keeps dependencies isolated
2. **Pin PyTorch version** - Use `torch==2.1.2` for compatibility
3. **GPU vs CPU** - CPU version is fine for small systems
4. **Check disk space** - Ensure 1 GB free for installation
