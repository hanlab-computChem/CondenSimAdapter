# Conda-forge Build Instructions

## Quick Start for Users

### One-Line Installation (After Published)

```bash
# Create and activate environment
conda create -n csa -c conda-forge conden-sim-adapter
conda activate csa

# Done! Test it
adapter --help
adapter --version
```

That's it! Everything is included:
- ✅ All Python dependencies
- ✅ OpenMM with CUDA 12.1 support
- ✅ PyTorch with GPU support
- ✅ GROMACS for AA simulations
- ✅ Neural network models (189MB)

## Build Locally (For Developers)

### Prerequisites

```bash
# Install conda-build
conda install conda-build

# Optional: Install mamba for faster builds
conda install mamba
```

### Build Steps

```bash
# Navigate to recipe directory
cd conda-forge-recipe

# Build the package
conda build .

# Or use mamba (faster)
mamba build .

# The built package will be in:
# ~/miniconda3/conda-bld/noarch/conden-sim-adapter-*.tar.bz2
```

### Test Local Build

```bash
# Create test environment with local package
conda create -n csa-test -c local conden-sim-adapter

# Or install directly
conda install ~/miniconda3/conda-bld/noarch/conden-sim-adapter-*.tar.bz2
```

### Convert to Other Platforms (Optional)

```bash
# Install conda-convert
conda install conda-convert

# Convert to other platforms
conda convert --platform all ~/miniconda3/conda-bld/noarch/conden-sim-adapter-*.tar.bz2 -o outputdir/
```

## Submit to Conda-forge

### Step 1: Fork and Clone

```bash
# Fork https://github.com/conda-forge/staged-recipes

git clone https://github.com/YOUR_USERNAME/staged-recipes.git
cd staged-recipes
```

### Step 2: Copy Recipe

```bash
mkdir -p recipes/conden-sim-adapter
cp /path/to/CondenSimAdapter/conda-forge-recipe/meta.yaml recipes/conden-sim-adapter/
cp /path/to/CondenSimAdapter/conda-forge-recipe/build.sh recipes/conden-sim-adapter/

# Add LICENSE if not in source
cp /path/to/CondenSimAdapter/LICENSE recipes/conden-sim-adapter/
```

### Step 3: Update meta.yaml for Release

Edit `recipes/conden-sim-adapter/meta.yaml`:

```yaml
source:
  url: https://github.com/hanlab-computChem/CondenSimAdapter/archive/refs/tags/v0.2.0.tar.gz
  sha256: REPLACE_WITH_ACTUAL_SHA256  # Get from: sha256sum v0.2.0.tar.gz
  # Remove: path: ../
```

### Step 4: Submit PR

```bash
git checkout -b add-conden-sim-adapter
git add recipes/conden-sim-adapter/
git commit -m "Add conden-sim-adapter recipe"
git push origin add-conden-sim-adapter
```

Create PR on GitHub. After merge, the package will be available on conda-forge.

## Package Size Breakdown

| Component | Size | Included |
|-----------|------|----------|
| Source code | ~1 MB | ✅ |
| Model files | ~189 MB | ✅ |
| Total conda package | ~200 MB | ✅ |
| Dependencies (downloaded) | ~2-3 GB | ✅ (CUDA, PyTorch, etc.) |

Total download: ~2.5 GB (one-time, includes everything)

## Troubleshooting

### Build Fails: "No module named 'conda_build'"

```bash
conda install conda-build
```

### Test Fails: Models Not Found

Ensure models are in the source:
```bash
ls CondenSimAdapter/backmap/cg2all/model/*.ckpt
# Should show: CalphaBasedModel.ckpt, Martini.ckpt, etc.
```

### CUDA Version Mismatch

The recipe pins `pytorch-cuda =12.1`. If you need different CUDA:
- For CUDA 11.8: Change to `pytorch-cuda =11.8`
- For CPU only: Remove `pytorch-cuda` line

### Size Limit Issues

Conda-forge has no strict size limit (unlike PyPI's 100MB). 
200MB is acceptable for a scientific package.

## Maintenance

### Update Version

1. Update version in `meta.yaml`
2. Update SHA256 hash
3. Reset build number to 0
4. Submit PR to conda-forge feedstock

### Add New Dependencies

Add to `requirements/run` in `meta.yaml`.

### Update Models

1. Update model files in source
2. Increase build number in `meta.yaml`
3. Rebuild and submit
