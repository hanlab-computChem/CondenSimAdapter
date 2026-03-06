# PyPI Release Guide

## Quick Start

```bash
# 1. Check package size
bash check_package_size.sh

# 2. If size > 100MB, request limit increase (see below)

# 3. Build and upload
python -m build
twine check dist/*
twine upload dist/*
```

## Package Size

| Component | Size |
|-----------|------|
| Source code + data | ~10 MB |
| Model files (.ckpt) | ~189 MB |
| **Total** | **~200 MB** |

PyPI default limit: **100 MB per file**

## Option 1: Request PyPI Size Limit Increase (Recommended)

This is the cleanest solution. PyPI supports size increases for legitimate use cases.

### Step 1: Create Issue

Go to: https://github.com/pypi/support/issues/new?select=limit-request.md

### Step 2: Fill in Template

```markdown
**Project:** CondenSimAdapter
**Size of release:** ~200 MB
**Which indexes:** PyPI (and TestPyPI if possible)

**Description:**

CondenSimAdapter is a scientific software package for protein condensate 
simulations, covering the workflow from coarse-grained (CG) to all-atom (AA).

**Why the large size:**

1. **Neural network models (~189 MB):** The package includes pre-trained 
   neural network checkpoint files for CG-to-AA backmapping. These are 
   inherently large due to the model architecture.

2. **Academic software conventions:** Users in the computational 
   chemistry/biology community expect a complete, working package that 
   can be installed with a single command.

3. **No viable alternatives:** 
   - External hosting adds complexity for academic users
   - Download-on-first-use breaks offline installations
   - Splitting the package would fragment the user experience

**Project details:**
- Repository: https://github.com/hanlab-computChem/CondenSimAdapter
- License: GPL-3.0
- Domain: Computational chemistry / Molecular simulation
- Target users: Academic researchers

**Request:** Please increase the size limit to 300 MB to accommodate 
future model updates.
```

### Step 3: Wait for Approval

- Usually takes 1-3 days
- PyPI maintainers are generally supportive of scientific packages
- Once approved, you can upload large files immediately

## Option 2: Upload to TestPyPI First

Even with size limit, you should test on TestPyPI:

```bash
# Build
python -m build

# Upload to TestPyPI
twine upload --repository testpypi dist/*

# Install from TestPyPI and verify
pip install --index-url https://test.pypi.org/simple/ CondenSimAdapter
adapter --version
```

## Option 3: Manual Upload (If GitHub Actions Fails)

```bash
# Set up API token first:
# https://pypi.org/manage/account/token/

# Configure twine
cat > ~/.pypirc << EOF
[pypi]
username = __token__
password = pypi-your-token-here
EOF

# Build and upload
python -m build
twine upload dist/*
```

## Complete Release Checklist

### Before Release

- [ ] Version updated in `pyproject.toml`
- [ ] CHANGELOG.md updated
- [ ] All tests passing: `pytest tests/ -v`
- [ ] Package size checked: `bash check_package_size.sh`
- [ ] If size > 100MB: PyPI limit increase requested

### Build

- [ ] Clean build: `rm -rf build/ dist/ *.egg-info/`
- [ ] Build: `python -m build`
- [ ] Check: `twine check dist/*`
- [ ] Test install: `pip install dist/*.whl`

### Upload

- [ ] TestPyPI upload: `twine upload --repository testpypi dist/*`
- [ ] Test from TestPyPI: `pip install --index-url https://test.pypi.org/simple/ CondenSimAdapter`
- [ ] PyPI upload: `twine upload dist/*`

### Verification

- [ ] Package visible on https://pypi.org/project/CondenSimAdapter/
- [ ] Installation works: `pip install CondenSimAdapter`
- [ ] CLI works: `adapter --version`
- [ ] Models included: `adapter models status`

### Announce

- [ ] GitHub Release created with notes
- [ ] Tag pushed: `git push origin v0.2.0`
- [ ] Documentation updated if needed

## Troubleshooting

### "File too large" Error

```
HTTPError: 400 Bad Request from https://upload.pypi.org/legacy/
File too large. Limit for project 'CondenSimAdapter' is 100 MB
```

**Solution:** Wait for PyPI size limit increase approval, or see alternatives below.

### Alternatives if PyPI Rejects Size Increase

If PyPI doesn't approve the increase, we have backup plans:

#### Plan B: GitHub Releases + pip

Host models on GitHub Releases, download on first use:

```bash
pip install CondenSimAdapter  # Small package
adapter models download       # Download models from GitHub
```

This is already implemented in `CondenSimAdapter/backmap/cg2all/model/__init__.py`

#### Plan C: conda-forge

conda-forge has no strict size limits:

```bash
conda install -c conda-forge conden-sim-adapter
```

Can include all 200MB without issues.

## File Structure in Package

```
CondenSimAdapter-0.2.0-py3-none-any.whl (~200 MB)
├── CondenSimAdapter/
│   ├── __init__.py
│   ├── core/
│   ├── backmap/
│   │   └── cg2all/
│   │       └── model/
│   │           ├── __init__.py
│   │           ├── CalphaBasedModel.ckpt      (~43 MB)
│   │           ├── CalphaBasedModel-FIX.ckpt  (~43 MB)
│   │           ├── ResidueBasedModel.ckpt     (~43 MB)
│   │           ├── Martini.ckpt               (~30 MB)
│   │           └── Martini3.ckpt              (~30 MB)
│   └── ...
└── ...
```

## Contact

- PyPI support: https://github.com/pypi/support/issues
- Conda-forge: https://github.com/conda-forge/staged-recipes
