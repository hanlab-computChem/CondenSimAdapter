#!/bin/bash
# ===========================================
# Multiscale² install script
# ===========================================

set -e  # Exit immediately on error

echo "==========================================="
echo "Multiscale² environment install script"
echo "==========================================="

# Check Python version (prefer current environment python)
PYTHON_BIN=${PYTHON_BIN:-python}
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    echo "❌ Python executable not found: ${PYTHON_BIN}"
    exit 1
fi
PYTHON_VERSION=$(${PYTHON_BIN} -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
echo "[1/6] Current Python version: ${PYTHON_VERSION}"

# Local torch wheel path
TORCH_WHEEL_PATH="${TORCH_WHEEL_PATH:-$(pwd)/Vender_packages/torch-2.1.2+cpu-cp311-cp311-linux_x86_64.whl}"

if [ "$PYTHON_VERSION" != "3.11" ]; then
    echo "⚠️  Warning: Python 3.11 is recommended, current version is ${PYTHON_VERSION}"
    echo "Continuing may cause compatibility issues."
    read -p "Continue? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Installation cancelled."
        exit 1
    fi
fi

# Create conda env (if not activated)
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "[2/6] No active conda env detected, creating CondenSimAdapter..."
    conda create -n CondenSimAdapter python=3.11 -y
    echo "Please run: conda activate CondenSimAdapter"
    exit 0
fi

if [ "$CONDA_DEFAULT_ENV" != "CondenSimAdapter" ]; then
    echo "[2/6] Current env: ${CONDA_DEFAULT_ENV} (not CondenSimAdapter)"
    echo "Installation will continue in the current env. Create CondenSimAdapter if needed."
fi

# Ensure mamba is available (faster conda installs)
if ! command -v mamba >/dev/null 2>&1; then
    echo "[2/6] Installing mamba..."
    conda install -c conda-forge -y mamba
fi

# Force CUDA version for conda packages
CUDA_VERSION=12.1

DRY_RUN=${DRY_RUN:-0}
if [ "${DRY_RUN}" = "1" ]; then
    cat <<'EOF'
Manual install commands (same as script)

# PyTorch / DGL (CPU)
pip install torch==2.1.2 torchvision==0.16.2 torchaudio==2.1.2 --index-url https://download.pytorch.org/whl/cpu
pip install numpy==1.26.4 scipy==1.14.1 click pydantic
pip install dgl==2.1.0 -f https://data.dgl.ai/wheels/torch-2.1/repo.html
pip install e3nn==0.5.1 torchdata==0.7.1 ml-collections==0.1.1
pip install git+https://github.com/huhlim/SE3Transformer.git

# Conda packages (OpenMM >= 8.2.0; CUDA pinned)
mamba install -c conda-forge -y "openmm>=8.2.0" "cuda-version=${CUDA_VERSION}" matplotlib networkx gromacswrapper parmed
pip install git+https://github.com/feiglab/mdsim.git

# Python packages
pip install mdanalysis==2.6.1 biopython==1.81 numba==0.60.0 tqdm pyyaml jinja2 localcider statsmodels
pip install mdtraj==1.10.0
EOF
    exit 0
fi

echo "[3/6] Installing PyTorch stack (CPU)..."

# Decide whether to use local torch wheel
USE_LOCAL_TORCH=0
if [ "$PYTHON_VERSION" = "3.11" ] && [ -f "$TORCH_WHEEL_PATH" ]; then
    echo "  Using local torch wheel: ${TORCH_WHEEL_PATH}"
    USE_LOCAL_TORCH=1
elif [ "$PYTHON_VERSION" = "3.11" ] && [ ! -f "$TORCH_WHEEL_PATH" ]; then
    echo "  ⚠️  Local torch wheel not found: ${TORCH_WHEEL_PATH}"
    echo "  Falling back to online install..."
fi

if [ "$USE_LOCAL_TORCH" = "1" ]; then
    # Install torch locally
    pip install "${TORCH_WHEEL_PATH}"
else
    # Install torch from index
    pip install \
        torch==2.1.2 \
        torchvision==0.16.2 \
        torchaudio==2.1.2 \
        --index-url https://download.pytorch.org/whl/cpu
fi

pip install numpy==1.26.4 scipy==1.14.1 click pydantic

pip install dgl==2.1.0 -f https://data.dgl.ai/wheels/torch-2.1/repo.html

pip install \
    e3nn==0.5.1 \
    torchdata==0.7.1 \
    ml-collections==0.1.1


echo "[4/6] Installing conda packages..."
mamba install -c conda-forge -y \
    "openmm>=8.2.0" \
    "cuda-version=${CUDA_VERSION}" \
    matplotlib \
    networkx \
    gromacswrapper \
    parmed 

#pip install git+https://github.com/feiglab/mdsim.git


echo "[5/6] Installing Python packages..."
pip install \
    mdanalysis==2.6.1 \
    biopython==1.81 \
    numba==0.60.0 \
    tqdm \
    pyyaml \
    jinja2 \
    localcider\
    statsmodels

echo "[6/6] Installing mdtraj..."
pip install mdtraj==1.10.0
pip install pytest
pip install PeptideConstructor


echo ""
echo "==========================================="
echo "✅ Installation complete!"
echo "==========================================="
echo ""
echo "Verify installation:"
echo "  python -c \"import torch; print(f'PyTorch: {torch.__version__}, CUDA: {torch.version.cuda}')\""
echo "  python -c \"import openmm; print(f'OpenMM: {openmm.__version__}')\""
echo "  python -c \"import MDAnalysis; print(f'MDAnalysis: {MDAnalysis.__version__}')\""
echo ""
echo "If you run into issues, run:"
echo "  python scripts/verify_env.py"

