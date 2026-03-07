#!/bin/bash
# ===========================================
# CondenSimAdapter - pip Installation Script
# For users who prefer pip over conda
# ===========================================

set -e

PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
INSTALL_TYPE="${1:-cpu}"  # cpu, gpu-cuda12, gpu-cuda13, gpu-hip6, gpu-hip7

echo "============================================"
echo "CondenSimAdapter - pip Installation"
echo "============================================"
echo ""
echo "Usage: bash install-pip.sh [INSTALL_TYPE]"
echo "  INSTALL_TYPE options:"
echo "    cpu          - CPU only (default)"
echo "    gpu          - NVIDIA GPU, CUDA 12.x (alias for gpu-cuda12)"
echo "    gpu-cuda12   - NVIDIA GPU, CUDA 12.x (recommended)"
echo "    gpu-cuda13   - NVIDIA GPU, CUDA 13.x"
echo "    gpu-hip6     - AMD GPU, ROCm/HIP 6"
echo "    gpu-hip7     - AMD GPU, ROCm/HIP 7"
echo ""

# Check Python version
CURRENT_PY=$(python -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
echo "Python version: $CURRENT_PY"

if [[ "$CURRENT_PY" != "$PYTHON_VERSION"* ]]; then
    echo "⚠️  Warning: Python $PYTHON_VERSION recommended, found $CURRENT_PY"
    read -p "Continue? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Step 1: Install PyTorch + DGL (CUDA 12.x ecosystem, targeting CUDA 12.4)
# dgl 2.1.0 is the last PyPI release and is API-compatible with this package.
echo ""
echo "📦 Step 1/4: Installing PyTorch + DGL ($INSTALL_TYPE)..."
if [ "$INSTALL_TYPE" = "gpu-cuda12" ] || [ "$INSTALL_TYPE" = "gpu" ]; then
    # CUDA 12.4 wheels; backward-compatible with CUDA 12.0-12.4 drivers
    pip install "torch>=2.4,<2.5" torchvision torchaudio \
        --index-url https://download.pytorch.org/whl/cu124
    pip install "dgl==2.1.0"
elif [ "$INSTALL_TYPE" = "gpu-cuda13" ]; then
    pip install "torch>=2.4,<2.5" torchvision torchaudio \
        --index-url https://download.pytorch.org/whl/cu124
    echo "  Note: using CUDA 12.4 wheels for torch; openmm will use cuda13 platform."
    pip install "dgl==2.1.0"
else
    pip install "torch>=2.4,<2.5" torchvision torchaudio \
        --index-url https://download.pytorch.org/whl/cpu
    pip install "dgl==2.1.0"
fi

# Step 2: Install OpenMM and scientific packages
# CondenSimAdapter depends on openmm[cuda12] by default.
# CPU-only users must override by installing plain openmm first, then the package with --no-deps.
echo ""
echo "📦 Step 2/4: Installing scientific packages..."
if [ "$INSTALL_TYPE" = "gpu-cuda13" ]; then
    OPENMM_PKG="openmm[cuda13]>=8.2"
elif [ "$INSTALL_TYPE" = "gpu-hip6" ]; then
    OPENMM_PKG="openmm[hip6]>=8.2"
elif [ "$INSTALL_TYPE" = "gpu-hip7" ]; then
    OPENMM_PKG="openmm[hip7]>=8.2"
elif [ "$INSTALL_TYPE" = "cpu" ]; then
    OPENMM_PKG="openmm>=8.2"
else
    # default: cuda12 (covers gpu, gpu-cuda12)
    OPENMM_PKG="openmm[cuda12]>=8.2"
fi

pip install \
    "${OPENMM_PKG}" \
    "mdtraj>=1.10" \
    "MDAnalysis>=2.6" \
    "biopython>=1.81" \
    parmed \
    gromacswrapper \
    "scipy>=1.10" \
    "matplotlib>=3.5" \
    "networkx>=2.8" \
    "numba>=0.60" \
    jinja2 \
    statsmodels \
    PeptideConstructor

# Step 3: Install CondenSimAdapter
echo ""
echo "📦 Step 3/4: Installing CondenSimAdapter..."
pip install CondenSimAdapter

# Step 4: Verify
echo ""
echo "📦 Step 4/4: Verifying installation..."
adapter --version
adapter --help

echo ""
echo "============================================"
echo "✅ Installation Complete!"
echo "============================================"
echo ""
echo "Quick start:"
echo "  adapter init my_project"
echo "  adapter cg -f config.yaml"
echo ""
echo "Check model status:"
echo "  adapter models status"
echo ""

# Check if models are included
python -c "
from CondenSimAdapter.backmap.cg2all.model import MODEL_PATHS
import sys

models = list(MODEL_PATHS.keys())
print('Models in package:', len(models))
for m in models:
    path = MODEL_PATHS[m]
    if path.exists():
        size_mb = path.stat().st_size / (1024*1024)
        print(f'  ✓ {m}: {size_mb:.1f} MB')
    else:
        print(f'  ✗ {m}: Not found')
        sys.exit(1)

print('')
print('✓ All models ready!')
" || {
    echo ""
    echo "⚠️  Some models may be missing."
    echo "   Run 'adapter models download' to download them."
}
