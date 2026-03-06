#!/bin/bash
# ===========================================
# CondenSimAdapter - pip Installation Script
# For users who prefer pip over conda
# ===========================================

set -e

PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
INSTALL_TYPE="${1:-cpu}"  # cpu or gpu

echo "============================================"
echo "CondenSimAdapter - pip Installation"
echo "============================================"
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

# Step 1: Install PyTorch
echo ""
echo "📦 Step 1/4: Installing PyTorch ($INSTALL_TYPE)..."
if [ "$INSTALL_TYPE" = "gpu" ]; then
    pip install torch==2.1.2 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
    pip install dgl -f https://data.dgl.ai/wheels/torch-2.1/cu121/repo.html
else
    pip install torch==2.1.2 torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
    pip install dgl -f https://data.dgl.ai/wheels/torch-2.1/repo.html
fi

# Step 2: Install OpenMM and scientific packages
echo ""
echo "📦 Step 2/4: Installing scientific packages..."
pip install \
    openmm>=8.2 \
    mdtraj>=1.10 \
    MDAnalysis>=2.6 \
    biopython>=1.81 \
    parmed \
    gromacswrapper \
    scipy>=1.10 \
    matplotlib>=3.5 \
    networkx>=2.8 \
    numba>=0.60 \
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
