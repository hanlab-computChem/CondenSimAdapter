#!/bin/bash
# ===========================================
# CondenSimAdapter - Conda Installation Script
# One-line setup: conda create -n csa && conda install -c conda-forge conden-sim-adapter
# ===========================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${1:-csa}"
PYTHON_VERSION="${PYTHON_VERSION:-3.11}"

echo "============================================"
echo "CondenSimAdapter - Conda Installation"
echo "============================================"

# Check conda
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install miniconda:"
    echo "   https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "📦 Step 1/4: Creating environment '${ENV_NAME}'..."
conda create -n "${ENV_NAME}" python=${PYTHON_VERSION} -y

echo "🚀 Step 2/4: Activating environment..."
echo "   Run: conda activate ${ENV_NAME}"
echo ""
echo "   Then install with:"
echo "   conda install -c conda-forge conden-sim-adapter"
echo ""

# Create activation script
cat > "${SCRIPT_DIR}/activate-${ENV_NAME}.sh" << EOF
#!/bin/bash
conda activate ${ENV_NAME}
echo "Environment '${ENV_NAME}' activated!"
echo "Test with: adapter --help"
EOF
chmod +x "${SCRIPT_DIR}/activate-${ENV_NAME}.sh"

echo "✅ Step 3/4: Activation script created: activate-${ENV_NAME}.sh"

echo ""
echo "============================================"
echo "Next Steps:"
echo "============================================"
echo "1. Activate environment:"
echo "   conda activate ${ENV_NAME}"
echo ""
echo "2. Install package:"
echo "   conda install -c conda-forge conden-sim-adapter"
echo ""
echo "   Or build locally:"
echo "   cd conda-forge-recipe && conda build ."
echo "   conda install --use-local conden-sim-adapter"
echo ""
echo "3. Test installation:"
echo "   adapter --help"
echo "   adapter init my_project"
echo ""
echo "============================================"
