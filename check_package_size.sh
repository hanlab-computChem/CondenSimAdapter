#!/bin/bash
# Check package size before publishing.
# Model checkpoints are downloaded at runtime and must not be included.

set -euo pipefail

echo "========================================"
echo "Package Size Checker"
echo "========================================"

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf build/ dist/ *.egg-info/

# Build package
echo "Building package..."
BUILD_ARGS=()
if [ "${CSA_BUILD_NO_ISOLATION:-0}" = "1" ]; then
    BUILD_ARGS+=(--no-isolation)
    BUILD_ARGS+=(--skip-dependency-check)
fi
python -m build "${BUILD_ARGS[@]}"

echo ""
echo "========================================"
echo "Package Sizes"
echo "========================================"

# Get sizes
WHEEL=$(ls dist/*.whl 2>/dev/null | head -1 || true)
SDIST=$(ls dist/*.tar.gz 2>/dev/null | head -1 || true)
WHEEL_SIZE=0

if [ -f "$WHEEL" ]; then
    WHEEL_SIZE=$(stat -c%s "$WHEEL")
    WHEEL_MB=$((WHEEL_SIZE / 1024 / 1024))
    echo "Wheel: $WHEEL"
    echo "Size: $WHEEL_SIZE bytes ($WHEEL_MB MB)"
    
    if [ $WHEEL_SIZE -gt 104857600 ]; then
        echo "⚠️  WARNING: Exceeds 100MB PyPI limit!"
        echo "   Need to request size limit increase."
    else
        echo "✓ Within PyPI 100MB limit"
    fi
fi

echo ""

if [ -f "$SDIST" ]; then
    SDIST_SIZE=$(stat -c%s "$SDIST")
    SDIST_MB=$((SDIST_SIZE / 1024 / 1024))
    echo "Sdist: $SDIST"
    echo "Size: $SDIST_SIZE bytes ($SDIST_MB MB)"
fi

echo ""
echo "========================================"
echo "Package Contents"
echo "========================================"

if [ -f "$WHEEL" ]; then
    BAD_WHEEL_CONTENT=$(python - "$WHEEL" << 'PY'
import sys
import zipfile

bad_suffixes = (".ckpt", ".pt", ".pth", ".pyc", ".pyo", ".o")
bad_names = []
with zipfile.ZipFile(sys.argv[1]) as zf:
    for name in zf.namelist():
        base = name.rsplit("/", 1)[-1]
        if name.endswith(bad_suffixes) or base == "genPairPACE" or "__pycache__" in name:
            bad_names.append(name)
print("\n".join(bad_names))
PY
)
    if [ -n "$BAD_WHEEL_CONTENT" ]; then
        echo "ERROR: Wheel contains generated or external model files:"
        echo "$BAD_WHEEL_CONTENT"
        exit 1
    fi
    echo "Wheel content check passed"
fi

echo ""
echo "========================================"
echo "Model Files"
echo "========================================"

MODEL_DIR="CondenSimAdapter/backmap/cg2all/model"
if [ -d "$MODEL_DIR" ]; then
    TOTAL_MODEL_SIZE=0
    for f in "$MODEL_DIR"/*.ckpt; do
        if [ -f "$f" ]; then
            SIZE=$(stat -c%s "$f")
            TOTAL_MODEL_SIZE=$((TOTAL_MODEL_SIZE + SIZE))
            echo "$(basename $f): $((SIZE / 1024 / 1024)) MB"
        fi
    done
    echo ""
    echo "Total model size: $((TOTAL_MODEL_SIZE / 1024 / 1024)) MB"
else
    echo "No model files found"
fi

echo ""
echo "========================================"
echo "Recommendations"
echo "========================================"

if [ -f "$WHEEL" ] && [ $WHEEL_SIZE -gt 104857600 ]; then
    echo "ERROR: Wheel exceeds 100MB PyPI limit."
    echo "Model checkpoints should be downloaded at runtime, not bundled."
    exit 1
else
    echo "Package size is acceptable for PyPI"
fi
