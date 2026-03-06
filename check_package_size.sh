#!/bin/bash
# Check package size before publishing

echo "========================================"
echo "Package Size Checker"
echo "========================================"

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf build/ dist/ *.egg-info/

# Build package
echo "Building package..."
python -m build

echo ""
echo "========================================"
echo "Package Sizes"
echo "========================================"

# Get sizes
WHEEL=$(ls dist/*.whl 2>/dev/null | head -1)
SDIST=$(ls dist/*.tar.gz 2>/dev/null | head -1)

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
    echo "1. Request PyPI size limit increase:"
    echo "   https://github.com/pypi/support/issues/new?template=limit-request.md"
    echo ""
    echo "2. Use this template for the request:"
    cat << 'EOF'

**Project:** CondenSimAdapter  
**Size of release:** ~200 MB  
**Which indexes:** PyPI  

**Description:**  
CondenSimAdapter is a scientific software package for protein condensate simulations.  
It includes neural network models (~189 MB) for CG-to-AA backmapping.  

The large file size is necessary because:  
1. Neural network checkpoint files are inherently large  
2. Users expect a complete, working package out-of-the-box  
3. External hosting adds complexity for academic users  

The package is:  
- Open source (GPL-3.0)  
- Actively maintained  
- Used by the computational chemistry/biology community  

Please increase the size limit to 300 MB.  

EOF
else
    echo "✓ Package size is acceptable for PyPI"
fi
