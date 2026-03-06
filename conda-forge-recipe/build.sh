#!/bin/bash

set -e

# Build the package
$PYTHON -m pip install . --no-deps --no-build-isolation -vv

# Verify installation
adapter --help || true

echo "Build completed successfully!"
