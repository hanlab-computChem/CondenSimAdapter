#!/bin/bash
# Create GitHub Release for model files
# Usage: ./scripts/create_github_release.sh [VERSION]
# Example: ./scripts/create_github_release.sh 0.1.0

set -e

VERSION=${1:-0.1.0}
REPO="hanlab-computChem/CondenSimAdapter"

echo "Creating GitHub Release v${VERSION} for models..."
echo "Repository: ${REPO}"

# Check if gh CLI is installed
if ! command -v gh &> /dev/null; then
    echo "Error: GitHub CLI (gh) is not installed."
    echo "Install from: https://cli.github.com/"
    exit 1
fi

# Check if logged in
if ! gh auth status &> /dev/null; then
    echo "Error: Not logged in to GitHub CLI."
    echo "Run: gh auth login"
    exit 1
fi

# Model files
MODEL_DIR="CondenSimAdapter/backmap/cg2all/model"
MODEL_FILES=(
    "${MODEL_DIR}/CalphaBasedModel.ckpt"
    "${MODEL_DIR}/CalphaBasedModel-FIX.ckpt"
    "${MODEL_DIR}/ResidueBasedModel.ckpt"
    "${MODEL_DIR}/Martini.ckpt"
    "${MODEL_DIR}/Martini3.ckpt"
)

# Check if all model files exist
echo ""
echo "Checking model files..."
MISSING=0
for file in "${MODEL_FILES[@]}"; do
    if [ -f "$file" ]; then
        SIZE=$(du -h "$file" | cut -f1)
        echo "  ✓ $(basename "$file") ($SIZE)"
    else
        echo "  ✗ $(basename "$file") NOT FOUND"
        MISSING=$((MISSING + 1))
    fi
done

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "Error: $MISSING model files are missing!"
    echo "Please ensure all models are in ${MODEL_DIR}/"
    exit 1
fi

echo ""
echo "Creating release v${VERSION}..."

# Create release (or use existing)
gh release create "v${VERSION}" \
    --title "Model Checkpoints v${VERSION}" \
    --notes "AI model checkpoints for cg2all backmapping.\n\nDownloaded automatically by the package on first use.\n\nModels included:\n- CalphaBasedModel.ckpt\n- CalphaBasedModel-FIX.ckpt\n- ResidueBasedModel.ckpt\n- Martini.ckpt\n- Martini3.ckpt" \
    || echo "Release v${VERSION} already exists, will upload assets..."

# Upload model files
echo ""
echo "Uploading model files..."
for file in "${MODEL_FILES[@]}"; do
    echo "  Uploading $(basename "$file")..."
    gh release upload "v${VERSION}" "$file" --clobber
done

echo ""
echo "✓ Release v${VERSION} created successfully!"
echo ""
echo "Release URL: https://github.com/${REPO}/releases/tag/v${VERSION}"
echo ""
echo "Users can now install the package and models will be downloaded automatically."
