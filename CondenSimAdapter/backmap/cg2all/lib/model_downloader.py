#!/usr/bin/env python
"""
Model downloader for cg2all checkpoints.

Downloads models from GitHub Release if not present locally.
"""

import os
import urllib.request
import urllib.error
from pathlib import Path
from typing import Optional
import sys

# GitHub Release URL template
# Format: https://github.com/{OWNER}/{REPO}/releases/download/{TAG}/{FILENAME}
GITHUB_RELEASE_URL = "https://github.com/hanlab-computChem/CondenSimAdapter/releases/download/v{version}/{filename}"

# Default version to download
DEFAULT_VERSION = "1.0.0-beta"

# Model file mapping
MODEL_FILES = {
    "CalphaBasedModel": "CalphaBasedModel.ckpt",
    "CalphaBasedModel-FIX": "CalphaBasedModel-FIX.ckpt",
    "ResidueBasedModel": "ResidueBasedModel.ckpt",
    "Martini": "Martini.ckpt",
    "Martini3": "Martini3.ckpt",
}

# Environment variable to override version
VERSION_ENV_VAR = "CONDENSIMADAPTER_MODEL_VERSION"


def _get_version() -> str:
    """Get the model version from environment or default."""
    return os.environ.get(VERSION_ENV_VAR, DEFAULT_VERSION)


def _download_file(url: str, dest_path: Path, show_progress: bool = True) -> None:
    """Download a file from URL to destination with progress bar."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    
    def report_hook(block_num, block_size, total_size):
        if show_progress and total_size > 0:
            downloaded = block_num * block_size
            percent = min(100, downloaded * 100 / total_size)
            sys.stdout.write(f"\r  Downloading: {percent:.1f}%")
            sys.stdout.flush()
    
    try:
        urllib.request.urlretrieve(url, dest_path, reporthook=report_hook)
        if show_progress:
            sys.stdout.write("\n")
            sys.stdout.flush()
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise FileNotFoundError(
                f"Model not found at {url}\n"
                f"Please ensure the release exists and contains the model files.\n"
                f"You can manually download models from: https://github.com/hanlab-computChem/CondenSimAdapter/releases"
            )
        raise


def download_model(model_name: str, model_home: Path, version: Optional[str] = None, force: bool = False) -> Path:
    """
    Download a model from GitHub Release if not present locally.
    
    Args:
        model_name: Name of the model (e.g., "CalphaBasedModel")
        model_home: Local directory to store models
        version: Release version to download (default: 0.1.0)
        force: If True, re-download even if file exists
        
    Returns:
        Path to the local model file
        
    Raises:
        FileNotFoundError: If model cannot be downloaded
    """
    model_filename = MODEL_FILES.get(model_name)
    if model_filename is None:
        raise ValueError(f"Unknown model: {model_name}. Available: {list(MODEL_FILES.keys())}")
    
    local_path = model_home / model_filename
    
    # Check if already exists and is valid
    if not force and local_path.exists():
        file_size = local_path.stat().st_size
        if file_size > 1024 * 1024:  # > 1MB indicates valid model
            return local_path
    
    # Determine version
    if version is None:
        version = _get_version()
    
    # Construct download URL
    url = GITHUB_RELEASE_URL.format(version=version, filename=model_filename)
    
    # Download
    print(f"Downloading model '{model_name}' (v{version}) from GitHub Release...")
    print(f"  URL: {url}")
    
    try:
        _download_file(url, local_path)
        print(f"  Saved to: {local_path}")
        return local_path
    except Exception as e:
        # Clean up partial download
        if local_path.exists():
            local_path.unlink()
        raise FileNotFoundError(
            f"Failed to download model '{model_name}'.\n"
            f"Error: {e}\n\n"
            f"You can manually download models from:\n"
            f"https://github.com/hanlab-computChem/CondenSimAdapter/releases\n\n"
            f"And place them in: {model_home}\n\n"
            f"Or set environment variable to use a different version:\n"
            f"  export {VERSION_ENV_VAR}=0.2.0"
        )


def ensure_model_available(model_name: str, model_home: Path, version: Optional[str] = None) -> Path:
    """
    Ensure a model is available locally, downloading if necessary.
    
    Args:
        model_name: Name of the model
        model_home: Local directory for models
        version: Release version (optional)
        
    Returns:
        Path to the local model file
    """
    model_filename = MODEL_FILES.get(model_name, f"{model_name}.ckpt")
    local_path = model_home / model_filename
    
    if local_path.exists() and local_path.stat().st_size > 1024 * 1024:
        return local_path
    
    # Try to download
    return download_model(model_name, model_home, version)


def list_available_models(model_home: Path) -> dict:
    """
    List all available models and their status.
    
    Args:
        model_home: Local directory for models
        
    Returns:
        Dict mapping model name to status dict with 'exists' and 'size' keys
    """
    result = {}
    for model_name, filename in MODEL_FILES.items():
        path = model_home / filename
        if path.exists():
            size_mb = path.stat().st_size / (1024 * 1024)
            result[model_name] = {"exists": True, "size_mb": round(size_mb, 2), "path": str(path)}
        else:
            result[model_name] = {"exists": False, "size_mb": 0, "path": str(path)}
    return result
