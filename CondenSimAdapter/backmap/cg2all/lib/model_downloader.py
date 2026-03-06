#!/usr/bin/env python
"""
Model downloader for cg2all checkpoints.

Downloads models from Hugging Face Hub if not present locally.
"""

import os
from pathlib import Path
from typing import Optional

try:
    from huggingface_hub import hf_hub_download, list_repo_files
    HF_AVAILABLE = True
except ImportError:
    HF_AVAILABLE = False

# HuggingFace repository containing the models
HF_REPO_ID = "hanlab/condensimadapter-cg2all-models"

# Model file mapping
MODEL_FILES = {
    "CalphaBasedModel": "CalphaBasedModel.ckpt",
    "CalphaBasedModel-FIX": "CalphaBasedModel-FIX.ckpt",
    "ResidueBasedModel": "ResidueBasedModel.ckpt",
    "Martini": "Martini.ckpt",
    "Martini3": "Martini3.ckpt",
}


def download_model(model_name: str, model_home: Path, force: bool = False) -> Path:
    """
    Download a model from Hugging Face Hub if not present locally.
    
    Args:
        model_name: Name of the model (e.g., "CalphaBasedModel")
        model_home: Local directory to store models
        force: If True, re-download even if file exists
        
    Returns:
        Path to the local model file
        
    Raises:
        FileNotFoundError: If model cannot be downloaded
        ImportError: If huggingface_hub is not installed
    """
    if not HF_AVAILABLE:
        raise ImportError(
            "huggingface_hub is required to download models.\n"
            "Install with: pip install huggingface_hub\n"
            "Or manually download models from: https://huggingface.co/hanlab/condensimadapter-cg2all-models"
        )
    
    model_filename = MODEL_FILES.get(model_name)
    if model_filename is None:
        raise ValueError(f"Unknown model: {model_name}. Available: {list(MODEL_FILES.keys())}")
    
    local_path = model_home / model_filename
    
    # Check if already exists and is valid
    if not force and local_path.exists():
        file_size = local_path.stat().st_size
        if file_size > 1024 * 1024:  # > 1MB indicates valid model
            return local_path
    
    # Download from Hugging Face
    print(f"Downloading model '{model_name}' from Hugging Face Hub...")
    try:
        downloaded_path = hf_hub_download(
            repo_id=HF_REPO_ID,
            filename=model_filename,
            local_dir=model_home,
            local_dir_use_symlinks=False,
        )
        print(f"Model downloaded successfully to: {downloaded_path}")
        return Path(downloaded_path)
    except Exception as e:
        raise FileNotFoundError(
            f"Failed to download model '{model_name}' from Hugging Face Hub.\n"
            f"Repository: https://huggingface.co/{HF_REPO_ID}\n"
            f"Error: {e}"
        )


def ensure_model_available(model_name: str, model_home: Path) -> Path:
    """
    Ensure a model is available locally, downloading if necessary.
    
    Args:
        model_name: Name of the model
        model_home: Local directory for models
        
    Returns:
        Path to the local model file
    """
    model_filename = MODEL_FILES.get(model_name, f"{model_name}.ckpt")
    local_path = model_home / model_filename
    
    if local_path.exists() and local_path.stat().st_size > 1024 * 1024:
        return local_path
    
    # Try to download
    return download_model(model_name, model_home)
