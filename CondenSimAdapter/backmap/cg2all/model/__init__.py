"""
Model management for cg2all backmapping.

Model files (~189MB total) are downloaded on first use to keep the pip package small.
Models are hosted on GitHub Releases.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Dict
import urllib.request
import hashlib
from tqdm import tqdm


MODEL_HOME = Path(__file__).parent

# Model download URLs (GitHub Releases)
# Replace with your actual release URLs after publishing
MODEL_URLS: Dict[str, str] = {
    "CalphaBasedModel.ckpt": "https://github.com/hanlab-computChem/CondenSimAdapter/releases/download/v1.0.0-beta/CalphaBasedModel.ckpt",
    "CalphaBasedModel-FIX.ckpt": "https://github.com/hanlab-computChem/CondenSimAdapter/releases/download/v1.0.0-beta/CalphaBasedModel-FIX.ckpt",
    "ResidueBasedModel.ckpt": "https://github.com/hanlab-computChem/CondenSimAdapter/releases/download/v1.0.0-beta/ResidueBasedModel.ckpt",
    "Martini.ckpt": "https://github.com/hanlab-computChem/CondenSimAdapter/releases/download/v1.0.0-beta/Martini.ckpt",
    "Martini3.ckpt": "https://github.com/hanlab-computChem/CondenSimAdapter/releases/download/v1.0.0-beta/Martini3.ckpt",
}

# Expected file sizes (bytes) for verification
MODEL_SIZES: Dict[str, int] = {
    "CalphaBasedModel.ckpt": 49_300_000,      # ~47 MB
    "CalphaBasedModel-FIX.ckpt": 49_300_000,  # ~47 MB
    "ResidueBasedModel.ckpt": 49_300_000,     # ~47 MB
    "Martini.ckpt": 49_300_000,               # ~47 MB
    "Martini3.ckpt": 49_300_000,              # ~47 MB
}


class DownloadProgressBar(tqdm):
    """Progress bar for downloads."""
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_file(url: str, output_path: Path, desc: str = "Downloading") -> None:
    """Download file with progress bar."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc=desc) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def get_model_path(model_name: str, auto_download: bool = True) -> Path:
    """Get path to model file, downloading if necessary.
    
    Args:
        model_name: Name of the model file (e.g., "CalphaBasedModel.ckpt")
        auto_download: Whether to download if missing
        
    Returns:
        Path to the model file
        
    Raises:
        FileNotFoundError: If model not found and auto_download=False or download fails
    """
    model_path = MODEL_HOME / model_name
    
    if model_path.exists():
        # Verify file size
        expected_size = MODEL_SIZES.get(model_name)
        actual_size = model_path.stat().st_size
        if expected_size and abs(actual_size - expected_size) > 1_000_000:  # 1MB tolerance
            print(f"Warning: {model_name} appears corrupted (size mismatch). Re-downloading...")
            model_path.unlink()
        else:
            return model_path
    
    if not auto_download:
        raise FileNotFoundError(
            f"Model file not found: {model_path}\n"
            f"Run 'adapter models download' to download models, or\n"
            f"download manually from: {MODEL_URLS.get(model_name, 'GitHub Releases')}"
        )
    
    # Download
    if model_name not in MODEL_URLS:
        raise FileNotFoundError(
            f"Unknown model: {model_name}. Available: {list(MODEL_URLS.keys())}"
        )
    
    url = MODEL_URLS[model_name]
    print(f"Downloading {model_name}...")
    print(f"URL: {url}")
    
    try:
        download_file(url, model_path, desc=model_name)
        print(f"✓ Saved to {model_path}")
        return model_path
    except Exception as e:
        if model_path.exists():
            model_path.unlink()  # Clean up partial download
        raise FileNotFoundError(f"Failed to download {model_name}: {e}")


def list_models() -> Dict[str, Dict]:
    """List all models and their status."""
    result = {}
    for name, url in MODEL_URLS.items():
        path = MODEL_HOME / name
        if path.exists():
            size_mb = path.stat().st_size / (1024 * 1024)
            status = "✓ Downloaded"
        else:
            size_mb = MODEL_SIZES.get(name, 0) / (1024 * 1024)
            status = "○ Not downloaded"
        
        result[name] = {
            "status": status,
            "size_mb": f"{size_mb:.1f}",
            "path": str(path) if path.exists() else None,
        }
    return result


def ensure_models(models: Optional[list[str]] = None) -> None:
    """Ensure specified models are downloaded. If None, download all.
    
    Args:
        models: List of model names to download, or None for all
    """
    if models is None:
        models = list(MODEL_URLS.keys())
    
    for model in models:
        get_model_path(model, auto_download=True)
    
    print("\n✓ All models ready!")


def print_model_status() -> None:
    """Print formatted model status table."""
    print("\n" + "="*60)
    print("Available Models")
    print("="*60)
    
    models = list_models()
    total_downloaded = 0
    total_size = 0
    
    for name, info in models.items():
        size = float(info["size_mb"])
        total_size += size
        if "✓" in info["status"]:
            total_downloaded += size
        print(f"{info['status']:15} {name:30} {info['size_mb']:>8} MB")
    
    print("-"*60)
    print(f"Downloaded: {total_downloaded:.1f} / {total_size:.1f} MB")
    print("="*60)


# For backwards compatibility - paths may not exist until downloaded
MODEL_PATHS: Dict[str, Path] = {
    name: MODEL_HOME / name for name in MODEL_URLS.keys()
}


def get_cached_model_path(model_name: str) -> Optional[Path]:
    """Get path only if model is already downloaded."""
    path = MODEL_HOME / model_name
    return path if path.exists() else None
