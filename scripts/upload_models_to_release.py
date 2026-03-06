#!/usr/bin/env python3
"""
Upload model files to GitHub Release using Python requests.
Usage: python upload_models_to_release.py [VERSION] [GITHUB_TOKEN]
"""

import os
import sys
import requests
from pathlib import Path

# Configuration
REPO = "hanlab-computChem/CondenSimAdapter"
MODEL_DIR = Path("CondenSimAdapter/backmap/cg2all/model")

# Models to upload (including FIX version)
MODELS = [
    "CalphaBasedModel.ckpt",
    "CalphaBasedModel-FIX.ckpt",
    "ResidueBasedModel.ckpt",
    "Martini.ckpt",
    "Martini3.ckpt",
]


def create_release(version: str, token: str):
    """Create a GitHub release."""
    url = f"https://api.github.com/repos/{REPO}/releases"
    headers = {
        "Authorization": f"token {token}",
        "Accept": "application/vnd.github.v3+json"
    }
    data = {
        "tag_name": f"v{version}",
        "name": f"Model Checkpoints v{version}",
        "body": """AI model checkpoints for cg2all backmapping.

Downloaded automatically by the package on first use.

Models included:
- CalphaBasedModel.ckpt
- CalphaBasedModel-FIX.ckpt (for fix_atom=True mode)
- ResidueBasedModel.ckpt
- Martini.ckpt
- Martini3.ckpt

Total size: ~235 MB
""",
        "draft": False,
        "prerelease": False
    }
    
    response = requests.post(url, headers=headers, json=data)
    if response.status_code == 201:
        print(f"✓ Release v{version} created")
        return response.json()["upload_url"].replace("{?name,label}", "")
    elif response.status_code == 422 and "already_exists" in str(response.json()):
        print(f"Release v{version} already exists, getting upload URL...")
        # Get existing release
        response = requests.get(url, headers=headers)
        for release in response.json():
            if release["tag_name"] == f"v{version}":
                return release["upload_url"].replace("{?name,label}", "")
        raise RuntimeError("Could not find existing release")
    else:
        raise RuntimeError(f"Failed to create release: {response.status_code} {response.text}")


def upload_asset(upload_url: str, token: str, file_path: Path):
    """Upload a file to the release."""
    filename = file_path.name
    
    # Check if asset already exists
    release_url = upload_url.rsplit("/", 1)[0]
    headers = {"Authorization": f"token {token}", "Accept": "application/vnd.github.v3+json"}
    response = requests.get(release_url, headers=headers)
    if response.status_code == 200:
        assets = response.json().get("assets", [])
        for asset in assets:
            if asset["name"] == filename:
                print(f"  {filename} already exists, deleting...")
                delete_url = f"https://api.github.com/repos/{REPO}/releases/assets/{asset['id']}"
                requests.delete(delete_url, headers=headers)
    
    # Upload file
    print(f"  Uploading {filename} ({file_path.stat().st_size / 1024 / 1024:.1f} MB)...")
    
    headers = {
        "Authorization": f"token {token}",
        "Content-Type": "application/octet-stream"
    }
    
    with open(file_path, "rb") as f:
        response = requests.post(
            f"{upload_url}?name={filename}",
            headers=headers,
            data=f
        )
    
    if response.status_code == 201:
        print(f"  ✓ {filename} uploaded")
        return True
    else:
        print(f"  ✗ Failed: {response.status_code} {response.text}")
        return False


def main():
    version = sys.argv[1] if len(sys.argv) > 1 else "0.1.0"
    token = sys.argv[2] if len(sys.argv) > 2 else os.environ.get("GITHUB_TOKEN")
    
    if not token:
        print("Error: GitHub token required")
        print("Usage: python upload_models_to_release.py [VERSION] [GITHUB_TOKEN]")
        print("Or set GITHUB_TOKEN environment variable")
        sys.exit(1)
    
    print(f"Creating release v{version} for {REPO}")
    print()
    
    # Create release
    upload_url = create_release(version, token)
    print()
    
    # Upload models
    success_count = 0
    for model in MODELS:
        file_path = MODEL_DIR / model
        if not file_path.exists():
            print(f"  ✗ {model} not found at {file_path}")
            continue
        
        if upload_asset(upload_url, token, file_path):
            success_count += 1
        print()
    
    print(f"✓ Done: {success_count}/{len(MODELS)} models uploaded")
    print(f"Release URL: https://github.com/{REPO}/releases/tag/v{version}")


if __name__ == "__main__":
    main()
