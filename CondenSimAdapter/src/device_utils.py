#!/usr/bin/env python3
"""
Device utilities for GPU/CPU selection.
"""

from typing import Dict, Optional, Tuple

SUPPORTED_DEVICES = ("cpu", "cuda", "opencl")


def normalize_device(device: Optional[str], default: str = "cpu") -> str:
    """Normalize device string and apply default."""
    if not device:
        device = default
    device = device.strip().lower()
    if device not in SUPPORTED_DEVICES:
        raise ValueError(f"Unsupported device: {device}. Supported: {SUPPORTED_DEVICES}")
    return device


def resolve_torch_device(device: str, gpu_id: Optional[int]) -> str:
    """Resolve torch device string with optional GPU index."""
    device = normalize_device(device, default="cpu")
    if device == "cuda" and gpu_id is not None:
        return f"cuda:{gpu_id}"
    return device


def resolve_openmm_platform(
    device: str,
    gpu_id: Optional[int],
    precision: str = "mixed",
) -> Tuple[str, Dict[str, str]]:
    """Resolve OpenMM platform name and properties."""
    device = normalize_device(device, default="cpu")
    platform_name = device.upper()
    properties: Dict[str, str] = {}

    if platform_name in ("CUDA", "OPENCL"):
        if gpu_id is None:
            gpu_id = 0
        properties = {
            "Precision": precision,
            "DeviceIndex": str(gpu_id),
        }

    return platform_name, properties
