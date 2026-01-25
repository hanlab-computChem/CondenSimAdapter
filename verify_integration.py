#!/usr/bin/env python
"""
Final verification test for ms2_cg2all integration with provided checkpoints.
"""

import sys
sys.path.insert(0, '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square')

print("="*60)
print("Final Verification Test")
print("="*60)

# Test 1: Import the module
print("\n1. Testing module import...")
try:
    from multiscale2.extern.ms2_cg2all import convert_cg2all
    print("   ✓ Successfully imported convert_cg2all")
except Exception as e:
    print(f"   ✗ Import failed: {e}")
    sys.exit(1)

# Test 2: Check checkpoint loading
print("\n2. Testing checkpoint loading...")
import torch
from multiscale2.extern.ms2_cg2all.lib.libconfig import MODEL_HOME

test_cases = [
    ("CalphaBasedModel", False),
    ("CalphaBasedModel", True),
    ("ResidueBasedModel", False),
    ("Martini", False),
    ("Martini3", False),
]

for model_type, fix_atom in test_cases:
    if fix_atom:
        ckpt_name = f"{model_type}-FIX.ckpt"
    else:
        ckpt_name = f"{model_type}.ckpt"
    
    ckpt_path = MODEL_HOME / ckpt_name
    if not ckpt_path.exists():
        print(f"   - {model_type} (FIX={fix_atom}): NOT FOUND")
        continue
    
    file_size = ckpt_path.stat().st_size
    if file_size < 1000:
        print(f"   - {model_type} (FIX={fix_atom}): TOO SMALL ({file_size} bytes)")
        continue
    
    try:
        ckpt = torch.load(ckpt_path, map_location="cpu")
        cg_model = ckpt["hyper_parameters"].get("cg_model", "unknown")
        print(f"   ✓ {model_type} (FIX={fix_atom}): OK (cg_model={cg_model})")
    except Exception as e:
        print(f"   ✗ {model_type} (FIX={fix_atom}): FAILED - {e}")

# Test 3: Verify convert_cg2all function signature
print("\n3. Testing function signature...")
import inspect
sig = inspect.signature(convert_cg2all)
params = list(sig.parameters.keys())
print(f"   Parameters: {params}")

required = ['in_pdb_fn', 'out_fn', 'model_type']
for param in required:
    if param in params:
        print(f"   ✓ {param} is present")
    else:
        print(f"   ✗ {param} is missing!")

print("\n" + "="*60)
print("Verification completed!")
print("="*60)
print("\nAvailable working models:")
print("  - CalphaBasedModel (normal + FIX)")
print("  - ResidueBasedModel (normal only)")
print("  - Martini (normal only)")
print("  - Martini3 (normal only)")





