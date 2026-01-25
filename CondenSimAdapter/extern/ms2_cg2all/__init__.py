#!/usr/bin/env python
"""
ms2_cg2all - Convert coarse-grained protein structures to all-atom models.

This module provides functionality to convert single CG PDB files to all-atom PDB files
using the cg2all neural network model.

Usage:
    from CondenSimAdapter.extern.ms2_cg2all import convert_cg2all
    
    # Convert a CA-trace to all-atom structure
    convert_cg2all(
        in_pdb_fn="input.ca.pdb",
        out_fn="output.all.pdb",
        model_type="CalphaBasedModel"
    )
"""

from .lib.snippets import convert_cg2all

__all__ = ["convert_cg2all"]
