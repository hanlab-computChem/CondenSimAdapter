"""
ms2_openabc - Internalized OpenABC package

OpenABC (OpenMM GPU-Accelerated simulations of Biomolecular Condensates)
is a flexible package that implements multiple popular coarse-grained force fields
for simulations, including:
- HPS (Hydropathy Scale) model
- MOFF (Maximum Entropy Optimized Force Field) model
- Mpipi model
- SMOG+3SPN2 model

This is an internalized version of OpenABC for use within CondenSimAdapter.
Original package: https://github.com/ZhangGroup-MITChemistry/OpenABC
"""

__version__ = '1.0.7'

# Expose main modules
from . import forcefields
from . import utils
__all__ = [
    'forcefields',
    'utils',
    'lib',
    '__version__',
]
