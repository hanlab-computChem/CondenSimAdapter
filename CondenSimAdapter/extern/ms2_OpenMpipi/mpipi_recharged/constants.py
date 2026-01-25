"""
Constants and platform configuration for OpenMpipi.
"""

import json
import os
import numpy as np
import openmm as mm

# Get the package base directory for data files
PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_files')

# Platform and precision configuration
try:
    PLATFORM = mm.Platform.getPlatformByName('CUDA')
    PROPERTIES = {'Precision': 'mixed'}
    print('CUDA platform is available. Using CUDA platform with mixed precision.')
except Exception as e:
    PLATFORM = None
    PROPERTIES = None
    print(f'CUDA platform is not available. Using the default platform. Error: {e}')

