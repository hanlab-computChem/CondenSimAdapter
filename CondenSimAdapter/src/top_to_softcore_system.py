"""
top_to_softcore_system.py: Used for loading Gromacs top files with three-stage optimization.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2024 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Jason Swails

Extended with three-stage optimization support:
- NONBONDED_STANDARD: Standard mode matching native gromacstopfile.py behavior
- NONBONDED_GAUSSIAN: Gaussian repulsion for geometric untangling
- NONBONDED_SOFTCORE: Softcore potential for free energy calculations (placeholder)
"""
from __future__ import absolute_import
__author__ = "Peter Eastman"
__version__ = "1.0"

from openmm.app import Topology
from openmm.app import PDBFile
from openmm.app import forcefield as ff
from openmm.app import element as elem
from openmm.app import amberprmtopfile as prmtop
from openmm.app.internal import amber_file_parser
from openmm.app.internal.customgbforces import GBSAGBn2Force
from openmm.app import gromacstopfile
import openmm.unit as unit
import openmm as mm
import math
import os
import re
import shutil
from collections import OrderedDict, defaultdict
from itertools import combinations, combinations_with_replacement
from copy import deepcopy

HBonds = ff.HBonds
AllBonds = ff.AllBonds
HAngles = ff.HAngles

OBC2 = prmtop.OBC2

# Three-stage optimization types
NONBONDED_GAUSSIAN = "gaussian"
NONBONDED_SOFTCORE = "softcore"
NONBONDED_STANDARD = "standard"
NONBONDED_STANDARD_CUSTOM = "standard_custom"  # Standard LJ+Coulomb using CustomNonbondedForce (for debugging)



# Implicit solvent models
IMPLICIT_NONE = None
IMPLICIT_OBC1 = "OBC1"
IMPLICIT_OBC2 = "OBC2"
IMPLICIT_GBN = "GBn"
IMPLICIT_GBN2 = "GBn2"

# AMBER99SB-ILDN atom type mapping for GBSA (Born radii in nm)
# Based on MBondi2/MBondi3 parameterization
AMBER99SB_ILDN_ATOM_MAPPING = {
    # Hydrogen atoms
    'H':   {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Amide hydrogen
    'H1':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Aliphatic hydrogen
    'HP':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Aromatic hydrogen
    'HC':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Alkyl hydrogen
    'HB1': {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Beta hydrogen (MET)
    'HB2': {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Beta hydrogen (MET)
    'HS':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Thiol hydrogen
    'HO':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Hydroxyl hydrogen
    'HZ':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Histidine hydrogen
    'HA':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Aromatic hydrogen (general)

    'H4':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Aromatic hydrogen (HIS C4, TRP N1-H)
    'H5':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Aromatic hydrogen (HIS C2)
    # Carbon atoms
    'CT':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon
    'C':   {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Carbonyl carbon
    'CA':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon
    'CB':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (general)
    'CC':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (general)
    'CD':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Alkene carbon
    'CE':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Alkene carbon
    'CF':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Alkyne carbon
    'CG':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (general)
    'CH':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (general)
    'CR':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (5-membered ring)
    'CW':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (5-membered ring)
    'CV':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (5-membered ring)
    'CZ':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Phenyl carbon
    'C*':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (sp2)
    'CN':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (TRP indole C2)

    # Oxygen atoms
    'O':   {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Carbonyl oxygen
    'O2':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Carboxyl oxygen
    'OH':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Hydroxyl oxygen
    'OS':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Ether oxygen
    'OW':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Water oxygen

    # Nitrogen atoms
    'N':   {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Amide nitrogen
    'N3':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Protonated amine
    'NT':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Terminal amine
    'N2':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Amine nitrogen
    'NA':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Pyrrole nitrogen
    'NB':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Pyridine nitrogen
    'NC':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Cyano nitrogen
    'ND':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Pyrazine nitrogen
    'NE':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Imidazole nitrogen
    'NF':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Pyrrole nitrogen (general)
    'NG':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Indole nitrogen
    'NH':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Pyridine nitrogen (general)
    'NI':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Isoquinoline nitrogen
    'NL':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Aliphatic nitrogen
    'NM':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Amide nitrogen (general)
    'NP':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Phosphate nitrogen
    'NQ':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Guanidinium nitrogen

    # Sulfur atoms
    'S':   {'element': 'S',  'radius': 0.180, 'screen': 0.85, 'mass': 32.06},  # Thioether sulfur
    'SH':  {'element': 'S',  'radius': 0.180, 'screen': 0.85, 'mass': 32.06},  # Thiol sulfur
    'S*':  {'element': 'S',  'radius': 0.180, 'screen': 0.85, 'mass': 32.06},  # Sulfur (general)
    'SM':  {'element': 'S',  'radius': 0.180, 'screen': 0.85, 'mass': 32.06},  # Sulfhydryl sulfur

    # Halogen atoms
    'F':   {'element': 'F',  'radius': 0.150, 'screen': 0.85, 'mass': 19.00},  # Fluorine
    'CL':  {'element': 'Cl', 'radius': 0.180, 'screen': 0.85, 'mass': 35.45},  # Chlorine
    'BR':  {'element': 'Br', 'radius': 0.200, 'screen': 0.85, 'mass': 79.90},  # Bromine
    'I':   {'element': 'I',  'radius': 0.220, 'screen': 0.85, 'mass': 126.90}, # Iodine

    # Phosphorus
    'P':   {'element': 'P',  'radius': 0.185, 'screen': 0.85, 'mass': 30.97},  # Phosphate phosphorus

    # Unknown/default (fallback values)
    'X':   {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # Wildcard type

    # =========================================================================
    # Additional atom types for compatibility with other force fields
    # Added based on missing_gbsa_type_mappings.txt analysis
    # =========================================================================

    # a99SBdisp force field types
    'C1':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (a99SBdisp)
    'C3':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (a99SBdisp)
    'C4':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (a99SBdisp)
    'C5':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (a99SBdisp)
    'C6':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (a99SBdisp)
    'C7':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (a99SBdisp)
    'C8':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (a99SBdisp)
    'C9':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (a99SBdisp)
    'HB':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Backbone amide hydrogen (a99SBdisp)
    'O3':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Hydroxyl oxygen (a99SBdisp)
    'OB':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Carbonyl oxygen (a99SBdisp)

    # amber03wsc force field types
    'CAx': {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (amber03wsc)
    'CTx': {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (amber03wsc)
    'Cx':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Carbonyl carbon (amber03wsc)
    'H1x': {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Aliphatic hydrogen (amber03wsc)
    'HCx': {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Alkyl hydrogen (amber03wsc)
    'HPx': {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Aromatic hydrogen (amber03wsc)
    'Hx':  {'element': 'H',  'radius': 0.100, 'screen': 0.85, 'mass': 1.008},  # Hydrogen (amber03wsc)
    'N2x': {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Amine nitrogen (amber03wsc)
    'N3x': {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Protonated amine (amber03wsc)
    'Nx':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Amide nitrogen (amber03wsc)
    'O2x': {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Carboxyl oxygen (amber03wsc)
    'Ox':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Carbonyl oxygen (amber03wsc)

    # amber14sb_parmbsc1 force field types
    '2C':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (amber14sb_parmbsc1)
    '3C':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (amber14sb_parmbsc1)
    'CO':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # Carbonyl carbon (amber14sb_parmbsc1)
    'CX':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (amber14sb_parmbsc1)

    # des-amber force field types (single-letter codes for residues)
    'AA':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (ALA in des-amber)
    'DD':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (ASP in des-amber)
    'EE':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (GLU in des-amber)
    'FF':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (PHE in des-amber)
    'GG':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (GLY in des-amber)
    'HE':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (HIS in des-amber)
    'II':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (ILE in des-amber)
    'KK':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (LYS in des-amber)
    'LL':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (LEU in des-amber)
    'MM':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (MET in des-amber)
    'NN':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (ASN in des-amber)
    'PP':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (PRO in des-amber)
    'QQ':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (GLN in des-amber)
    'RR':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # sp2 carbon (ARG in des-amber)
    'SS':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (SER in des-amber)
    'TT':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (THR in des-amber)
    'VV':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # sp3 carbon (VAL in des-amber)
    'WW':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (TRP in des-amber)
    'YY':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Aromatic carbon (TYR in des-amber)

    # des-amber special atom types
    'C&':  {'element': 'C',  'radius': 0.175, 'screen': 0.75, 'mass': 12.01},  # Special carbon (des-amber)
    'CT_CT': {'element': 'C', 'radius': 0.170, 'screen': 0.72, 'mass': 12.01}, # sp3 carbon (des-amber)
    'H1_H1B': {'element': 'H', 'radius': 0.100, 'screen': 0.85, 'mass': 1.008}, # Aliphatic hydrogen (des-amber)
    'H_HN': {'element': 'H', 'radius': 0.100, 'screen': 0.85, 'mass': 1.008},   # Amide hydrogen (des-amber)
    'N_N':  {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Amide nitrogen (des-amber)
    'N_N2': {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.01},  # Amine nitrogen (des-amber)
    'OHS':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Hydroxyl oxygen (des-amber)
    'O_O':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 16.00},  # Carbonyl oxygen (des-amber)
}

# CHARMM36 atom type mapping for GBSA implicit solvent calculations
CHARMM36_ATOM_MAPPING = {
    # Hydrogen atoms
    'H':   {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Backbone amide hydrogen
    'HA':  {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Alpha hydrogen
    'HB1': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Beta hydrogen (methyl)
    'HB2': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Beta hydrogen (methylene)
    'HB3': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Beta hydrogen
    'HC':  {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Alkyl hydrogen
    'HA1': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Alpha hydrogen (ILE)
    'HA2': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Alpha hydrogen
    'HA3': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Methyl hydrogen
    'HP':  {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Aromatic hydrogen
    'HR1': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # HIS ring hydrogen
    'HR2': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # HIS ring hydrogen
    'HR3': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # HIS ring hydrogen
    'HS':  {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Thiol hydrogen
    'HZ':  {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # HIS ring hydrogen
    'HZ1': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Lysine hydrogen
    'HZ2': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Lysine hydrogen
    'HZ3': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Lysine hydrogen
    'HD1': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # HIS ring hydrogen
    'HD2': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # HIS ring hydrogen
    'HE1': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # HIS/PHE/TRP ring hydrogen
    'HE2': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # HIS/PHE ring hydrogen
    'HE3': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # MET methyl hydrogen
    'HG1': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Sidechain hydrogen
    'HG2': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # Sidechain hydrogen
    'HG3': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # MET methyl hydrogen
    'HH':  {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # TYR hydroxyl hydrogen
    'HH1': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # ARG hydrogen
    'HH2': {'element': 'H',  'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # ARG hydrogen
    'HG21': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},  # THR/VAL/ILE methyl hydrogen
    'HG22': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},
    'HG23': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},
    'HD11': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},  # LEU methyl hydrogen
    'HD12': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},
    'HD13': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},
    'HD21': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},  # ASN/GLN hydrogen
    'HD22': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},
    'HE1': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # TRP/PHE hydrogen
    'HE2': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},
    'HE3': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # MET methyl
    'HZ2': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # TRP hydrogen
    'HZ3': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # TRP hydrogen
    'HH2': {'element': 'H', 'radius': 0.105, 'screen': 0.85, 'mass': 1.008},   # TRP hydrogen

    # Carbon atoms
    'CT1': {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # Alpha carbon
    'CT2': {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # Beta carbon
    'CT3': {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # Methyl carbon
    'C':   {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Carbonyl carbon
    'CA':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Aromatic carbon
    'CB':  {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # Beta carbon (aromatic)
    'CC':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Carbonyl carbon (carboxyl)
    'CD':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Carbonyl carbon (carboxyl)
    'CE':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Alkene carbon
    'CF':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Aromatic carbon
    'CG':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Aromatic carbon
    'CH':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Aromatic carbon
    'CP1': {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # PRO alpha carbon
    'CP2': {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # PRO beta carbon
    'CP3': {'element': 'C',  'radius': 0.170, 'screen': 0.72, 'mass': 12.01},  # PRO delta carbon
    'CY':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # TRP carbon
    'CAI': {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # TRP indole carbon
    'CPT': {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # TRP carbon
    'CPH': {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # HIS carbon
    'CPH1': {'element': 'C', 'radius': 0.175, 'screen': 0.72, 'mass': 12.01}, # HIS carbon
    'CPH2': {'element': 'C', 'radius': 0.175, 'screen': 0.72, 'mass': 12.01}, # HIS carbon
    'CS':  {'element': 'C',  'radius': 0.175, 'screen': 0.72, 'mass': 12.01},  # Sulfonyl carbon
    'CT2A': {'element': 'C', 'radius': 0.170, 'screen': 0.72, 'mass': 12.01}, # ASP/GLU beta carbon
    'CT3A': {'element': 'C', 'radius': 0.170, 'screen': 0.72, 'mass': 12.01}, # Gamma carbon

    # Oxygen atoms
    'O':   {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 15.999},  # Carbonyl oxygen
    'O2':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 15.999},  # Carboxyl oxygen
    'OH1': {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 15.999},  # Hydroxyl oxygen
    'OC':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 15.999},  # Carboxyl oxygen
    'ON':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 15.999},  # Nitro oxygen
    'OS':  {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 15.999},  # Ether oxygen
    'OT1': {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 15.999},  # Terminal oxygen
    'OT2': {'element': 'O',  'radius': 0.150, 'screen': 0.85, 'mass': 15.999},  # Terminal oxygen

    # Nitrogen atoms
    'N':   {'element': 'N',  'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # Amide nitrogen
    'NH1': {'element': 'N', 'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # Amide nitrogen
    'NH2': {'element': 'N', 'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # Amide nitrogen (GLN/ASN)
    'NH3': {'element': 'N', 'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # Terminal amine
    'NC2': {'element': 'N', 'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # Guanidinium nitrogen
    'NR1': {'element': 'N', 'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # HIS ring nitrogen
    'NR2': {'element': 'N', 'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # HIS ring nitrogen
    'NR3': {'element': 'N', 'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # Protonated amine
    'NY':  {'element': 'N', 'radius': 0.155, 'screen': 0.85, 'mass': 14.007},  # TRP nitrogen

    # Sulfur atoms
    'S':   {'element': 'S',  'radius': 0.180, 'screen': 0.85, 'mass': 32.06},  # Thioether sulfur
    'SH1': {'element': 'S', 'radius': 0.180, 'screen': 0.85, 'mass': 32.06},  # Thiol sulfur
    'SM':  {'element': 'S', 'radius': 0.180, 'screen': 0.85, 'mass': 32.06},  # Sulfhydryl sulfur

    # Phosphorus
    'P':   {'element': 'P',  'radius': 0.185, 'screen': 0.85, 'mass': 30.974}, # Phosphate phosphorus

    # Halogen atoms
    'F':   {'element': 'F',  'radius': 0.150, 'screen': 0.85, 'mass': 18.998}, # Fluorine
    'CL':  {'element': 'Cl', 'radius': 0.180, 'screen': 0.85, 'mass': 35.45},  # Chlorine
    'BR':  {'element': 'Br', 'radius': 0.200, 'screen': 0.85, 'mass': 79.904}, # Bromine
    'I':   {'element': 'I',  'radius': 0.220, 'screen': 0.85, 'mass': 126.90}, # Iodine
}

novarcharre = re.compile(r'\W')

def _find_all_instances_in_string(string, substr):
    """ Find indices of all instances of substr in string """
    indices = []
    idx = string.find(substr, 0)
    while idx > -1:
        indices.append(idx)
        idx = string.find(substr, idx+1)
    return indices

def _replace_defines(line, defines):
    """ Replaces defined tokens in a given line """
    if not defines: return line
    for define in reversed(defines):
        value = defines[define]
        indices = _find_all_instances_in_string(line, define)
        if not indices: continue
        # Check to see if it's inside of quotes
        inside = ''
        idx = 0
        n_to_skip = 0
        new_line = []
        for i, char in enumerate(line):
            if n_to_skip:
                n_to_skip -= 1
                continue
            if char in ('\'"'):
                if not inside:
                    inside = char
                else:
                    if inside == char:
                        inside = ''
            if idx < len(indices) and i == indices[idx]:
                if inside:
                    new_line.append(char)
                    idx += 1
                    continue
                if i == 0 or novarcharre.match(line[i-1]):
                    endidx = indices[idx] + len(define)
                    if endidx >= len(line) or novarcharre.match(line[endidx]):
                        new_line.extend(list(value))
                        n_to_skip = len(define) - 1
                        idx += 1
                        continue
                idx += 1
            new_line.append(char)
        line = ''.join(new_line)

    return line

class GromacsTopFileWithSoftcore(object):
    """GromacsTopFile parses a Gromacs top file and constructs a Topology and (optionally) an OpenMM System from it.

    Extended with three-stage optimization support:
    - NONBONDED_STANDARD: Standard mode matching native gromacstopfile.py behavior
    - NONBONDED_GAUSSIAN: Gaussian repulsion for geometric untangling
    - NONBONDED_SOFTCORE: Softcore potential for free energy calculations (placeholder)
    """

    class _MoleculeType(object):
        """Inner class to store information about a molecule type."""
        def __init__(self, name, nrexcl):
            self.name = name
            self.nrexcl = nrexcl
            self.atoms = []
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.exclusions = []
            self.pairs = []
            self.constraints = []
            self.cmaps = []
            self.vsites2 = []
            self.vsites3 = []
            self.has_virtual_sites = False
            self.has_nbfix_terms = False

        def findExclusionsFromBonds(self, genpairs):
            """Find exclusions between atoms separated by up to nrexcl bonds if genpairs is false,
               or up to 2 bonds if genpairs is true.
            """
            bondedTo = [set() for i in range(len(self.atoms))]
            for fields in self.bonds:
                i = int(fields[0])-1
                j = int(fields[1])-1
                bondedTo[i].add(j)
                bondedTo[j].add(i)

            # Identify all neighbors of each atom with each separation.

            bondedWithSeparation = [bondedTo]
            maxBonds = self.nrexcl
            if genpairs:
                maxBonds = min(maxBonds, 2)
            for i in range(maxBonds-1):
                lastBonds = bondedWithSeparation[-1]
                newBonds = deepcopy(lastBonds)
                for atom in range(len(self.atoms)):
                    for a1 in lastBonds[atom]:
                        for a2 in bondedTo[a1]:
                            newBonds[atom].add(a2)
                bondedWithSeparation.append(newBonds)

            # Build the list of pairs.

            pairs = []
            for atom in range(len(self.atoms)):
                for otherAtom in bondedWithSeparation[-1][atom]:
                    if otherAtom > atom:
                        pairs.append((atom, otherAtom))
            return pairs

    def _processFile(self, file):
        append = ''
        for line in open(file):
            if line.strip().endswith('\\'):
                append = '%s %s' % (append, line[:line.rfind('\\')])
            else:
                self._processLine(append+' '+line, file)
                append = ''

    def _processLine(self, line, file):
        """Process one line from a file."""
        if ';' in line:
            line = line[:line.index(';')]
        stripped = line.strip()
        ignore = not all(self._ifStack)
        if stripped.startswith('*') or len(stripped) == 0:
            # A comment or empty line.
            return

        elif stripped.startswith('[') and not ignore:
            # The start of a category.
            if not stripped.endswith(']'):
                raise ValueError('Illegal line in .top file: '+line)
            self._currentCategory = stripped[1:-1].strip()

        elif stripped.startswith('#'):
            # A preprocessor command.
            fields = stripped.split()
            command = fields[0]
            if len(self._ifStack) != len(self._elseStack):
                raise RuntimeError('#if/#else stack out of sync')

            if command == '#include' and not ignore:
                # Locate the file to include
                name = stripped[len(command):].strip(' \t"<>')
                searchDirs = self._includeDirs+(os.path.dirname(file),)
                for dir in searchDirs:
                    file = os.path.join(dir, name)
                    if os.path.isfile(file):
                        # We found the file, so process it.
                        self._processFile(file)
                        break
                else:
                    raise ValueError('Could not locate #include file: '+name)
            elif command == '#define' and not ignore:
                # Add a value to our list of defines.
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                name = fields[1]
                valueStart = stripped.find(name, len(command))+len(name)+1
                value = line[valueStart:].strip()
                value = value or '1' # Default define is 1
                self._defines[name] = value
            elif command == '#ifdef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                name = fields[1]
                self._ifStack.append(name in self._defines)
                self._elseStack.append(False)
            elif command == '#undef':
                # Un-define a variable
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                if fields[1] in self._defines:
                    self._defines.pop(fields[1])
            elif command == '#ifndef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                name = fields[1]
                self._ifStack.append(name not in self._defines)
                self._elseStack.append(False)
            elif command == '#endif':
                # Pop an entry off the if stack.
                if len(self._ifStack) == 0:
                    raise ValueError('Unexpected line in .top file: '+line)
                del(self._ifStack[-1])
                del(self._elseStack[-1])
            elif command == '#else':
                # Reverse the last entry on the if stack
                if len(self._ifStack) == 0:
                    raise ValueError('Unexpected line in .top file: '+line)
                if self._elseStack[-1]:
                    raise ValueError('Unexpected line in .top file: '
                                     '#else has already been used ' + line)
                self._ifStack[-1] = (not self._ifStack[-1])
                self._elseStack[-1] = True

        elif not ignore:
            # Gromacs occasionally uses #define's to introduce specific
            # parameters for individual terms (for instance, this is how
            # ff99SB-ILDN is implemented). So make sure we do the appropriate
            # pre-processor replacements necessary
            line = _replace_defines(line, self._defines)
            # A line of data for the current category
            if self._currentCategory is None:
                raise ValueError('Unexpected line in .top file: '+line)
            if self._currentCategory == 'defaults':
                self._processDefaults(line)
            elif self._currentCategory == 'moleculetype':
                self._processMoleculeType(line)
            elif self._currentCategory == 'molecules':
                self._processMolecule(line)
            elif self._currentCategory == 'atoms':
                self._processAtom(line)
            elif self._currentCategory == 'bonds':
                self._processBond(line)
            elif self._currentCategory == 'angles':
                self._processAngle(line)
            elif self._currentCategory == 'dihedrals':
                self._processDihedral(line)
            elif self._currentCategory == 'exclusions':
                self._processExclusion(line)
            elif self._currentCategory == 'pairs':
                self._processPair(line)
            elif self._currentCategory == 'constraints':
                self._processConstraint(line)
            elif self._currentCategory == 'settles':
                self._processSettles(line)
            elif self._currentCategory == 'cmap':
                self._processCmap(line)
            elif self._currentCategory == 'atomtypes':
                self._processAtomType(line)
            elif self._currentCategory == 'bondtypes':
                self._processBondType(line)
            elif self._currentCategory == 'angletypes':
                self._processAngleType(line)
            elif self._currentCategory == 'dihedraltypes':
                self._processDihedralType(line)
            elif self._currentCategory == 'pairtypes':
                self._processPairType(line)
            elif self._currentCategory == 'cmaptypes':
                self._processCmapType(line)
            elif self._currentCategory == 'nonbond_params':
                self._processNonbondType(line)
            elif self._currentCategory == 'virtual_sites2' or self._currentCategory == 'dummies2':
                self._processVirtualSites2(line)
            elif self._currentCategory == 'virtual_sites3' or self._currentCategory == 'dummies3':
                self._processVirtualSites3(line)
            elif self._currentCategory.startswith('virtual_sites') or self._currentCategory.startswith('dummies'):
                if self._currentMoleculeType is None:
                    raise ValueError('Found %s before [ moleculetype ]' %
                                     self._currentCategory)
                self._currentMoleculeType.has_virtual_sites = True

    def _processDefaults(self, line):
        """Process the [ defaults ] line."""
        fields = line.split()
        if len(fields) < 5:
            # fudgeLJ and fudgeQQ not specified, assumed 1.0 by default
            if len(fields) == 3:
                fields.append(1.0)
                fields.append(1.0)
            else:
                raise ValueError('Too few fields in [ defaults ] line: '+line)
        if fields[0] != '1':
            raise ValueError('Unsupported nonbonded type: '+fields[0])
        if not fields[1] in ('1', '2', '3'):
            raise ValueError('Unsupported combination rule: '+fields[1])
        if fields[2].lower() == 'no':
            self._genpairs = False
        self._defaults = fields

    def _processMoleculeType(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            raise ValueError('Too few fields in [ moleculetypes ] line: '+line)
        type = gromacstopfile.GromacsTopFile._MoleculeType(fields[0], int(fields[1]))
        self._moleculeTypes[fields[0]] = type
        self._currentMoleculeType = type

    def _processMolecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            raise ValueError('Too few fields in [ molecules ] line: '+line)
        self._molecules.append((fields[0], int(fields[1])))

    def _processAtom(self, line):
        """Process a line in the [ atoms ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ atoms ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ atoms ] line: '+line)
        self._currentMoleculeType.atoms.append(fields)

    def _processBond(self, line):
        """Process a line in the [ bonds ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ bonds ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 3:
            raise ValueError('Too few fields in [ bonds ] line: '+line)
        if fields[2] not in ('1', '2'):
                raise ValueError('Unsupported function type in [ bonds ] line: '+line)
        self._currentMoleculeType.bonds.append(fields)

    def _processAngle(self, line):
        """Process a line in the [ angles ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ angles ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 4:
            raise ValueError('Too few fields in [ angles ] line: '+line)
        if fields[3] not in ('1', '2', '5'):
            raise ValueError('Unsupported function type in [ angles ] line: '+line)
        self._currentMoleculeType.angles.append(fields)

    def _processDihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ dihedrals ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ dihedrals ] line: '+line)
        if fields[4] not in ('1', '2', '3', '4', '5', '9'):
            raise ValueError('Unsupported function type in [ dihedrals ] line: '+line)
        self._currentMoleculeType.dihedrals.append(fields)

    def _processExclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ exclusions ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 2:
            raise ValueError('Too few fields in [ exclusions ] line: '+line)
        self._currentMoleculeType.exclusions.append(fields)

    def _processPair(self, line):
        """Process a line in the [ pairs ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ pairs ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 3:
            raise ValueError('Too few fields in [ pairs ] line: '+line)
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ pairs ] line: '+line)
        self._currentMoleculeType.pairs.append(fields)

    def _processConstraint(self, line):
        """Process a line in the [ constraints ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ constraints ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 4:
            raise ValueError('Too few fields in [ constraints ] line: '+line)
        self._currentMoleculeType.constraints.append(fields)

    def _processSettles(self, line):
        """Process a line in the [ settles ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ settles ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 4:
            raise ValueError('Too few fields in [ settles ] line: '+line)
        atom = int(fields[0])
        self._currentMoleculeType.constraints.append([str(atom), str(atom+1), fields[1], fields[2]])
        self._currentMoleculeType.constraints.append([str(atom), str(atom+2), fields[1], fields[2]])
        self._currentMoleculeType.constraints.append([str(atom+1), str(atom+2), fields[1], fields[3]])

    def _processCmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ cmap ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ cmap ] line: '+line)
        self._currentMoleculeType.cmaps.append(fields)

    def _processAtomType(self, line):
        """Process a line in the [ atomtypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ atomtypes ] line: '+line)
        if len(fields[3]) == 1:
            # Bonded type and atomic number are both missing.
            fields.insert(1, None)
            fields.insert(1, None)
        elif len(fields[4]) == 1 and fields[4].isalpha():
            if fields[1][0].isalpha():
                # Atomic number is missing.
                fields.insert(2, None)
            else:
                # Bonded type is missing.
                fields.insert(1, None)
        self._atomTypes[fields[0]] = fields

    def _processBondType(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ bondtypes ] line: '+line)
        if fields[2] not in ('1', '2'):
            raise ValueError('Unsupported function type in [ bondtypes ] line: '+line)
        self._bondTypes[tuple(fields[:3])] = fields

    def _processAngleType(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ angletypes ] line: '+line)
        if fields[3] not in ('1', '2', '5'):
            raise ValueError('Unsupported function type in [ angletypes ] line: '+line)
        self._angleTypes[tuple(fields[:3])] = fields

    def _processDihedralType(self, line):
        """Process a line in the [ dihedraltypes ] category."""
        fields = line.split()
        if len(fields[2]) == 1 and fields[2].isdigit():
            # The third field contains the function type, meaning only two atom types are specified.
            # Interpret them as the two inner ones.
            fields = ['X', fields[0], fields[1], 'X']+fields[2:]
        if len(fields) < 7:
            raise ValueError('Too few fields in [ dihedraltypes ] line: '+line)
        if fields[4] not in ('1', '2', '3', '4', '5', '9'):
            raise ValueError('Unsupported function type in [ dihedraltypes ] line: '+line)
        key = tuple(fields[:5])
        if fields[4] == '9' and key in self._dihedralTypes:
            # There are multiple dihedrals defined for these atom types.
            self._dihedralTypes[key].append(fields)
        else:
            self._dihedralTypes[key] = [fields]

    def _processPairType(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ pairtypes] line: '+line)
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ pairtypes ] line: '+line)
        self._pairTypes[tuple(fields[:2])] = fields

    def _processCmapType(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8+int(fields[6])*int(fields[7]):
            raise ValueError('Too few fields in [ cmaptypes ] line: '+line)
        if fields[5] != '1':
            raise ValueError('Unsupported function type in [ cmaptypes ] line: '+line)
        self._cmapTypes[tuple(fields[:5])] = fields

    def _processNonbondType(self, line):
        """Process a line in the [ nonbond_params ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ nonbond_params ] line: '+line)
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ nonbond_params ] line: '+line)
        self._nonbondTypes[tuple(sorted(fields[:2]))] = fields

    def _processVirtualSites2(self, line):
        """Process a line in the [ virtual_sites2 ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ virtual_sites2 ] line: ' + line)
        if fields[3] != '1':
            raise ValueError('Unsupported function type in [ virtual_sites2 ] line: '+line)
        self._currentMoleculeType.vsites2.append(fields[:5])

    def _processVirtualSites3(self, line):
        """Process a line in the [ virtual_sites3 ] category."""
        fields = line.split()
        if len(fields) < 7:
            raise ValueError('Too few fields in [ virtual_sites3 ] line: ' + line)
        if fields[4] not in ('1', '4'):
            raise ValueError('Unsupported function type in [ virtual_sites3 ] line: '+line)
        self._currentMoleculeType.vsites3.append(fields)

    def __init__(self, file, periodicBoxVectors=None, unitCellDimensions=None, includeDir=None, defines=None, forcefield_type='AMBER'):
        """Load a top file.

        Parameters
        ----------
        file : str
            the name of the file to load
        periodicBoxVectors : tuple of Vec3=None
            the vectors defining the periodic box
        unitCellDimensions : Vec3=None
            the dimensions of the crystallographic unit cell.  For
            non-rectangular unit cells, specify periodicBoxVectors instead.
        includeDir : string=None
            A directory in which to look for other files included from the
            top file. If not specified, we will attempt to locate a gromacs
            installation on your system. When gromacs is installed in
            /usr/local, this will resolve to /usr/local/gromacs/share/gromacs/top
        defines : dict={}
            preprocessor definitions that should be predefined when parsing the file
        forcefield_type : str='AMBER'
            Force field type for GBSA calculations: 'AMBER' or 'CHARMM'
         """
        # Validate forcefield_type
        if forcefield_type.upper() not in ['AMBER', 'CHARMM']:
            raise ValueError(f"Unsupported forcefield_type: {forcefield_type}. "
                           f"Supported types: 'AMBER', 'CHARMM'")
        self._forcefield_type = forcefield_type.upper()
        if includeDir is None:
            includeDir = _defaultGromacsIncludeDir()
        self._includeDirs = (os.path.dirname(file), includeDir)
        # Most of the gromacs water itp files for different forcefields,
        # unless the preprocessor #define FLEXIBLE is given, don't define
        # bonds between the water hydrogen and oxygens, but only give the
        # constraint distances and exclusions.
        self._defines = OrderedDict()
        self._defines['FLEXIBLE'] = True
        self._genpairs = True
        if defines is not None:
            for define, value in defines.items():
                self._defines[define] = value

        # Parse the file.

        self._currentCategory = None
        self._ifStack = []
        self._elseStack = []
        self._moleculeTypes = {}
        self._molecules = []
        self._currentMoleculeType = None
        self._atomTypes = {}
        self._bondTypes= {}
        self._angleTypes = {}
        self._dihedralTypes = {}
        self._pairTypes = {}
        self._cmapTypes = {}
        self._nonbondTypes = {}
        self._processFile(file)

        # Create the Topology from it.

        top = Topology()
        ## The Topology read from the prmtop file
        self.topology = top
        if periodicBoxVectors is not None:
            if unitCellDimensions is not None:
                raise ValueError("specify either periodicBoxVectors or unitCellDimensions, but not both")
            top.setPeriodicBoxVectors(periodicBoxVectors)
        else:
            top.setUnitCellDimensions(unitCellDimensions)
        # DISABLED: PDBFile._loadNameReplacementTables()
        # We use original atom names from GROMACS topology to preserve naming consistency
        for moleculeName, moleculeCount in self._molecules:
            if moleculeName not in self._moleculeTypes:
                raise ValueError("Unknown molecule type: "+moleculeName)
            moleculeType = self._moleculeTypes[moleculeName]
            if moleculeCount > 0 and moleculeType.has_virtual_sites:
                raise ValueError('Virtual sites not yet supported by Gromacs parsers')

            # Create the specified number of molecules of this type.

            for i in range(moleculeCount):
                atoms = []
                lastResidue = None
                c = top.addChain()
                for index, fields in enumerate(moleculeType.atoms):
                    resNumber = fields[2]
                    if resNumber != lastResidue:
                        lastResidue = resNumber
                        resName = fields[3]
                        # DISABLED: Skip residue name replacement to preserve GROMACS naming
                        # if resName in PDBFile._residueNameReplacements:
                        #     resName = PDBFile._residueNameReplacements[resName]
                        r = top.addResidue(resName, c)
                        # DISABLED: Skip atom name replacement to preserve GROMACS naming
                        # if resName in PDBFile._atomNameReplacements:
                        #     atomReplacements = PDBFile._atomNameReplacements[resName]
                        # else:
                        #     atomReplacements = {}
                    atomName = fields[4]
                    # DISABLED: Skip atom name replacement to preserve GROMACS naming
                    # if atomName in atomReplacements:
                    #     atomName = atomReplacements[atomName]

                    # Try to determine the element.

                    atomicNumber = self._atomTypes[fields[1]][2]
                    if atomicNumber is None:
                        # Try to guess the element from the name.
                        upper = atomName.upper()
                        if upper.startswith('CL'):
                            element = elem.chlorine
                        elif upper.startswith('NA'):
                            element = elem.sodium
                        elif upper.startswith('MG'):
                            element = elem.magnesium
                        else:
                            try:
                                element = elem.get_by_symbol(atomName[0])
                            except KeyError:
                                element = None
                    elif atomicNumber == '0':
                        element = None
                    else:
                        element = elem.Element.getByAtomicNumber(int(atomicNumber))
                    atoms.append(top.addAtom(atomName, element, r))

                # Add bonds to the topology

                for fields in moleculeType.bonds:
                    top.addBond(atoms[int(fields[0])-1], atoms[int(fields[1])-1])

    def createSystem(self, nonbondedMethod=ff.CutoffPeriodic, nonbondedCutoff=1.0*unit.nanometer, constraints=None,
                     rigidWater=True, ewaldErrorTolerance=0.0005, removeCMMotion=True, hydrogenMass=None, switchDistance=None,
                     nonbonded_type=NONBONDED_STANDARD, add_implicit_solvent=False, gb_model='OBC2', salt_conc=0.0,
                     gaussian_width=0.085*unit.nanometer,
                     soft_lambda=0.85, soft_alpha_lj=0.85, soft_alpha_coul=0.3, soft_sigma_coul=1.0):
        """Construct an OpenMM System representing the topology described by this
        top file.

        Parameters
        ----------
        nonbondedMethod : object=CutoffPeriodic
            The method to use for nonbonded interactions. Hardcoded to CutoffPeriodic.
            Allowed values are CutoffPeriodic (only supported value).
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with
            constraints. Allowed values are None, HBonds, AllBonds, or HAngles.
            Regardless of this value, constraints that are explicitly specified
            in the top file will always be included.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        ewaldErrorTolerance : float=0.0005
            The error tolerance for CutoffPeriodic (not used for Ewald/PME since those methods are not supported)
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep their
            total mass the same.  If rigidWater is used to make water molecules
            rigid, then water hydrogens are not altered.
        switchDistance : float=None
            The distance at which the potential energy switching function is turned on for
            Lennard-Jones interactions. If this is None, no switching function will be used.
        nonbonded_type : str="standard"
            Type of nonbonded force to use for three-stage optimization:
            - "gaussian": Gaussian repulsion (Stage 1 - initial untangling, no chemistry)
            - "softcore": Gapsys linearized soft-core (Stage 2 - stable EM, eliminates singularities)
            - "standard": Standard potential (Stage 3 - final sampling, identical to native gromacstopfile.py)
        add_implicit_solvent : bool=False
            Whether to add implicit solvent (only valid when nonbonded_type="standard")
        gb_model : str="OBC2"
            Implicit solvent model: "OBC1", "OBC2", "GBn", "GBn2" [NOT YET IMPLEMENTED]
        salt_conc : float=0.0
            Salt concentration (M) for Debye-Hückel screening [NOT YET IMPLEMENTED]
        gaussian_width : Quantity=0.085*nanometer
            Width parameter (ga_w) for Gaussian repulsion force in nm.
            Controls the "stiffness" of the repulsion.
        soft_lambda : float=0.85
            Softness parameter for softcore potential (Gapsys approach).
            Values:
            - 1.0: Hard core (no softening, like standard LJ)
            - 0.0: Maximum softness (atoms can overlap completely)
            - 0.85: Recommended for stable EM (default)
        soft_alpha_lj : float=0.85
            LJ soft-core alpha parameter (Gapsys paper recommends 0.85).
            Controls the switching distance for LJ: r_sw_lj = alpha * 1.244 * sigma * (1-lambda)^1/6
        soft_alpha_coul : float=0.3
            Coulomb soft-core alpha parameter (Gapsys paper recommends 0.3).
            Controls the switching distance for Coulomb: r_sw_q = alpha * (1 + sigma_coul*|q1q2|) * (1-lambda)^1/6
        soft_sigma_coul : float=1.0
            Coulomb soft-core sigma parameter (Gapsys paper recommends 1.0).
            Scales the charge contribution to switching distance.

        Returns
        -------
        System
             the newly created System
        """

        # Build a list of atom types for NBFIX

        atom_types = []
        for moleculeName, moleculeCount in self._molecules:
            moleculeType = self._moleculeTypes[moleculeName]
            for _ in range(moleculeCount):
                for atom in moleculeType.atoms:
                    atom_types.append(atom[1])
        has_nbfix_terms = any([pair in self._nonbondTypes for pair in combinations_with_replacement(sorted(set(atom_types)), 2)])

        # Create the System.

        sys = mm.System()
        boxVectors = self.topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(*boxVectors)
        elif nonbondedMethod == ff.CutoffPeriodic and boxVectors is None:
            raise ValueError('Illegal nonbonded method for a non-periodic system')

        # Hardcoded to CutoffPeriodic for all modes
        use_three_stage = nonbonded_type != NONBONDED_STANDARD or add_implicit_solvent

        if not use_three_stage:
            # STANDARD MODE: Use identical logic to native gromacstopfile.py
            # This is a direct copy from original implementation
            nb = mm.NonbondedForce()
            sys.addForce(nb)
            lj = None
            if has_nbfix_terms:
                lj = mm.CustomNonbondedForce('(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2)')
                lj.addPerParticleParameter('type')
                sys.addForce(lj)
            elif self._defaults[1] in ('1', '3'):
                lj = mm.CustomNonbondedForce('A1*A2/r^12-C1*C2/r^6')
                lj.addPerParticleParameter('C')
                lj.addPerParticleParameter('A')
                sys.addForce(lj)
            # For combination rule 2: lj remains None, use NonbondedForce only
        else:
            # THREE-STAGE OPTIMIZATION MODE
            if nonbonded_type == NONBONDED_GAUSSIAN:
                # Stage 1: Gaussian repulsion for geometric untangling
                self.nb_force = _create_gaussian_nb_force(nonbondedCutoff, gaussian_width.value_in_unit(unit.nanometer))
                sys.addForce(self.nb_force)
                nb = mm.NonbondedForce()  # For completeness
                sys.addForce(nb)
                lj = None  # Not used in gaussian mode
            elif nonbonded_type == NONBONDED_SOFTCORE:
                # Stage 2: Gapsys linearized soft-core potential
                # Based on Gapsys et al., JCTC 2015, 11, 5920-5930
                # Eliminates singularities for stable Energy Minimization
                
                # Create the softcore force - handles BOTH LJ and Coulomb
                self.nb_force = _create_gapsys_linearized_nb_force(
                    nonbondedCutoff, 
                    current_lambda=soft_lambda,
                    alpha_lj=soft_alpha_lj, 
                    alpha_coul=soft_alpha_coul, 
                    sigma_coul=soft_sigma_coul,
                    switchDistance=switchDistance,
                    has_nbfix_terms=has_nbfix_terms
                )
                sys.addForce(self.nb_force)
                
                # No NonbondedForce needed - softcore CustomNonbondedForce handles everything
                nb = None
                lj = None

                # Add GBSA implicit solvent if requested (skip for CHARMM)
                if add_implicit_solvent and self._forcefield_type.upper() != 'CHARMM':
                    _add_gbsa_solvent(sys, self, gb_model, salt_conc, nonbondedCutoff, self._forcefield_type)
            elif nonbonded_type == NONBONDED_STANDARD_CUSTOM:
                # Standard LJ+Coulomb using CustomNonbondedForce (for debugging/comparison)
                # This implements the exact same formula as standard NonbondedForce
                # but using CustomNonbondedForce architecture to match softcore mode
                
                # Create CustomNonbondedForce with standard LJ+Coulomb formula
                self.nb_force = _create_standard_custom_nb_force(
                    nonbondedCutoff,
                    switchDistance=switchDistance
                )
                sys.addForce(self.nb_force)
                
                # No NonbondedForce needed
                nb = None
                lj = None

                # Add GBSA implicit solvent if requested (skip for CHARMM)
                if add_implicit_solvent and self._forcefield_type.upper() != 'CHARMM':
                    _add_gbsa_solvent(sys, self, gb_model, salt_conc, nonbondedCutoff, self._forcefield_type)
            elif nonbonded_type == NONBONDED_STANDARD and add_implicit_solvent:
                # Stage 3: Standard potential with implicit solvent
                # Use NonbondedForce for standard electrostatics
                # Use CustomNonbondedForce for LJ (if needed)
                # Add GBSAGBn2Force for implicit solvent

                # Create NonbondedForce for electrostatics
                nb = mm.NonbondedForce()
                sys.addForce(nb)

                # Create CustomNonbondedForce for LJ if needed
                lj = None
                if has_nbfix_terms:
                    lj = mm.CustomNonbondedForce('(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2)')
                    lj.addPerParticleParameter('type')
                    sys.addForce(lj)
                elif self._defaults[1] in ('1', '3'):
                    lj = mm.CustomNonbondedForce('A1*A2/r^12-C1*C2/r^6')
                    lj.addPerParticleParameter('C')
                    lj.addPerParticleParameter('A')
                    sys.addForce(lj)

                # Add GBn2 implicit solvent (skip for CHARMM)
                if self._forcefield_type.upper() != 'CHARMM':
                    _add_gbsa_solvent(sys, self, gb_model, salt_conc, nonbondedCutoff, self._forcefield_type)
            else:
                # Fallback (shouldn't happen)
                raise ValueError(f"Invalid combination: nonbonded_type='{nonbonded_type}', add_implicit_solvent={add_implicit_solvent}")
        bonds = {}
        angles = {}
        periodic = None
        rb = None
        harmonicTorsion = None
        cmap = None
        mapIndices = {}
        bondIndices = []
        topologyAtoms = list(self.topology.atoms())
        exclusions = []
        pairs = []
        particle_params = []  # Store [q, sigma, epsilon] for softcore mode
        fudgeQQ = float(self._defaults[4])
        fudgeLJ = float(self._defaults[3])

        # Build a lookup table to let us process dihedrals more quickly.

        dihedralTypeTable = {}
        for key in self._dihedralTypes:
            if key[1] != 'X' and key[2] != 'X':
                if (key[1], key[2]) not in dihedralTypeTable:
                    dihedralTypeTable[(key[1], key[2])] = []
                dihedralTypeTable[(key[1], key[2])].append(key)
                if (key[2], key[1]) not in dihedralTypeTable:
                    dihedralTypeTable[(key[2], key[1])] = []
                dihedralTypeTable[(key[2], key[1])].append(key)
        wildcardDihedralTypes = []
        for key in self._dihedralTypes:
            if key[1] == 'X' or key[2] == 'X':
                wildcardDihedralTypes.append(key)
                for types in dihedralTypeTable.values():
                    types.append(key)

        if has_nbfix_terms:
            # Build a lookup table and angle/dihedral indices list to
            # let us handle exclusion manually.
            angleIndices = []
            torsionIndices = []
            atom_partners = defaultdict(lambda : defaultdict(set))
        else:
            atom_partners = None

        # atom_charges is always needed in three-stage mode for dihedral processing
        atom_charges = []

        # Loop over molecules and create the specified number of each type.

        for moleculeName, moleculeCount in self._molecules:
            moleculeType = self._moleculeTypes[moleculeName]
            exclusionsFromBonds = moleculeType.findExclusionsFromBonds(self._genpairs)
            for i in range(moleculeCount):

                # Record the types of all atoms.

                baseAtomIndex = sys.getNumParticles()
                atomTypes = [atom[1] for atom in moleculeType.atoms]
                try:
                    bondedTypes = [self._atomTypes[t][1] for t in atomTypes]
                except KeyError as e:
                    raise ValueError('Unknown atom type: ' + e.message)
                bondedTypes = [b if b is not None else a for a, b in zip(atomTypes, bondedTypes)]

                # Add atoms.

                for fields in moleculeType.atoms:
                    if len(fields) >= 8:
                        mass = float(fields[7])
                    else:
                        mass = float(self._atomTypes[fields[1]][3])
                    sys.addParticle(mass)

                # Add bonds.

                atomBonds = [{} for x in range(len(moleculeType.atoms))]
                for fields in moleculeType.bonds:
                    atoms = [int(x)-1 for x in fields[:2]]
                    types = tuple(bondedTypes[i] for i in atoms)
                    bondType = fields[2]
                    reversedTypes = types[::-1]+(bondType,)
                    types = types+(bondType,)
                    if len(fields) >= 5:
                        params = fields[3:5]
                    elif types in self._bondTypes:
                        params = self._bondTypes[types][3:5]
                    elif reversedTypes in self._bondTypes:
                        params = self._bondTypes[reversedTypes][3:5]
                    else:
                        raise ValueError('No parameters specified for bond: '+fields[0]+', '+fields[1])
                    # Decide whether to use a constraint or a bond.
                    useConstraint = False
                    if rigidWater and topologyAtoms[baseAtomIndex+atoms[0]].residue.name == 'HOH':
                        useConstraint = True
                    if constraints in (AllBonds, HAngles):
                        useConstraint = True
                    elif constraints is HBonds:
                        elements = [topologyAtoms[baseAtomIndex+i].element for i in atoms]
                        if elem.hydrogen in elements:
                            useConstraint = True
                    # Add the bond or constraint.
                    length = float(params[0])
                    if useConstraint:
                        sys.addConstraint(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], length)
                    elif bondType == '1':
                        if bondType not in bonds:
                            bonds[bondType] = mm.HarmonicBondForce()
                            sys.addForce(bonds[bondType])
                        bonds[bondType].addBond(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], length, float(params[1]))
                    elif bondType == '2':
                        if bondType not in bonds:
                            bonds[bondType] = mm.CustomBondForce('0.25*k*(r^2-r0^2)^2')
                            bonds[bondType].addPerBondParameter('r0')
                            bonds[bondType].addPerBondParameter('k')
                            bonds[bondType].setName('GROMOSBondForce')
                            sys.addForce(bonds[bondType])
                        bonds[bondType].addBond(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], (length, float(params[1])))
                    else:
                        raise ValueError('Internal error: bondType has unexpected value: '+bondType)
                    # Record information that will be needed for constraining angles.
                    atomBonds[atoms[0]][atoms[1]] = length
                    atomBonds[atoms[1]][atoms[0]] = length

                # Add angles.

                degToRad = math.pi/180
                for fields in moleculeType.angles:
                    atoms = [int(x)-1 for x in fields[:3]]
                    types = tuple(bondedTypes[i] for i in atoms)
                    angleType = fields[3]
                    if len(fields) >= 6:
                        params = fields[4:]
                    elif types in self._angleTypes:
                        params = self._angleTypes[types][4:]
                    elif types[::-1] in self._angleTypes:
                        params = self._angleTypes[types[::-1]][4:]
                    else:
                        raise ValueError('No parameters specified for angle: '+fields[0]+', '+fields[1]+', '+fields[2])
                    # Decide whether to use a constraint or a bond.
                    useConstraint = False
                    if rigidWater and topologyAtoms[baseAtomIndex+atoms[0]].residue.name == 'HOH':
                        useConstraint = True
                    if constraints is HAngles:
                        elements = [topologyAtoms[baseAtomIndex+i].element for i in atoms]
                        if elements[0] == elem.hydrogen and elements[2] == elem.hydrogen:
                            useConstraint = True
                        elif elements[1] == elem.oxygen and (elements[0] == elem.hydrogen or elements[2] == elem.hydrogen):
                            useConstraint = True
                    # Add the bond or constraint.
                    theta = float(params[0])*degToRad
                    if useConstraint:
                        # Compute the distance between atoms and add a constraint
                        if atoms[0] in atomBonds[atoms[1]] and atoms[2] in atomBonds[atoms[1]]:
                            l1 = atomBonds[atoms[1]][atoms[0]]
                            l2 = atomBonds[atoms[1]][atoms[2]]
                            length = math.sqrt(l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta))
                            sys.addConstraint(baseAtomIndex+atoms[0], baseAtomIndex+atoms[2], length)
                    else:
                        if angleType in ('1', '5'):
                            if angleType not in angles:
                                angles[angleType] = mm.HarmonicAngleForce()
                                sys.addForce(angles[angleType])
                            angles[angleType].addAngle(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], theta, float(params[1]))
                            if angleType == '5':
                                # This is a Urey-Bradley term, so also add the bond.
                                if '1' not in bonds:
                                    bonds['1'] = mm.HarmonicBondForce()
                                    sys.addForce(bonds['1'])
                                k = float(params[3])
                                if k != 0:
                                    bonds['1'].addBond(baseAtomIndex + atoms[0], baseAtomIndex + atoms[2], float(params[2]), k)
                        elif angleType == '2':
                            if angleType not in angles:
                                angles[angleType] = mm.CustomAngleForce('0.5*k*(cos(theta)-cos(theta0))^2')
                                angles[angleType].addPerAngleParameter('theta0')
                                angles[angleType].addPerAngleParameter('k')
                                angles[angleType].setName('GROMOSAngleForce')
                                sys.addForce(angles[angleType])
                            angles[angleType].addAngle(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], (theta, float(params[1])))
                        else:
                            raise ValueError('Internal error: angleType has unexpected value: '+angleType)

                # Add torsions.

                for fields in moleculeType.dihedrals:
                    atoms = [int(x)-1 for x in fields[:4]]
                    types = tuple(bondedTypes[i] for i in atoms)
                    dihedralType = fields[4]
                    reversedTypes = types[::-1]+(dihedralType,)
                    types = types+(dihedralType,)
                    if (dihedralType in ('1', '4', '5', '9') and len(fields) > 7) or (dihedralType == '3' and len(fields) > 10) or (dihedralType == '2' and len(fields) > 6):
                        paramsList = [fields]
                    else:
                        # Look for a matching dihedral type.
                        paramsList = None
                        if (types[1], types[2]) in dihedralTypeTable:
                            dihedralTypes = dihedralTypeTable[(types[1], types[2])]
                        else:
                            dihedralTypes = wildcardDihedralTypes
                        for key in dihedralTypes:
                            if all(a == b or a == 'X' for a, b in zip(key, types)) or all(a == b or a == 'X' for a, b in zip(key, reversedTypes)):
                                paramsList = self._dihedralTypes[key]
                                if 'X' not in key:
                                    break
                        if paramsList is None:
                            raise ValueError('No parameters specified for dihedral: '+fields[0]+', '+fields[1]+', '+fields[2]+', '+fields[3])
                    for params in paramsList:
                        if dihedralType in ('1', '4', '9'):
                            # Periodic torsion
                            k = float(params[6])
                            if k != 0:
                                if periodic is None:
                                    periodic = mm.PeriodicTorsionForce()
                                    sys.addForce(periodic)
                                periodic.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], int(float(params[7])), float(params[5])*degToRad, k)
                        elif dihedralType == '2':
                            # Harmonic torsion
                            k = float(params[6])
                            phi0 = float(params[5])
                            if k != 0:
                                if harmonicTorsion is None:
                                    harmonicTorsion = mm.CustomTorsionForce('0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g' % math.pi)
                                    harmonicTorsion.addPerTorsionParameter('theta0')
                                    harmonicTorsion.addPerTorsionParameter('k')
                                    harmonicTorsion.setName('HarmonicTorsionForce')
                                    sys.addForce(harmonicTorsion)
                                # map phi0 into correct space
                                phi0 = phi0 - 360 if phi0 > 180 else phi0
                                harmonicTorsion.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], (phi0*degToRad, k))
                        else:
                            # RB Torsion
                            c = [float(x) for x in params[5:11]]
                            if any(x != 0 for x in c):
                                if rb is None:
                                    rb = mm.RBTorsionForce()
                                    sys.addForce(rb)
                                if dihedralType == '5':
                                    # Convert Fourier coefficients to RB coefficients.
                                    c = [c[1]+0.5*(c[0]+c[2]), 0.5*(-c[0]+3*c[2]), -c[1]+4*c[3], -2*c[2], -4*c[3], 0]
                                rb.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], c[0], c[1], c[2], c[3], c[4], c[5])

                # Add CMAP terms.

                for fields in moleculeType.cmaps:
                    atoms = [int(x)-1 for x in fields[:5]]
                    types = tuple(bondedTypes[i] for i in atoms)
                    if len(fields) >= 8 and len(fields) >= 8+int(fields[6])*int(fields[7]):
                        params = fields
                    elif types in self._cmapTypes:
                        params = self._cmapTypes[types]
                    elif types[::-1] in self._cmapTypes:
                        params = self._cmapTypes[types[::-1]]
                    else:
                        raise ValueError('No parameters specified for cmap: '+fields[0]+', '+fields[1]+', '+fields[2]+', '+fields[3]+', '+fields[4])
                    if cmap is None:
                        cmap = mm.CMAPTorsionForce()
                        sys.addForce(cmap)
                    mapSize = int(params[6])
                    if mapSize != int(params[7]):
                        raise ValueError('Non-square CMAPs are not supported')
                    map = []
                    for i in range(mapSize):
                        for j in range(mapSize):
                            map.append(float(params[8+mapSize*((j+mapSize//2)%mapSize)+((i+mapSize//2)%mapSize)]))
                    map = tuple(map)
                    if map not in mapIndices:
                        mapIndices[map] = cmap.addMap(mapSize, map)
                    cmap.addTorsion(mapIndices[map], baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3],
                                 baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], baseAtomIndex+atoms[4])

                # Set nonbonded parameters for particles.

                for fields in moleculeType.atoms:
                    params = self._atomTypes[fields[1]]

                    if len(fields) > 6:
                        q = float(fields[6])
                    else:
                        q = float(params[4])

                    # SOFTCORE/STANDARD_CUSTOM mode: doesn't use nbfix or standard nb/lj
                    # Check this BEFORE has_nbfix_terms since nb and lj are None in this mode
                    if nonbonded_type in (NONBONDED_SOFTCORE, NONBONDED_STANDARD_CUSTOM) and use_three_stage:
                        # Softcore/Standard_Custom mode: add parameters to CustomNonbondedForce
                        # Store parameters for later use in pair/exception processing
                        
                        if has_nbfix_terms and nonbonded_type == NONBONDED_SOFTCORE:
                            # NBFIX mode with softcore: Use type index + charge
                            # Type index will be set later after building type map
                            self.nb_force.addParticle([0, q])  # [type, charge]
                            atom_charges.append(q)
                            # Store for later type assignment
                            particle_params.append([0, q])  # Placeholder for type
                        else:
                            # Non-NBFIX mode: Convert from combination rule format to sigma/epsilon
                            if self._defaults[1] == '1':
                                # Rule 1: C6, C12 -> sigma, epsilon
                                C6 = float(params[6])
                                C12 = float(params[7])
                                sigma = (C12 / C6) ** (1.0 / 6.0) if C6 > 0 else 0.0
                                epsilon = C6 ** 2 / (4.0 * C12) if C12 > 0 else 0.0
                            elif self._defaults[1] == '2':
                                # Rule 2: sigma, epsilon directly
                                sigma = float(params[6])
                                epsilon = float(params[7])
                            elif self._defaults[1] == '3':
                                # Rule 3: sigma, epsilon -> convert to C6, C12, then back
                                # Actually for rule 3: C6 = 4*epsilon*sigma^6, C12 = 4*epsilon*sigma^12
                                # We need sigma, epsilon for the softcore force
                                sigma = float(params[6])
                                epsilon = float(params[7])
                            self.nb_force.addParticle([q, sigma, epsilon])
                            particle_params.append([q, sigma, epsilon])
                    # GAUSSIAN mode: doesn't need LJ parameters or nbfix
                    elif nonbonded_type == NONBONDED_GAUSSIAN:
                        # Gaussian repulsion mode: add dummy particles to CustomNonbondedForce
                        # The gaussian force only cares about positions, not particle properties
                        self.nb_force.addParticle([0, 0])
                        # Also add to NonbondedForce for completeness (though not used)
                        nb.addParticle(0, 0, 0)
                        # But still need to track charges for dihedral processing
                        atom_charges.append(q)
                    elif has_nbfix_terms:
                        nb.addParticle(q, 1.0, 0.0)
                        atom_charges.append(q)
                        lj.addParticle([0])
                    elif self._defaults[1] == '1':
                        nb.addParticle(q, 1.0, 0.0)
                        lj.addParticle([math.sqrt(float(params[6])), math.sqrt(float(params[7]))])
                    elif self._defaults[1] == '2':
                        nb.addParticle(q, float(params[6]), float(params[7]))
                    elif self._defaults[1] == '3':
                        nb.addParticle(q, 1.0, 0.0)
                        sigma = float(params[6])
                        epsilon = float(params[7])
                        lj.addParticle([math.sqrt(4*epsilon*sigma**6), math.sqrt(4*epsilon*sigma**12)])

                for fields in moleculeType.bonds:
                    atoms = [int(x)-1 for x in fields[:2]]
                    bondIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1]))
                for fields in moleculeType.constraints:
                    if fields[2] == '1':
                        atoms = [int(x)-1 for x in fields[:2]]
                        bondIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1]))
                if has_nbfix_terms:
                    for fields in moleculeType.bonds:
                        atoms = [int(x)-1 for x in fields[:2]]
                        atom_partners[baseAtomIndex+atoms[0]]['bond'].add(baseAtomIndex+atoms[1])
                        atom_partners[baseAtomIndex+atoms[1]]['bond'].add(baseAtomIndex+atoms[0])
                    for fields in moleculeType.angles:
                        atoms = [int(x)-1 for x in fields[:3]]
                        angleIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2]))
                        for pair in combinations(atoms, 2):
                            atom_partners[baseAtomIndex+pair[0]]['angle'].add(baseAtomIndex+pair[1])
                            atom_partners[baseAtomIndex+pair[1]]['angle'].add(baseAtomIndex+pair[0])
                    for fields in moleculeType.dihedrals:
                        atoms = [int(x)-1 for x in fields[:4]]
                        torsionIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3]))
                        for pair in combinations(atoms, 2):
                            atom_partners[baseAtomIndex+pair[0]]['torsion'].add(baseAtomIndex+pair[1])
                            atom_partners[baseAtomIndex+pair[1]]['torsion'].add(baseAtomIndex+pair[0])

                # Record nonbonded exceptions.

                for fields in moleculeType.pairs:
                    atoms = [int(x)-1 for x in fields[:2]]
                    types = tuple(atomTypes[i] for i in atoms)
                    # Get particle parameters from nb or particle_params (softcore mode)
                    if nb is not None:
                        atom1params = nb.getParticleParameters(baseAtomIndex+atoms[0])
                        atom2params = nb.getParticleParameters(baseAtomIndex+atoms[1])
                        atom1params = [x.value_in_unit_system(unit.md_unit_system) for x in atom1params]
                        atom2params = [x.value_in_unit_system(unit.md_unit_system) for x in atom2params]
                    else:
                        # Softcore mode: use stored particle_params
                        atom1params = particle_params[baseAtomIndex+atoms[0]]
                        atom2params = particle_params[baseAtomIndex+atoms[1]]
                    if len(fields) >= 5:
                        params = [float(x) for x in fields[3:5]]
                    elif types in self._pairTypes:
                        params = [float(x) for x in self._pairTypes[types][3:5]]
                    elif types[::-1] in self._pairTypes:
                        params = [float(x) for x in self._pairTypes[types[::-1]][3:5]]
                    elif not self._genpairs:
                        raise ValueError('No pair parameters defined for atom '
                                         'types %s and gen-pairs is "no"' % types)
                    elif has_nbfix_terms:
                        continue
                    else:
                        # Generate the parameters based on the atom parameters.
                        if self._defaults[1] == '2':
                            params = [0.5*(atom1params[1]+atom2params[1]), fudgeLJ*math.sqrt(atom1params[2]*atom2params[2])]
                        elif nb is not None:
                            atom1lj = lj.getParticleParameters(baseAtomIndex+atoms[0])
                            atom2lj = lj.getParticleParameters(baseAtomIndex+atoms[1])
                            params = [fudgeLJ*atom1lj[0]*atom2lj[0], fudgeLJ*atom1lj[1]*atom2lj[1]]
                        else:
                            # Softcore mode: use particle_params directly
                            params = [fudgeLJ*math.sqrt(atom1params[2]*atom2params[2]), 0.0]
                    pairs.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], atom1params[0]*atom2params[0]*fudgeQQ, params[0], params[1]))
                for fields in moleculeType.exclusions:
                    atoms = [int(x)-1 for x in fields]
                    for atom in atoms[1:]:
                        exclusions.append((baseAtomIndex+atoms[0], baseAtomIndex+atom))
                for atoms in exclusionsFromBonds:
                    exclusions.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1]))

                # Record virtual sites

                for fields in moleculeType.vsites2:
                    atoms = [int(x)-1 for x in fields[:3]]
                    c1 = float(fields[4])
                    vsite = mm.TwoParticleAverageSite(baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], (1-c1), c1)
                    sys.setVirtualSite(baseAtomIndex+atoms[0], vsite)
                for fields in moleculeType.vsites3:
                    atoms = [int(x)-1 for x in fields[:4]]
                    vsiteType = fields[4]
                    c1 = float(fields[5])
                    c2 = float(fields[6])
                    if vsiteType == '1':
                        vsite = mm.ThreeParticleAverageSite(baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], 1-c1-c2, c1, c2)
                    elif vsiteType == '4':
                        c3 = float(fields[7])
                        vsite = mm.OutOfPlaneSite(baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], c1, c2, c3)
                    else:
                        raise ValueError('Internal error: vsites3 has unexpected type: '+vsiteType)
                    sys.setVirtualSite(baseAtomIndex+atoms[0], vsite)

                # Add explicitly specified constraints.

                for fields in moleculeType.constraints:
                    atoms = [int(x)-1 for x in fields[:2]]
                    length = float(fields[3])
                    sys.addConstraint(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], length)

        # Create nonbonded exceptions.

        if not has_nbfix_terms:
            if nb is not None:
                nb.createExceptionsFromBonds(bondIndices, fudgeQQ, fudgeLJ)
            else:
                # Softcore mode: exclusions will be added in the loop below
                # No additional action needed here
                pass
        else:
            # For SOFTCORE/GAUSSIAN/STANDARD_CUSTOM modes or when has_nbfix_terms=False,
            # atom_partners is None, so skip the nbfix-specific exclusion handling.
            # In these modes, exclusions are handled directly by the CustomNonbondedForce
            # or through the pair processing above.
            excluded_atom_pairs = set()
            # Skip nbfix-specific torsion and exclusion processing for these modes
            pass

        # For SOFTCORE/GAUSSIAN/STANDARD_CUSTOM modes without nbfix,
        # add 1-4 pair exclusions directly to CustomNonbondedForce
        if has_nbfix_terms and atom_partners is None and nonbonded_type in (NONBONDED_GAUSSIAN, NONBONDED_SOFTCORE, NONBONDED_STANDARD_CUSTOM):
            # Add 1-4 pair exclusions directly
            for tor in torsionIndices:
                key = min((tor[0], tor[3]), (tor[3], tor[0]))
                if key in excluded_atom_pairs: continue
                excluded_atom_pairs.add(key)
                self.nb_force.addExclusion(tor[0], tor[3])

            # Add bond and angle exclusions
            for bond in bondIndices:
                if bond[0] < bond[1]:
                    self.nb_force.addExclusion(bond[0], bond[1])
                    excluded_atom_pairs.add((bond[0], bond[1]))

            # Add excluded atoms
            # Skip this block when atom_partners is None (SOFTCORE/GAUSSIAN modes without nbfix)
            if atom_partners is not None:
                for atom_idx, atom in atom_partners.items():
                    # Exclude all bonds and angles
                    for atom2 in atom['bond']:
                        if atom2 > atom_idx:
                            if nb is not None:
                                nb.addException(atom_idx, atom2, 0.0, 1.0, 0.0)
                            else:
                                self.nb_force.addExclusion(atom_idx, atom2)
                            excluded_atom_pairs.add((atom_idx, atom2))
                    for atom2 in atom['angle']:
                        if ((atom_idx, atom2) in excluded_atom_pairs):
                            continue
                        if atom2 > atom_idx:
                            if nb is not None:
                                nb.addException(atom_idx, atom2, 0.0, 1.0, 0.0)
                            else:
                                self.nb_force.addExclusion(atom_idx, atom2)
                            excluded_atom_pairs.add((atom_idx, atom2))
                    for atom2 in atom['dihedral']:
                        if atom2 <= atom_idx: continue
                        if ((atom_idx, atom2) in excluded_atom_pairs):
                            continue
                        if nb is not None:
                            nb.addException(atom_idx, atom2, 0.0, 1.0, 0.0)
                        else:
                            self.nb_force.addExclusion(atom_idx, atom2)

        for exclusion in exclusions:
            if nb is not None:
                nb.addException(exclusion[0], exclusion[1], 0.0, 1.0, 0.0, True)
            else:
                self.nb_force.addExclusion(exclusion[0], exclusion[1])

        if lj is not None:
            # We're using a CustomNonbondedForce for LJ interactions, so also create a CustomBondForce
            # to handle the exceptions.

            pair_bond = mm.CustomBondForce('-C/r^6+A/r^12')
            pair_bond.addPerBondParameter('C')
            pair_bond.addPerBondParameter('A')
            pair_bond.setName('LennardJonesExceptions')
            sys.addForce(pair_bond)
            for pair in pairs:
                nb.addException(pair[0], pair[1], pair[2], 1.0, 0.0, True)
                pair_bond.addBond(pair[0], pair[1], [pair[3], pair[4]])
            for i in range(nb.getNumExceptions()):
                ii, jj, q, eps, sig = nb.getExceptionParameters(i)
                lj.addExclusion(ii, jj)
        elif nonbonded_type == NONBONDED_GAUSSIAN:
            # For gaussian mode, also add exclusions to the CustomNonbondedForce
            for i in range(nb.getNumExceptions()):
                ii, jj, q, eps, sig = nb.getExceptionParameters(i)
                self.nb_force.addExclusion(ii, jj)
        elif nonbonded_type == NONBONDED_SOFTCORE:
            # For softcore mode: also add exclusions to the CustomNonbondedForce
            # Check if nb is not None before getting exceptions
            if nb is not None:
                for i in range(nb.getNumExceptions()):
                    ii, jj, q, eps, sig = nb.getExceptionParameters(i)
                    self.nb_force.addExclusion(ii, jj)
        elif nonbonded_type == NONBONDED_STANDARD_CUSTOM:
            # For standard_custom mode: create CustomBondForce for 1-4 pair interactions
            # using standard LJ+Coulomb formula (same as NonbondedForce exceptions)
            
            if len(pairs) > 0:
                # Create energy expression using standard formula
                energy_expr = _create_standard_pair_force_expression()
                
                # Create CustomBondForce with standard LJ+Coulomb formula
                pair_force = mm.CustomBondForce(energy_expr)
                pair_force.setName('StandardCustomPairs')
                
                # Add global parameters
                pair_force.addGlobalParameter("ONE_4PI_EPS0", 138.935456)  # kJ/mol/nm/e^2
                
                # Add per-bond parameters: [charge_prod, sigma, epsilon]
                pair_force.addPerBondParameter("charge_prod")  # Already scaled by fudgeQQ
                pair_force.addPerBondParameter("sigma")        # Combined sigma
                pair_force.addPerBondParameter("epsilon")      # Combined and scaled epsilon
                
                # Add bonds for each pair
                for pair in pairs:
                    i, j = pair[0], pair[1]
                    
                    # Get particle parameters
                    q1, sigma1, eps1 = particle_params[i]
                    q2, sigma2, eps2 = particle_params[j]
                    
                    # Calculate pair parameters with fudge factors (same as softcore)
                    charge_prod = q1 * q2 * fudgeQQ  # Scaled charge product
                    sigma_pair = 0.5 * (sigma1 + sigma2)  # Combined sigma
                    epsilon_pair = math.sqrt(eps1 * eps2) * fudgeLJ  # Combined and scaled epsilon
                    
                    # Add bond with parameters
                    pair_force.addBond(i, j, [charge_prod, sigma_pair, epsilon_pair])
                
                # Add force to system
                sys.addForce(pair_force)
                
                # Add 1-4 pairs to CustomNonbondedForce exclusion list
                # This prevents double-counting of 1-4 interactions
                for pair in pairs:
                    self.nb_force.addExclusion(pair[0], pair[1])
        elif nonbonded_type == NONBONDED_SOFTCORE:
            # For softcore mode: create CustomBondForce for 1-4 pair interactions
            # Pairs are excluded from CustomNonbondedForce, so we need a separate force
            # that uses the full Gapsys linearized formula to match the nonbonded force
            
            if len(pairs) > 0:
                # Get softcore parameters from the nb_force (they were set when creating the force)
                # We need to extract lambda_val, alpha_lj, alpha_coul, sigma_coul from the force
                # But we can't query them, so we'll use the same values that were passed
                # These should match what was used in _create_gapsys_linearized_nb_force
                
                # Create energy expression using helper function
                energy_expr = _create_gapsys_pair_force_expression(has_nbfix_terms=has_nbfix_terms)
                
                # Create CustomBondForce with Gapsys formula
                pair_force = mm.CustomBondForce(energy_expr)
                pair_force.setName('GapsysSoftcorePairs')
                
                # Add global parameters (must match _create_gapsys_linearized_nb_force)
                # We need to get these from the system or use the same defaults
                # Since we can't query the nb_force, we'll use the values from createSystem parameters
                # But we don't have direct access here, so we need to store them or pass them
                # For now, we'll use reasonable defaults and note that they should match
                pair_force.addGlobalParameter("lambda_val", soft_lambda)
                pair_force.addGlobalParameter("alpha_lj", soft_alpha_lj)
                pair_force.addGlobalParameter("alpha_coul", soft_alpha_coul)
                pair_force.addGlobalParameter("sigma_coul", soft_sigma_coul)
                pair_force.addGlobalParameter("ONE_4PI_EPS0", 138.935456)  # kJ/mol/nm/e^2
                
                if has_nbfix_terms:
                    # NBFIX mode: Add per-bond parameters [charge_prod, c6, c12]
                    pair_force.addPerBondParameter("charge_prod")  # Already scaled by fudgeQQ
                    pair_force.addPerBondParameter("c6")           # C6 parameter (scaled)
                    pair_force.addPerBondParameter("c12")          # C12 parameter (scaled)
                    
                    # Need to build NBFIX lookup for pair types
                    # We need access to the acoef/bcoef tables or recompute C6/C12 for each pair
                    # For now, we'll compute C6/C12 directly from atom types
                    
                    # Add bonds for each pair
                    for pair in pairs:
                        i, j = pair[0], pair[1]
                        
                        # Get atom types for this pair
                        atom_type_i = atom_types[i]
                        atom_type_j = atom_types[j]
                        
                        # Get charges
                        q_i = atom_charges[i] if i < len(atom_charges) else 0.0
                        q_j = atom_charges[j] if j < len(atom_charges) else 0.0
                        charge_prod = q_i * q_j * fudgeQQ
                        
                        # Look up C6/C12 from nonbond_params or compute from atom types
                        try:
                            # Try to find NBFIX parameters
                            nbfix_key = tuple(sorted((atom_type_i, atom_type_j)))
                            if nbfix_key in self._nonbondTypes:
                                types = self._nonbondTypes[nbfix_key]
                                params = (float(types[3]), float(types[4]))
                                if self._defaults[1] == '2':
                                    c6 = 4 * params[1] * params[0]**6
                                    c12 = 4 * params[1] * params[0]**12
                                else:
                                    c6 = params[0]
                                    c12 = params[1]
                            else:
                                # Use combination rules
                                atom_i = self._atomTypes[atom_type_i]
                                atom_j = self._atomTypes[atom_type_j]
                                params1 = (float(atom_i[6]), float(atom_i[7]))
                                params2 = (float(atom_j[6]), float(atom_j[7]))
                                if self._defaults[1] == '1':
                                    c6 = math.sqrt(params1[0] * params2[0])
                                    c12 = math.sqrt(params1[1] * params2[1])
                                else:
                                    if self._defaults[1] == '2':
                                        sigma = (params1[0] + params2[0]) / 2
                                    else:
                                        sigma = math.sqrt(params1[0] * params2[0])
                                    epsilon = math.sqrt(params1[1] * params2[1])
                                    c6 = 4 * epsilon * sigma**6
                                    c12 = 4 * epsilon * sigma**12
                            
                            # Apply fudge factor to LJ parameters
                            c6 *= fudgeLJ
                            c12 *= fudgeLJ
                            
                        except (KeyError, IndexError) as e:
                            # Fallback: use zero LJ parameters
                            c6 = 0.0
                            c12 = 0.0
                        
                        # Add bond with parameters
                        pair_force.addBond(i, j, [charge_prod, c6, c12])
                else:
                    # Non-NBFIX mode: Add per-bond parameters [charge_prod, sigma, epsilon]
                    pair_force.addPerBondParameter("charge_prod")  # Already scaled by fudgeQQ
                    pair_force.addPerBondParameter("sigma")        # Combined sigma
                    pair_force.addPerBondParameter("epsilon")      # Combined and scaled epsilon
                    
                    # Add bonds for each pair
                    for pair in pairs:
                        i, j = pair[0], pair[1]
                        
                        # Get particle parameters
                        q1, sigma1, eps1 = particle_params[i]
                        q2, sigma2, eps2 = particle_params[j]
                        
                        # Calculate pair parameters with fudge factors
                        charge_prod = q1 * q2 * fudgeQQ  # Scaled charge product
                        sigma_pair = 0.5 * (sigma1 + sigma2)  # Combined sigma (unchanged)
                        epsilon_pair = math.sqrt(eps1 * eps2) * fudgeLJ  # Combined and scaled epsilon
                        
                        # Add bond with parameters
                        pair_force.addBond(i, j, [charge_prod, sigma_pair, epsilon_pair])
                
                # Add force to system
                sys.addForce(pair_force)
                
                # CRITICAL FIX: Add 1-4 pairs to CustomNonbondedForce exclusion list
                # This prevents double-counting of 1-4 interactions
                # (similar to GAUSSIAN mode handling at lines 1411-1415)
                for pair in pairs:
                    self.nb_force.addExclusion(pair[0], pair[1])
        elif self._defaults[1] == '2':
            for pair in pairs:
                nb.addException(pair[0], pair[1], pair[2], pair[3], pair[4], True)

        # Finish configuring the NonbondedForce.

        # Hardcoded to CutoffPeriodic
        methodMap = {ff.CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic}
        if nb is not None:
            nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            nb.setCutoffDistance(nonbondedCutoff)
            nb.setEwaldErrorTolerance(ewaldErrorTolerance)
            if switchDistance is not None:
                nb.setUseSwitchingFunction(True)
                nb.setSwitchingDistance(switchDistance)
        if lj is not None:
            # Hardcoded to CutoffPeriodic
            lj.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            lj.setCutoffDistance(nonbondedCutoff)
            # Long range correction disabled for CutoffPeriodic
            if switchDistance is not None:
                lj.setUseSwitchingFunction(True)
                lj.setSwitchingDistance(switchDistance)
            lj.setName('LennardJonesForce')

        if has_nbfix_terms:
            # NBFIX processing: Build type tables and assign particle types
            if nonbonded_type == NONBONDED_SOFTCORE:
                # SOFTCORE mode with NBFIX: Build acoef/bcoef tables for CustomNonbondedForce
                atom_nbfix_types = set([])
                for pair in self._nonbondTypes:
                    atom_nbfix_types.add(pair[0])
                    atom_nbfix_types.add(pair[1])

                lj_idx_list = [0 for _ in atom_types]
                lj_radii, lj_depths = [], []
                atom_params_list = []
                num_lj_types = 0
                lj_type_list = []
                
                for i, atom_type in enumerate(atom_types):
                    atom = self._atomTypes[atom_type]
                    if lj_idx_list[i]: continue  # already assigned
                    ljtype = (float(atom[6]), float(atom[7]))
                    atom_params_list.append(ljtype)
                    num_lj_types += 1
                    lj_idx_list[i] = num_lj_types
                    lj_type_list.append(atom)
                    for j in range(i+1, len(atom_types)):
                        atom_type2 = atom_types[j]
                        if lj_idx_list[j] > 0: continue  # already assigned
                        atom2 = self._atomTypes[atom_type2]
                        ljtype2 = (float(atom2[6]), float(atom2[7]))
                        if atom2 is atom:
                            lj_idx_list[j] = num_lj_types
                        elif atom_type not in atom_nbfix_types:
                            # Only non-NBFIXed atom types can be compressed
                            if ljtype == ljtype2:
                                lj_idx_list[j] = num_lj_types

                # Now everything is assigned. Create the A-coefficient and B-coefficient arrays
                acoef = [0.0 for i in range(num_lj_types * num_lj_types)]
                bcoef = [0.0 for i in range(num_lj_types * num_lj_types)]
                
                for i in range(num_lj_types):
                    namei = lj_type_list[i][0]
                    for j in range(num_lj_types):
                        namej = lj_type_list[j][0]
                        try:
                            types = self._nonbondTypes[tuple(sorted((namei, namej)))]
                            params = (float(types[3]), float(types[4]))
                            if self._defaults[1] == '2':
                                c6 = 4 * params[1] * params[0]**6
                                c12 = 4 * params[1] * params[0]**12
                            else:
                                c6 = params[0]
                                c12 = params[1]
                        except KeyError:
                            params1 = atom_params_list[i]
                            params2 = atom_params_list[j]
                            if self._defaults[1] == '1':
                                c6 = math.sqrt(params1[0]*params2[0])
                                c12 = math.sqrt(params1[1]*params2[1])
                            else:
                                if self._defaults[1] == '2':
                                    sigma = (params1[0] + params2[0]) / 2
                                else:
                                    sigma = math.sqrt(params1[0] * params2[0])
                                epsilon = math.sqrt(params1[1] * params2[1])
                                c6 = 4 * epsilon * sigma**6
                                c12 = 4 * epsilon * sigma**12
                        acoef[i + num_lj_types*j] = math.sqrt(c12)
                        bcoef[i + num_lj_types*j] = c6
                
                # Add tabulated functions to the softcore force
                self.nb_force.addTabulatedFunction('acoef', mm.Discrete2DFunction(num_lj_types, num_lj_types, acoef))
                self.nb_force.addTabulatedFunction('bcoef', mm.Discrete2DFunction(num_lj_types, num_lj_types, bcoef))
                
                # Update particle parameters with correct type indices
                for i, idx in enumerate(lj_idx_list):
                    # particle_params[i] = [type, charge] (placeholder was set earlier)
                    type_idx = idx - 1  # adjust for indexing from 0
                    charge = atom_charges[i] if i < len(atom_charges) else 0.0
                    self.nb_force.setParticleParameters(i, [type_idx, charge])
                    # Update particle_params for pair processing
                    if i < len(particle_params):
                        particle_params[i] = [type_idx, charge]
            
            elif nonbonded_type in (NONBONDED_GAUSSIAN, NONBONDED_STANDARD_CUSTOM):
                # These modes: skip nbfix processing entirely (they don't use LJ parameters)
                pass
            
            else:
                # Standard mode with NBFIX: Use existing logic for NonbondedForce + CustomNonbondedForce
                atom_nbfix_types = set([])
                for pair in self._nonbondTypes:
                    atom_nbfix_types.add(pair[0])
                    atom_nbfix_types.add(pair[1])

                lj_idx_list = [0 for _ in atom_types]
                lj_radii, lj_depths = [], []
                atom_params_list = []
                num_lj_types = 0
                lj_type_list = []
                for i,atom_type in enumerate(atom_types):
                    atom = self._atomTypes[atom_type]
                    if lj_idx_list[i]: continue # already assigned
                    ljtype = (float(atom[6]), float(atom[7]))
                    atom_params_list.append(ljtype)
                    num_lj_types += 1
                    lj_idx_list[i] = num_lj_types
                    lj_type_list.append(atom)
                    for j in range(i+1, len(atom_types)):
                        atom_type2 = atom_types[j]
                        if lj_idx_list[j] > 0: continue # already assigned
                        atom2 = self._atomTypes[atom_type2]
                        ljtype2 = (float(atom2[6]), float(atom2[7]))
                        if atom2 is atom:
                            lj_idx_list[j] = num_lj_types
                        elif atom_type not in atom_nbfix_types:
                            # Only non-NBFIXed atom types can be compressed
                            if ljtype == ljtype2:
                                lj_idx_list[j] = num_lj_types

                # Now everything is assigned. Create the A-coefficient and
                # B-coefficient arrays
                acoef = [0 for i in range(num_lj_types*num_lj_types)]
                bcoef = acoef[:]
                for i in range(num_lj_types):
                    namei = lj_type_list[i][0]
                    for j in range(num_lj_types):
                        namej = lj_type_list[j][0]
                        try:
                            types = self._nonbondTypes[tuple(sorted((namei, namej)))]
                            params = (float(types[3]), float(types[4]))
                            if self._defaults[1] == '2':
                                c6 = 4 * params[1] * params[0]**6
                                c12 = 4 * params[1] * params[0]**12
                            else:
                                c6 = params[0]
                                c12 = params[1]
                        except KeyError:
                            params1 = atom_params_list[i]
                            params2 = atom_params_list[j]
                            if self._defaults[1] == '1':
                                c6 = math.sqrt(params1[0]*params2[0])
                                c12 = math.sqrt(params1[1]*params2[1])
                            else:
                                if self._defaults[1] == '2':
                                    sigma = (params1[0] + params2[0]) / 2
                                else:
                                    sigma = math.sqrt(params1[0] + params2[0])
                                epsilon = math.sqrt(params1[1] * params2[1])
                                c6 = 4 * epsilon * sigma**6
                                c12 = 4 * epsilon * sigma**12
                        acoef[i+num_lj_types*j] = math.sqrt(c12)
                        bcoef[i+num_lj_types*j] = c6
                lj.addTabulatedFunction('acoef', mm.Discrete2DFunction(num_lj_types, num_lj_types, acoef))
                lj.addTabulatedFunction('bcoef', mm.Discrete2DFunction(num_lj_types, num_lj_types, bcoef))
                for i, idx in enumerate(lj_idx_list):
                    lj.setParticleParameters(i, [idx-1]) # adjust for indexing from 0

        # Adjust masses.

        if hydrogenMass is not None:
            for atom1, atom2 in self.topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if rigidWater and atom2.residue.name == 'HOH':
                    continue
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Add a CMMotionRemover.

        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())
        return sys

def _defaultGromacsIncludeDir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) for the gromacs binary
    'pdb2gmx' or 'gmx' in the PATH, or (3) just using the default gromacs
    install location, /usr/local/gromacs/share/gromacs/top """
    if 'GMXDATA' in os.environ:
        return os.path.join(os.environ['GMXDATA'], 'top')
    if 'GMXBIN' in os.environ:
        return os.path.abspath(os.path.join(os.environ['GMXBIN'], '..', 'share', 'gromacs', 'top'))

    pdb2gmx_path = shutil.which('pdb2gmx')
    if pdb2gmx_path is not None:
        return os.path.abspath(os.path.join(os.path.dirname(pdb2gmx_path), '..', 'share', 'gromacs', 'top'))
    else:
        gmx_path = shutil.which('gmx')
        if gmx_path is not None:
            return os.path.abspath(os.path.join(os.path.dirname(gmx_path), '..', 'share', 'gromacs', 'top'))

    return '/usr/local/gromacs/share/gromacs/top'

# ============================================================================
# Three-stage optimization methods
# ============================================================================

def _create_gaussian_nb_force(nonbonded_cutoff, gaussian_width=0.085):
    """Create Gaussian repulsion force for initial untangling (Stage 1).

    Formula: E = ga_k * ga_h * exp(-(r/ga_w)^2)

    This provides pure geometric repulsion without chemical information,
    useful for quickly removing atomic overlaps and untangling knots.

    Parameters
    ----------
    nonbonded_cutoff : Quantity
        Cutoff distance for nonbonded interactions
    gaussian_width : float=0.085
        Width parameter (ga_w) in nm. Default 0.085 nm (0.85 Å).

    Returns
    -------
    CustomNonbondedForce
        Gaussian repulsion force object
    """
    nb_force = mm.CustomNonbondedForce(
        "ga_k * ga_h * exp(-(r/ga_w)^2);"
    )
    nb_force.addGlobalParameter("ga_k", 1.0)   # ON/OFF switch
    nb_force.addGlobalParameter("ga_h", 800.0) # Height (kJ/mol)
    nb_force.addGlobalParameter("ga_w", gaussian_width)  # Width (nm)
    nb_force.addPerParticleParameter("dummy_type")
    nb_force.addPerParticleParameter("dummy_q")
    nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
    return nb_force

def _add_gbsa_solvent(sys, top, gb_model='GBn2', salt_conc=0.0, nonbondedCutoff=None, forcefield_type='AMBER'):
    """Add GBSA implicit solvent to the system.

    Parameters
    ----------
    sys : mm.System
        OpenMM system object
    top : GromacsTopFileWithSoftcore
        Topology object with _molecules and _moleculeTypes attributes
    gb_model : str='GBn2'
        GBSA model: 'OBC2' or 'GBn2'
    salt_conc : float=0.0
        Salt concentration (M), for Debye-Hückel screening
    nonbondedCutoff : Quantity=None
        Cutoff distance for nonbonded interactions. If None, uses 2.0 nm default.
    forcefield_type : str='AMBER'
        Force field type: 'AMBER' or 'CHARMM'. Determines which atom type mapping to use for GBSA parameters.
    """
    from math import sqrt

    # Validate and select forcefield mapping
    if forcefield_type.upper() not in ['AMBER', 'CHARMM']:
        raise ValueError(f"Unsupported forcefield_type: {forcefield_type}. "
                        f"Supported types: 'AMBER', 'CHARMM'")

    # Select the appropriate atom type mapping based on forcefield
    if forcefield_type.upper() == 'AMBER':
        atom_mapping = AMBER99SB_ILDN_ATOM_MAPPING
    else:  # CHARMM
        atom_mapping = CHARMM36_ATOM_MAPPING

    # Validate gb_model
    if gb_model not in ['OBC2', 'GBn2']:
        raise ValueError(f"Unsupported GB model: {gb_model}. "
                       f"Supported models: OBC2, GBn2")

    # Set dielectric constants
    solute_dielectric = 1.0
    solvent_dielectric = 78.5

    # Set cutoff distance (use 2.0 nm as default if not provided)
    if nonbondedCutoff is not None:
        if unit.is_quantity(nonbondedCutoff):
            cutoff = nonbondedCutoff.value_in_unit(unit.nanometer)
        else:
            cutoff = float(nonbondedCutoff)
    else:
        cutoff = 2.0  # nm default

    # Compute kappa from salt concentration if needed
    kappa = 0.0
    if salt_conc > 0:
        # The constant matches Amber's kappa conversion factor
        kappa = 50.33355 * sqrt(salt_conc / solvent_dielectric / 300.0)  # Assuming 300K
        # Multiply by 0.73 to account for ion exclusions, and multiply by 10
        # to convert to 1/nm from 1/angstroms
        kappa *= 7.3

    # Create the appropriate GB force
    if gb_model == 'OBC2':
        # Use native GBSAOBCForce
        gb_force = mm.GBSAOBCForce()
        gb_force.setSurfaceAreaEnergy(0.0)  # No SASA contribution
        gb_force.setSoluteDielectric(solute_dielectric)
        gb_force.setSolventDielectric(solvent_dielectric)

        # Set nonbonded method to CutoffPeriodic
        try:
            gb_force.setNonbondedMethod(mm.GBSAOBCForce.CutoffPeriodic)
            gb_force.setCutoffDistance(cutoff)
        except Exception:
            pass  # Some methods may not be available

        # Set kappa if salt is present
        if kappa > 0 and hasattr(gb_force, 'setKappa'):
            gb_force.setKappa(kappa)

    elif gb_model == 'GBn2':
        # Use custom GBn2 force from OpenMM internal
        gb_force = GBSAGBn2Force(
            solventDielectric=solvent_dielectric,
            soluteDielectric=solute_dielectric,
            SA=None,  # No surface area contribution
            cutoff=cutoff,
            kappa=kappa
        )
        # Explicitly set CutoffPeriodic to match NonbondedForce
        gb_force.setNonbondedMethod(GBSAGBn2Force.CutoffPeriodic)
        gb_force.setCutoffDistance(cutoff)

    # Get charges from the topology
    # Build a list of charges from the atom information
    charges = []
    for moleculeName, moleculeCount in top._molecules:
        moleculeType = top._moleculeTypes[moleculeName]
        for _ in range(moleculeCount):
            for atom_fields in moleculeType.atoms:
                # atom_fields format: [nr, type, resnr, residue, atom, cgnr, charge, mass]
                if len(atom_fields) > 6:
                    charge = float(atom_fields[6])
                else:
                    charge = 0.0
                charges.append(charge)

    # Add particles to the GB force
    if gb_model == 'OBC2':
        # Use atom_mapping for both radius and screening factors
        for i, atom in enumerate(top.topology.atoms()):
            if i < len(charges):
                charge = charges[i]
            else:
                charge = 0.0

            # Get atom type and look up in mapping
            atom_type = atom.name
            if atom_type in atom_mapping:
                atom_info = atom_mapping[atom_type]
                radius = atom_info['radius']
                screening = atom_info.get('screen', 0.8)
            elif atom.element and atom.element.symbol in ['H', 'C', 'N', 'O', 'S', 'P']:
                # Fallback to element-based defaults
                element_defaults = {
                    'H': (0.100, 0.85), 'C': (0.170, 0.72),
                    'N': (0.155, 0.79), 'O': (0.150, 0.85),
                    'S': (0.180, 0.96), 'P': (0.185, 0.86)
                }
                defaults = element_defaults.get(atom.element.symbol, (0.170, 0.8))
                radius, screening = defaults
            else:
                radius, screening = 0.170, 0.8

            gb_force.addParticle(charge, radius, screening)

    elif gb_model == 'GBn2':
        # For GBn2, get parameters from getStandardParameters
        gb_parms = GBSAGBn2Force.getStandardParameters(top.topology)

        # Add particles: [charge, or, sr, alpha, beta, gamma]
        # or = offset radius = radius + 0.0195141 nm
        OFFSET = 0.0195141

        for i, atom in enumerate(top.topology.atoms()):
            if i < len(charges) and i < len(gb_parms):
                charge = charges[i]
                radius = gb_parms[i][0]  # This is in nm
                sr = gb_parms[i][1]
                alpha = gb_parms[i][2]
                beta = gb_parms[i][3]
                gamma = gb_parms[i][4]

                # Compute offset radius
                or_value = radius + OFFSET

                # Add particle with parameters: [charge, or, sr, alpha, beta, gamma]
                gb_force.addParticle([charge, or_value, sr, alpha, beta, gamma])
            else:
                # Fallback for missing data
                gb_force.addParticle([0.0, 0.19, 0.8, 1.0, 1.0, 1.0])

    # Finalize and add to system
    if hasattr(gb_force, 'finalize') and gb_model == 'GBn2':
        gb_force.finalize()

    sys.addForce(gb_force)

def _create_standard_custom_nb_force(nonbonded_cutoff, switchDistance=None):
    """Create standard LJ+Coulomb CustomNonbondedForce for debugging/comparison.
    
    This implements the exact same formula as NonbondedForce but using
    CustomNonbondedForce to match the architecture of softcore mode.
    
    Parameters
    ----------
    nonbonded_cutoff : Quantity
        Cutoff distance (e.g. 1.0 * unit.nanometer)
    switchDistance : Quantity or None
        Switching distance for smoothing
    
    Returns
    -------
    nb_force : openmm.CustomNonbondedForce
        CustomNonbondedForce implementing standard LJ+Coulomb potential
    """
    # Standard LJ + Coulomb formula
    # LJ: 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
    # Coulomb: (1/4πε0) * q1*q2/r = 138.935456 * q1*q2/r (kJ/mol/nm/e^2)
    
    energy_expr = """
4.0*sqrt(epsilon1*epsilon2)*((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6)
+ ONE_4PI_EPS0*q1*q2/r
"""
    
    nb_force = mm.CustomNonbondedForce(energy_expr)
    
    # Global parameters
    nb_force.addGlobalParameter("ONE_4PI_EPS0", 138.935456)  # kJ/mol/nm/e^2
    
    # Per-particle parameters
    nb_force.addPerParticleParameter("q")       # Charge
    nb_force.addPerParticleParameter("sigma")   # LJ sigma
    nb_force.addPerParticleParameter("epsilon") # LJ epsilon
    
    # Set cutoff
    nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
    
    # Use switching function if provided
    if switchDistance is not None:
        nb_force.setUseSwitchingFunction(True)
        if unit.is_quantity(switchDistance):
            nb_force.setSwitchingDistance(switchDistance.value_in_unit(unit.nanometer))
        else:
            nb_force.setSwitchingDistance(switchDistance)
    else:
        nb_force.setUseSwitchingFunction(False)
    
    # Disable LRC (to match standard mode behavior with CutoffPeriodic)
    nb_force.setUseLongRangeCorrection(False)
    
    return nb_force

def _create_standard_pair_force_expression():
    """Create energy expression for standard 1-4 pair interactions.
    
    This implements the standard LJ+Coulomb formula for CustomBondForce,
    matching what NonbondedForce does for exceptions.
    
    For CustomBondForce, we use per-bond parameters:
    - charge_prod: q1 * q2 (already scaled by fudgeQQ)
    - sigma: 0.5 * (sigma1 + sigma2)
    - epsilon: sqrt(epsilon1 * epsilon2) (already scaled by fudgeLJ)
    
    Returns
    -------
    str
        Energy expression string for CustomBondForce
    """
    # Standard LJ + Coulomb (parameters are already combined and scaled)
    energy_expr = """
4.0*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*charge_prod/r
"""
    return energy_expr

def _create_gapsys_pair_force_expression(has_nbfix_terms=False):
    """Create energy expression for Gapsys softcore pair interactions (1-4 pairs).
    
    This function generates the same Gapsys linearized formula used in
    _create_gapsys_linearized_nb_force, but adapted for CustomBondForce.
    
    For CustomBondForce, we use per-bond parameters:
    - charge_prod: q1 * q2 (already scaled by fudgeQQ)
    - sigma: 0.5 * (sigma1 + sigma2) OR derived from C6/C12 for NBFIX
    - epsilon: sqrt(epsilon1 * epsilon2) (already scaled by fudgeLJ) OR derived from C6/C12 for NBFIX
    
    Parameters
    ----------
    has_nbfix_terms : bool=False
        If True, use C6/C12-based formula compatible with NBFIX parameters
    
    Returns
    -------
    str
        Energy expression string for CustomBondForce
    """
    if has_nbfix_terms:
        # NBFIX mode: Parameters are [charge_prod, c6, c12]
        # Formula: E_LJ = C12/r^12 - C6/r^6
        # Need to compute sigma_eff and epsilon_eff for softcore switching
        
        # For NBFIX pairs, we store C6 and C12 directly (scaled by fudge factors)
        # sigma_eff = (C12/C6)^(1/6)
        # epsilon_eff = C6^2 / (4*C12)
        
        q_prod = "charge_prod"
        soft = "(1.0 - lambda_val)^0.1666667"
        c6_expr = "c6"
        c12_expr = "c12"
        
        # Compute sigma_eff and epsilon_eff
        sig_eff = f"(({c12_expr})/({c6_expr}))^0.1666667"
        eps_eff = f"({c6_expr})*({c6_expr})/(4.0*({c12_expr}))"
        
        # LJ switching parameters
        rsw_lj_inline = f"{sig_eff} * alpha_lj * 1.244 * {soft}"
        rsw_lj_safe_inline = f"max({rsw_lj_inline}, 1.0e-6)"
        rsw_lj_sq_inline = f"{rsw_lj_safe_inline}^2"
        u_lj_inline = f"{sig_eff} / {rsw_lj_safe_inline}"
        V_lj_sw_inline = f"4.0 * {eps_eff} * ({u_lj_inline}^12 - {u_lj_inline}^6)"
        F_lj_sw_inline = f"(24.0 * {eps_eff} / {rsw_lj_safe_inline}) * (2.0 * {u_lj_inline}^12 - {u_lj_inline}^6)"
        dF_lj_sw_inline = f"(24.0 * {eps_eff} / {rsw_lj_sq_inline}) * (-26.0 * {u_lj_inline}^12 + 7.0 * {u_lj_inline}^6)"
        
        # Standard LJ from C6/C12
        lj_standard = f"({c12_expr})/(r^12) - ({c6_expr})/(r^6)"
        
        # Coulomb parameters
        rsw_q_inline = f"alpha_coul * (1.0 + sigma_coul * abs({q_prod})) * {soft}"
        rsw_q_safe_inline = f"max({rsw_q_inline}, 1.0e-6)"
        rsw_q_sq_inline = f"{rsw_q_safe_inline}^2"
        C_pre_inline = f"ONE_4PI_EPS0 * {q_prod}"
        V_q_sw_inline = f"{C_pre_inline} / {rsw_q_safe_inline}"
        F_q_sw_inline = f"{C_pre_inline} / {rsw_q_sq_inline}"
        dF_q_sw_inline = f"-2.0 * {C_pre_inline} / ({rsw_q_sq_inline} * {rsw_q_safe_inline})"
        
        energy_expr = f"""
lambda_val * (
    select(step(r - {rsw_lj_safe_inline}), 
           {lj_standard},
           {V_lj_sw_inline} + {F_lj_sw_inline} * ({rsw_lj_inline} - r) - 0.5 * {dF_lj_sw_inline} * ({rsw_lj_inline} - r)^2)
    +
    select(step(r - {rsw_q_safe_inline}), 
           ONE_4PI_EPS0 * {q_prod} / r,
           {V_q_sw_inline} + {F_q_sw_inline} * ({rsw_q_inline} - r) - 0.5 * {dF_q_sw_inline} * ({rsw_q_inline} - r)^2)
);
"""
    else:
        # Standard mode: Parameters are [charge_prod, sigma, epsilon]
        # For CustomBondForce, we have per-bond parameters: charge_prod, sigma, epsilon
        # These are already combined and scaled appropriately
        
        # Helper expressions for inline use
        e = "epsilon"  # Already combined and scaled
        sig = "sigma"  # Already combined
        q_prod = "charge_prod"  # Already scaled by fudgeQQ
        soft = "(1.0 - lambda_val)^0.1666667"
        
        # LJ parameters (inline)
        rsw_lj_inline = f"{sig} * alpha_lj * 1.244 * {soft}"
        rsw_lj_safe_inline = f"max({rsw_lj_inline}, 1.0e-6)"
        rsw_lj_sq_inline = f"{rsw_lj_safe_inline}^2"
        u_lj_inline = f"{sig} / {rsw_lj_safe_inline}"
        V_lj_sw_inline = f"4.0 * {e} * ({u_lj_inline}^12 - {u_lj_inline}^6)"
        F_lj_sw_inline = f"(24.0 * {e} / {rsw_lj_safe_inline}) * (2.0 * {u_lj_inline}^12 - {u_lj_inline}^6)"
        dF_lj_sw_inline = f"(24.0 * {e} / {rsw_lj_sq_inline}) * (-26.0 * {u_lj_inline}^12 + 7.0 * {u_lj_inline}^6)"
        
        # Coulomb parameters (inline)
        rsw_q_inline = f"alpha_coul * (1.0 + sigma_coul * abs({q_prod})) * {soft}"
        rsw_q_safe_inline = f"max({rsw_q_inline}, 1.0e-6)"
        rsw_q_sq_inline = f"{rsw_q_safe_inline}^2"
        C_pre_inline = f"ONE_4PI_EPS0 * {q_prod}"
        V_q_sw_inline = f"{C_pre_inline} / {rsw_q_safe_inline}"
        F_q_sw_inline = f"{C_pre_inline} / {rsw_q_sq_inline}"
        dF_q_sw_inline = f"-2.0 * {C_pre_inline} / ({rsw_q_sq_inline} * {rsw_q_safe_inline})"
        
        energy_expr = f"""
lambda_val * (
    select(step(r - {rsw_lj_safe_inline}), 
           4.0 * {e} * ({sig}/r)^12 - 4.0 * {e} * ({sig}/r)^6,
           {V_lj_sw_inline} + {F_lj_sw_inline} * ({rsw_lj_inline} - r) - 0.5 * {dF_lj_sw_inline} * ({rsw_lj_inline} - r)^2)
    +
    select(step(r - {rsw_q_safe_inline}), 
           ONE_4PI_EPS0 * {q_prod} / r,
           {V_q_sw_inline} + {F_q_sw_inline} * ({rsw_q_inline} - r) - 0.5 * {dF_q_sw_inline} * ({rsw_q_inline} - r)^2)
);
"""
    return energy_expr

def _create_gapsys_linearized_nb_force(nonbonded_cutoff, current_lambda=1.0,
                                        alpha_lj=0.85, alpha_coul=0.3, sigma_coul=1.0,
                                        use_implicit_solvent=False, switchDistance=None,
                                        has_nbfix_terms=False):
    """Create Gapsys "New Soft-Core" (Linearized Force) CustomNonbondedForce.

    Based on Gapsys et al., J. Chem. Theory Comput. 2015, 11, 11, 5920–5930
    "Interaction of Legolane with the Outer Membrane of Gram-Negative Bacteria"

    This implementation linearizes both LJ and Coulomb forces near r=0,
    eliminating singularities and enabling stable Energy Minimization.

    Formula Overview:
    - For r > r_sw (switching distance): Use standard LJ/Coulomb
    - For r <= r_sw: Use Taylor-expanded linearized force

    Parameters
    ----------
    nonbonded_cutoff : Quantity
        Cutoff distance (e.g. 1.0 * unit.nanometer)
    current_lambda : float=1.0
        Decoupling parameter (1.0 = fully coupled, 0.0 = decoupled).
        NOTE: For softcore EM, we use a "softness" lambda (default 0.85)
        to control the softening, NOT the decoupling.
    alpha_lj : float=0.85
        LJ soft-core control parameter (Gapsys paper recommends 0.85)
    alpha_coul : float=0.3
        Coulomb soft-core control parameter (Gapsys paper recommends 0.3)
    sigma_coul : float=1.0
        Coulomb soft-core charge scaling parameter (Gapsys paper recommends 1.0)
    use_implicit_solvent : bool=False
        If True, use CutoffPeriodic with GBSA forces
    has_nbfix_terms : bool=False
        If True, use NBFIX-compatible formula with type parameters

    Returns
    -------
    nb_force : openmm.CustomNonbondedForce
        CustomNonbondedForce implementing linearized soft-core potential
    """
    # NOTE: Main expression must use ONLY built-in variables (r, q1, q2, sigma1, sigma2, epsilon1, epsilon2)
    # and global parameters. Custom intermediate variables CANNOT appear in the main expression.
    # All intermediate variables are defined AFTER the main expression for documentation only.
    
    if has_nbfix_terms:
        # NBFIX mode: Use tabulated functions for LJ parameters
        # Formula: E_LJ = C12/r^12 - C6/r^6
        # where C6 = bcoef(type1, type2) and C12 = acoef(type1, type2)^2
        #
        # CRITICAL: OpenMM's backward variable definition (defining variables after using them)
        # only works for INDEPENDENT variables. Variables CANNOT reference other backward-defined
        # variables. This was discovered through extensive testing (tests/test_sonnet/charmm/nbfix_gapsys/).
        #
        # Solution: Only use backward definition for c6 and c12 (which depend only on built-in functions).
        # All other calculations must be inlined into the main expression.
        
        # Helper: inline calculations
        q_prod = "q1 * q2"
        soft_factor = "(1.0 - lambda_val)^0.1666667"
        
        # sig_eff and eps_eff in terms of c6/c12
        sig_eff = "(c12/c6)^0.1666667"
        eps_eff = "(c6*c6)/(4.0*c12)"
        
        # LJ softcore switching distance
        rsw_lj = f"({sig_eff} * alpha_lj * 1.244 * {soft_factor})"
        rsw_lj_safe = f"max({rsw_lj}, 1.0e-6)"
        rsw_lj_sq = f"({rsw_lj_safe} * {rsw_lj_safe})"
        u_lj = f"({sig_eff} / {rsw_lj_safe})"
        
        # LJ softcore parameters at r_sw
        V_lj_sw = f"(4.0 * {eps_eff} * ({u_lj}^12 - {u_lj}^6))"
        F_lj_sw = f"((24.0 * {eps_eff} / {rsw_lj_safe}) * (2.0 * {u_lj}^12 - {u_lj}^6))"
        dF_lj_sw = f"((24.0 * {eps_eff} / {rsw_lj_sq}) * (-26.0 * {u_lj}^12 + 7.0 * {u_lj}^6))"
        
        # Coulomb softcore switching distance
        rsw_q = f"(alpha_coul * (1.0 + sigma_coul * abs({q_prod})) * {soft_factor})"
        rsw_q_safe = f"max({rsw_q}, 1.0e-6)"
        rsw_q_sq = f"({rsw_q_safe} * {rsw_q_safe})"
        C_pre = f"(ONE_4PI_EPS0 * {q_prod})"
        
        # Coulomb softcore parameters at r_sw
        V_q_sw = f"({C_pre} / {rsw_q_safe})"
        F_q_sw = f"({C_pre} / {rsw_q_sq})"
        dF_q_sw = f"(-2.0 * {C_pre} / ({rsw_q_sq} * {rsw_q_safe}))"
        
        # Build fully inlined energy expression
        energy_expr = f"""
lambda_val * (
    select(step(r - {rsw_lj_safe}), 
           c12/(r^12) - c6/(r^6),
           {V_lj_sw} + {F_lj_sw} * ({rsw_lj} - r) - 0.5 * {dF_lj_sw} * ({rsw_lj} - r)^2)
    +
    select(step(r - {rsw_q_safe}), 
           ONE_4PI_EPS0 * {q_prod} / r,
           {V_q_sw} + {F_q_sw} * ({rsw_q} - r) - 0.5 * {dF_q_sw} * ({rsw_q} - r)^2)
);

c6 = bcoef(type1, type2);
c12 = (acoef(type1, type2))^2
"""
        
        nb_force = mm.CustomNonbondedForce(energy_expr)
        
        # Global parameters
        nb_force.addGlobalParameter("lambda_val", current_lambda)
        nb_force.addGlobalParameter("alpha_lj", alpha_lj)
        nb_force.addGlobalParameter("alpha_coul", alpha_coul)
        nb_force.addGlobalParameter("sigma_coul", sigma_coul)
        nb_force.addGlobalParameter("ONE_4PI_EPS0", 138.935456)  # kJ/mol/nm/e^2
        
        # Per-particle parameters for NBFIX mode
        nb_force.addPerParticleParameter("type")   # Type index for NBFIX lookup
        nb_force.addPerParticleParameter("q")      # Charge
        
    else:
        # Standard mode: Use combination rules
        # Helper expressions for inline use
        e = "sqrt(epsilon1 * epsilon2)"
        sig = "0.5 * (sigma1 + sigma2)"
        q_prod = "q1 * q2"
        soft = "(1.0 - lambda_val)^0.1666667"
        
        # LJ parameters (inline)
        rsw_lj_inline = f"{sig} * alpha_lj * 1.244 * {soft}"
        rsw_lj_safe_inline = f"max({rsw_lj_inline}, 1.0e-6)"
        rsw_lj_sq_inline = f"{rsw_lj_safe_inline}^2"
        u_lj_inline = f"{sig} / {rsw_lj_safe_inline}"
        V_lj_sw_inline = f"4.0 * {e} * ({u_lj_inline}^12 - {u_lj_inline}^6)"
        F_lj_sw_inline = f"(24.0 * {e} / {rsw_lj_safe_inline}) * (2.0 * {u_lj_inline}^12 - {u_lj_inline}^6)"
        dF_lj_sw_inline = f"(24.0 * {e} / {rsw_lj_sq_inline}) * (-26.0 * {u_lj_inline}^12 + 7.0 * {u_lj_inline}^6)"
        
        # Coulomb parameters (inline)
        rsw_q_inline = f"alpha_coul * (1.0 + sigma_coul * abs({q_prod})) * {soft}"
        rsw_q_safe_inline = f"max({rsw_q_inline}, 1.0e-6)"
        rsw_q_sq_inline = f"{rsw_q_safe_inline}^2"
        C_pre_inline = f"ONE_4PI_EPS0 * {q_prod}"
        V_q_sw_inline = f"{C_pre_inline} / {rsw_q_safe_inline}"
        F_q_sw_inline = f"{C_pre_inline} / {rsw_q_sq_inline}"
        dF_q_sw_inline = f"-2.0 * {C_pre_inline} / ({rsw_q_sq_inline} * {rsw_q_safe_inline})"
        
        energy_expr = f"""
lambda_val * (
    select(step(r - {rsw_lj_safe_inline}), 
           4.0 * {e} * ({sig}/r)^12 - 4.0 * {e} * ({sig}/r)^6,
           {V_lj_sw_inline} + {F_lj_sw_inline} * ({rsw_lj_inline} - r) - 0.5 * {dF_lj_sw_inline} * ({rsw_lj_inline} - r)^2)
    +
    select(step(r - {rsw_q_safe_inline}), 
           ONE_4PI_EPS0 * {q_prod} / r,
           {V_q_sw_inline} + {F_q_sw_inline} * ({rsw_q_inline} - r) - 0.5 * {dF_q_sw_inline} * ({rsw_q_inline} - r)^2)
);

/* Documentation: LJ parameters (NOT used in main expression above) */
e = {e};
sig = {sig};
q_prod = {q_prod};
soft_factor = {soft};

rsw_lj = alpha_lj * 1.244 * sig * soft_factor;
rsw_lj_safe = max(rsw_lj, 1.0e-6);
rsw_lj_sq = rsw_lj_safe * rsw_lj_safe;
u_lj = sig / rsw_lj_safe;
V_lj_sw  = 4.0 * e * (u_lj^12 - u_lj^6);
F_lj_sw  = (24.0 * e / rsw_lj_safe) * (2.0 * u_lj^12 - u_lj^6);
dF_lj_sw = (24.0 * e / rsw_lj_sq) * (-26.0 * u_lj^12 + 7.0 * u_lj^6);

/* Documentation: Coulomb parameters (NOT used in main expression above) */
rsw_q = alpha_coul * (1.0 + sigma_coul * abs(q_prod)) * soft_factor;
rsw_q_safe = max(rsw_q, 1.0e-6);
rsw_q_sq = rsw_q_safe * rsw_q_safe;
C_pre = ONE_4PI_EPS0 * q_prod;
V_q_sw  = C_pre / rsw_q_safe;
F_q_sw  = C_pre / rsw_q_sq;
dF_q_sw = -2.0 * C_pre / (rsw_q_sq * rsw_q_safe);
"""
        
        nb_force = mm.CustomNonbondedForce(energy_expr)
        
        # Global parameters
        nb_force.addGlobalParameter("lambda_val", current_lambda)
        nb_force.addGlobalParameter("alpha_lj", alpha_lj)
        nb_force.addGlobalParameter("alpha_coul", alpha_coul)
        nb_force.addGlobalParameter("sigma_coul", sigma_coul)
        nb_force.addGlobalParameter("ONE_4PI_EPS0", 138.935456)  # kJ/mol/nm/e^2
        
        # Per-particle parameters (sigma1, epsilon1 for each particle)
        # Note: Using sigma/epsilon from combination rule 2 (sigma/epsilon directly)
        # For combination rule 1/3, we need to convert C6/C12 to sigma/epsilon
        nb_force.addPerParticleParameter("q")      # Charge
        nb_force.addPerParticleParameter("sigma")  # LJ sigma (for rule 2) or derived from C6/C12
        nb_force.addPerParticleParameter("epsilon") # LJ epsilon

    # Set cutoff
    nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
    
    # Use switching function only if switchDistance is provided (matching standard mode behavior)
    if switchDistance is not None:
        nb_force.setUseSwitchingFunction(True)
        if unit.is_quantity(switchDistance):
            nb_force.setSwitchingDistance(switchDistance.value_in_unit(unit.nanometer))
        else:
            nb_force.setSwitchingDistance(switchDistance)
    else:
        nb_force.setUseSwitchingFunction(False)
    
    # LRC is disabled for Gapsys linearized potential
    # Gapsys linearized potential doesn't support analytical long-range correction
    nb_force.setUseLongRangeCorrection(False)
    
    return nb_force
