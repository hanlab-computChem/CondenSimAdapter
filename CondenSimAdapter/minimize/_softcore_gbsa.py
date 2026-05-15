"""Internal module: GBSA implicit solvent helpers for softcore optimization."""

from __future__ import absolute_import

import openmm as mm
import openmm.unit as unit
from openmm.app.internal.customgbforces import GBSAGBn2Force

# AMBER99SB-ILDN atom type mapping for GBSA (Born radii in nm)
# Based on MBondi2/MBondi3 parameterization
AMBER99SB_ILDN_ATOM_MAPPING = {
    # Hydrogen atoms
    "H": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Amide hydrogen
    "H1": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Aliphatic hydrogen
    "HP": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Aromatic hydrogen
    "HC": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Alkyl hydrogen
    "HB1": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Beta hydrogen (MET)
    "HB2": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Beta hydrogen (MET)
    "HS": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Thiol hydrogen
    "HO": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Hydroxyl hydrogen
    "HZ": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Histidine hydrogen
    "HA": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Aromatic hydrogen (general)
    "H4": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Aromatic hydrogen (HIS C4, TRP N1-H)
    "H5": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Aromatic hydrogen (HIS C2)
    # Carbon atoms
    "CT": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # sp3 carbon
    "C": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # Carbonyl carbon
    "CA": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # Aromatic carbon
    "CB": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (general)
    "CC": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (general)
    "CD": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # Alkene carbon
    "CE": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # Alkene carbon
    "CF": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # Alkyne carbon
    "CG": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # sp2 carbon (general)
    "CH": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # sp2 carbon (general)
    "CR": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (5-membered ring)
    "CW": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (5-membered ring)
    "CV": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (5-membered ring)
    "CZ": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # Phenyl carbon
    "C*": {"element": "C", "radius": 0.175, "screen": 0.75, "mass": 12.01},  # Aromatic carbon (sp2)
    "CN": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (TRP indole C2)
    # Oxygen atoms
    "O": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 16.00},  # Carbonyl oxygen
    "O2": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 16.00},  # Carboxyl oxygen
    "OH": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 16.00},  # Hydroxyl oxygen
    "OS": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 16.00},  # Ether oxygen
    "OW": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 16.00},  # Water oxygen
    # Nitrogen atoms
    "N": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Amide nitrogen
    "N3": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Protonated amine
    "NT": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Terminal amine
    "N2": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Amine nitrogen
    "NA": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Pyrrole nitrogen
    "NB": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Pyridine nitrogen
    "NC": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Cyano nitrogen
    "ND": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Pyrazine nitrogen
    "NE": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Imidazole nitrogen
    "NF": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.01,
    },  # Pyrrole nitrogen (general)
    "NG": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Indole nitrogen
    "NH": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.01,
    },  # Pyridine nitrogen (general)
    "NI": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Isoquinoline nitrogen
    "NL": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Aliphatic nitrogen
    "NM": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.01,
    },  # Amide nitrogen (general)
    "NP": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Phosphate nitrogen
    "NQ": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.01},  # Guanidinium nitrogen
    # Sulfur atoms
    "S": {"element": "S", "radius": 0.180, "screen": 0.85, "mass": 32.06},  # Thioether sulfur
    "SH": {"element": "S", "radius": 0.180, "screen": 0.85, "mass": 32.06},  # Thiol sulfur
    "S*": {"element": "S", "radius": 0.180, "screen": 0.85, "mass": 32.06},  # Sulfur (general)
    "SM": {"element": "S", "radius": 0.180, "screen": 0.85, "mass": 32.06},  # Sulfhydryl sulfur
    # Halogen atoms
    "F": {"element": "F", "radius": 0.150, "screen": 0.85, "mass": 19.00},  # Fluorine
    "CL": {"element": "Cl", "radius": 0.180, "screen": 0.85, "mass": 35.45},  # Chlorine
    "BR": {"element": "Br", "radius": 0.200, "screen": 0.85, "mass": 79.90},  # Bromine
    "I": {"element": "I", "radius": 0.220, "screen": 0.85, "mass": 126.90},  # Iodine
    # Phosphorus
    "P": {"element": "P", "radius": 0.185, "screen": 0.85, "mass": 30.97},  # Phosphate phosphorus
    # Unknown/default (fallback values)
    "X": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # Wildcard type
    # =========================================================================
    # Additional atom types for compatibility with other force fields
    # Added based on missing_gbsa_type_mappings.txt analysis
    # =========================================================================
    # a99SBdisp force field types
    "C1": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (a99SBdisp)
    "C3": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (a99SBdisp)
    "C4": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (a99SBdisp)
    "C5": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (a99SBdisp)
    "C6": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (a99SBdisp)
    "C7": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (a99SBdisp)
    "C8": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (a99SBdisp)
    "C9": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (a99SBdisp)
    "HB": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Backbone amide hydrogen (a99SBdisp)
    "O3": {
        "element": "O",
        "radius": 0.150,
        "screen": 0.85,
        "mass": 16.00,
    },  # Hydroxyl oxygen (a99SBdisp)
    "OB": {
        "element": "O",
        "radius": 0.150,
        "screen": 0.85,
        "mass": 16.00,
    },  # Carbonyl oxygen (a99SBdisp)
    # amber03wsc force field types
    "CAx": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (amber03wsc)
    "CTx": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (amber03wsc)
    "Cx": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Carbonyl carbon (amber03wsc)
    "H1x": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Aliphatic hydrogen (amber03wsc)
    "HCx": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Alkyl hydrogen (amber03wsc)
    "HPx": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Aromatic hydrogen (amber03wsc)
    "Hx": {"element": "H", "radius": 0.100, "screen": 0.85, "mass": 1.008},  # Hydrogen (amber03wsc)
    "N2x": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.01,
    },  # Amine nitrogen (amber03wsc)
    "N3x": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.01,
    },  # Protonated amine (amber03wsc)
    "Nx": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.01,
    },  # Amide nitrogen (amber03wsc)
    "O2x": {
        "element": "O",
        "radius": 0.150,
        "screen": 0.85,
        "mass": 16.00,
    },  # Carboxyl oxygen (amber03wsc)
    "Ox": {
        "element": "O",
        "radius": 0.150,
        "screen": 0.85,
        "mass": 16.00,
    },  # Carbonyl oxygen (amber03wsc)
    # amber14sb_parmbsc1 force field types
    "2C": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (amber14sb_parmbsc1)
    "3C": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (amber14sb_parmbsc1)
    "CO": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # Carbonyl carbon (amber14sb_parmbsc1)
    "CX": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (amber14sb_parmbsc1)
    # des-amber force field types (single-letter codes for residues)
    "AA": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (ALA in des-amber)
    "DD": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # sp2 carbon (ASP in des-amber)
    "EE": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # sp2 carbon (GLU in des-amber)
    "FF": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (PHE in des-amber)
    "GG": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (GLY in des-amber)
    "HE": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # sp2 carbon (HIS in des-amber)
    "II": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (ILE in des-amber)
    "KK": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # sp2 carbon (LYS in des-amber)
    "LL": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (LEU in des-amber)
    "MM": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (MET in des-amber)
    "NN": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # sp2 carbon (ASN in des-amber)
    "PP": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (PRO in des-amber)
    "QQ": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # sp2 carbon (GLN in des-amber)
    "RR": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # sp2 carbon (ARG in des-amber)
    "SS": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (SER in des-amber)
    "TT": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (THR in des-amber)
    "VV": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (VAL in des-amber)
    "WW": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (TRP in des-amber)
    "YY": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Aromatic carbon (TYR in des-amber)
    # des-amber special atom types
    "C&": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.75,
        "mass": 12.01,
    },  # Special carbon (des-amber)
    "CT_CT": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # sp3 carbon (des-amber)
    "H1_H1B": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Aliphatic hydrogen (des-amber)
    "H_HN": {
        "element": "H",
        "radius": 0.100,
        "screen": 0.85,
        "mass": 1.008,
    },  # Amide hydrogen (des-amber)
    "N_N": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.01,
    },  # Amide nitrogen (des-amber)
    "N_N2": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.01,
    },  # Amine nitrogen (des-amber)
    "OHS": {
        "element": "O",
        "radius": 0.150,
        "screen": 0.85,
        "mass": 16.00,
    },  # Hydroxyl oxygen (des-amber)
    "O_O": {
        "element": "O",
        "radius": 0.150,
        "screen": 0.85,
        "mass": 16.00,
    },  # Carbonyl oxygen (des-amber)
}

# CHARMM36 atom type mapping for GBSA implicit solvent calculations
CHARMM36_ATOM_MAPPING = {
    # Hydrogen atoms
    "H": {
        "element": "H",
        "radius": 0.105,
        "screen": 0.85,
        "mass": 1.008,
    },  # Backbone amide hydrogen
    "HA": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Alpha hydrogen
    "HB1": {
        "element": "H",
        "radius": 0.105,
        "screen": 0.85,
        "mass": 1.008,
    },  # Beta hydrogen (methyl)
    "HB2": {
        "element": "H",
        "radius": 0.105,
        "screen": 0.85,
        "mass": 1.008,
    },  # Beta hydrogen (methylene)
    "HB3": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Beta hydrogen
    "HC": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Alkyl hydrogen
    "HA1": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Alpha hydrogen (ILE)
    "HA2": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Alpha hydrogen
    "HA3": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Methyl hydrogen
    "HP": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Aromatic hydrogen
    "HR1": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # HIS ring hydrogen
    "HR2": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # HIS ring hydrogen
    "HR3": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # HIS ring hydrogen
    "HS": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Thiol hydrogen
    "HZ": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # HIS ring hydrogen
    "HZ1": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Lysine hydrogen
    "HZ2": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Lysine hydrogen
    "HZ3": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Lysine hydrogen
    "HD1": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # HIS ring hydrogen
    "HD2": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # HIS ring hydrogen
    "HE1": {
        "element": "H",
        "radius": 0.105,
        "screen": 0.85,
        "mass": 1.008,
    },  # HIS/PHE/TRP ring hydrogen
    "HE2": {
        "element": "H",
        "radius": 0.105,
        "screen": 0.85,
        "mass": 1.008,
    },  # HIS/PHE ring hydrogen
    "HE3": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # MET methyl hydrogen
    "HG1": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Sidechain hydrogen
    "HG2": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # Sidechain hydrogen
    "HG3": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # MET methyl hydrogen
    "HH": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # TYR hydroxyl hydrogen
    "HH1": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # ARG hydrogen
    "HH2": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # ARG hydrogen
    "HG21": {
        "element": "H",
        "radius": 0.105,
        "screen": 0.85,
        "mass": 1.008,
    },  # THR/VAL/ILE methyl hydrogen
    "HG22": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},
    "HG23": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},
    "HD11": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # LEU methyl hydrogen
    "HD12": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},
    "HD13": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},
    "HD21": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # ASN/GLN hydrogen
    "HD22": {"element": "H", "radius": 0.105, "screen": 0.85, "mass": 1.008},  # TRP hydrogen
    # Carbon atoms
    "CT1": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # Alpha carbon
    "CT2": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # Beta carbon
    "CT3": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # Methyl carbon
    "C": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # Carbonyl carbon
    "CA": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # Aromatic carbon
    "CB": {
        "element": "C",
        "radius": 0.170,
        "screen": 0.72,
        "mass": 12.01,
    },  # Beta carbon (aromatic)
    "CC": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.72,
        "mass": 12.01,
    },  # Carbonyl carbon (carboxyl)
    "CD": {
        "element": "C",
        "radius": 0.175,
        "screen": 0.72,
        "mass": 12.01,
    },  # Carbonyl carbon (carboxyl)
    "CE": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # Alkene carbon
    "CF": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # Aromatic carbon
    "CG": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # Aromatic carbon
    "CH": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # Aromatic carbon
    "CP1": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # PRO alpha carbon
    "CP2": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # PRO beta carbon
    "CP3": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # PRO delta carbon
    "CY": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # TRP carbon
    "CAI": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # TRP indole carbon
    "CPT": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # TRP carbon
    "CPH": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # HIS carbon
    "CPH1": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # HIS carbon
    "CPH2": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # HIS carbon
    "CS": {"element": "C", "radius": 0.175, "screen": 0.72, "mass": 12.01},  # Sulfonyl carbon
    "CT2A": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # ASP/GLU beta carbon
    "CT3A": {"element": "C", "radius": 0.170, "screen": 0.72, "mass": 12.01},  # Gamma carbon
    # Oxygen atoms
    "O": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 15.999},  # Carbonyl oxygen
    "O2": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 15.999},  # Carboxyl oxygen
    "OH1": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 15.999},  # Hydroxyl oxygen
    "OC": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 15.999},  # Carboxyl oxygen
    "ON": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 15.999},  # Nitro oxygen
    "OS": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 15.999},  # Ether oxygen
    "OT1": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 15.999},  # Terminal oxygen
    "OT2": {"element": "O", "radius": 0.150, "screen": 0.85, "mass": 15.999},  # Terminal oxygen
    # Nitrogen atoms
    "N": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.007},  # Amide nitrogen
    "NH1": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.007},  # Amide nitrogen
    "NH2": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.007,
    },  # Amide nitrogen (GLN/ASN)
    "NH3": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.007},  # Terminal amine
    "NC2": {
        "element": "N",
        "radius": 0.155,
        "screen": 0.85,
        "mass": 14.007,
    },  # Guanidinium nitrogen
    "NR1": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.007},  # HIS ring nitrogen
    "NR2": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.007},  # HIS ring nitrogen
    "NR3": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.007},  # Protonated amine
    "NY": {"element": "N", "radius": 0.155, "screen": 0.85, "mass": 14.007},  # TRP nitrogen
    # Sulfur atoms
    "S": {"element": "S", "radius": 0.180, "screen": 0.85, "mass": 32.06},  # Thioether sulfur
    "SH1": {"element": "S", "radius": 0.180, "screen": 0.85, "mass": 32.06},  # Thiol sulfur
    "SM": {"element": "S", "radius": 0.180, "screen": 0.85, "mass": 32.06},  # Sulfhydryl sulfur
    # Phosphorus
    "P": {"element": "P", "radius": 0.185, "screen": 0.85, "mass": 30.974},  # Phosphate phosphorus
    # Halogen atoms
    "F": {"element": "F", "radius": 0.150, "screen": 0.85, "mass": 18.998},  # Fluorine
    "CL": {"element": "Cl", "radius": 0.180, "screen": 0.85, "mass": 35.45},  # Chlorine
    "BR": {"element": "Br", "radius": 0.200, "screen": 0.85, "mass": 79.904},  # Bromine
    "I": {"element": "I", "radius": 0.220, "screen": 0.85, "mass": 126.90},  # Iodine
}


def _add_gbsa_solvent(
    sys, top, gb_model="GBn2", salt_conc=0.0, nonbondedCutoff=None, forcefield_type="AMBER"
):
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
    if forcefield_type.upper() not in ["AMBER", "CHARMM"]:
        raise ValueError(
            f"Unsupported forcefield_type: {forcefield_type}. Supported types: 'AMBER', 'CHARMM'"
        )

    # Select the appropriate atom type mapping based on forcefield
    if forcefield_type.upper() == "AMBER":
        atom_mapping = AMBER99SB_ILDN_ATOM_MAPPING
    else:  # CHARMM
        atom_mapping = CHARMM36_ATOM_MAPPING

    # Validate gb_model
    if gb_model not in ["OBC2", "GBn2"]:
        raise ValueError(f"Unsupported GB model: {gb_model}. Supported models: OBC2, GBn2")

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
    if gb_model == "OBC2":
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
        if kappa > 0 and hasattr(gb_force, "setKappa"):
            gb_force.setKappa(kappa)

    elif gb_model == "GBn2":
        # Use custom GBn2 force from OpenMM internal
        gb_force = GBSAGBn2Force(
            solventDielectric=solvent_dielectric,
            soluteDielectric=solute_dielectric,
            SA=None,  # No surface area contribution
            cutoff=cutoff,
            kappa=kappa,
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
    if gb_model == "OBC2":
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
                radius = atom_info["radius"]
                screening = atom_info.get("screen", 0.8)
            elif atom.element and atom.element.symbol in ["H", "C", "N", "O", "S", "P"]:
                # Fallback to element-based defaults
                element_defaults = {
                    "H": (0.100, 0.85),
                    "C": (0.170, 0.72),
                    "N": (0.155, 0.79),
                    "O": (0.150, 0.85),
                    "S": (0.180, 0.96),
                    "P": (0.185, 0.86),
                }
                defaults = element_defaults.get(atom.element.symbol, (0.170, 0.8))
                radius, screening = defaults
            else:
                radius, screening = 0.170, 0.8

            gb_force.addParticle(charge, radius, screening)

    elif gb_model == "GBn2":
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
    if hasattr(gb_force, "finalize") and gb_model == "GBn2":
        gb_force.finalize()

    sys.addForce(gb_force)
