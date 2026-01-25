import math
import os
import re
import shutil
from collections import OrderedDict, defaultdict

import simtk.openmm as mm
import simtk.unit as unit
from simtk.openmm.app import PDBFile, Topology
from simtk.openmm.app import amberprmtopfile as prmtop
from simtk.openmm.app import element as elem
from simtk.openmm.app import forcefield as ff


HBonds = ff.HBonds
AllBonds = ff.AllBonds
HAngles = ff.HAngles

# NO NEEDED as we do not use implicit solvent
#OBC2 = prmtop.OBC2


novarcharre = re.compile(r"\W")

# Ring constraints data for aromatic residues (rigid triangle constraints)
# Force constants:
# - Distance constraints (12, 13, 14): 2500000 kJ/mol/nm^2
# - Angle constraints (13): 5000 kJ/mol/rad^2
# Distances in nm, angles in degrees
RING_CONSTRAINTS = {
    "PHE": {
        "constraints_12": [
            {"atom1": "CG", "atom2": "CD1", "distance_nm": 0.138885},
            {"atom1": "CG", "atom2": "CD2", "distance_nm": 0.138517},
            {"atom1": "CD1", "atom2": "CE1", "distance_nm": 0.138293},
            {"atom1": "CE1", "atom2": "CZ", "distance_nm": 0.13886},
            {"atom1": "CZ", "atom2": "CE2", "distance_nm": 0.139295},
            {"atom1": "CE2", "atom2": "CD2", "distance_nm": 0.139449},
        ],
        "constraints_13": [
            {"atom1": "CG", "atom2": "CE1", "distance_nm": 0.239904},
            {"atom1": "CG", "atom2": "CE2", "distance_nm": 0.240817},
            {"atom1": "CD1", "atom2": "CZ", "distance_nm": 0.240385},
            {"atom1": "CD1", "atom2": "CD2", "distance_nm": 0.240479},
            {"atom1": "CE1", "atom2": "CE2", "distance_nm": 0.240875},
            {"atom1": "CZ", "atom2": "CD2", "distance_nm": 0.240842},
        ],
        "constraints_14": [
            {"atom1": "CG", "atom2": "CZ", "distance_nm": 0.277557},
            {"atom1": "CD1", "atom2": "CE2", "distance_nm": 0.278262},
            {"atom1": "CE1", "atom2": "CD2", "distance_nm": 0.277469},
        ],
        "angles_13": [
            {"atom1": "CG", "atom2": "CD1", "atom3": "CE1", "angle_deg": 119.8849},
            {"atom1": "CG", "atom2": "CD2", "atom3": "CE2", "angle_deg": 120.0752},
            {"atom1": "CD1", "atom2": "CE1", "atom3": "CZ", "angle_deg": 120.3013},
            {"atom1": "CD1", "atom2": "CG", "atom3": "CD2", "angle_deg": 120.1999},
            {"atom1": "CE1", "atom2": "CZ", "atom3": "CE2", "angle_deg": 119.989},
            {"atom1": "CZ", "atom2": "CE2", "atom3": "CD2", "angle_deg": 119.5439},
        ]
    },
    "TYR": {
        "constraints_12": [
            {"atom1": "CG", "atom2": "CD1", "distance_nm": 0.138401},
            {"atom1": "CG", "atom2": "CD2", "distance_nm": 0.138874},
            {"atom1": "CD1", "atom2": "CE1", "distance_nm": 0.138196},
            {"atom1": "CE1", "atom2": "CZ", "distance_nm": 0.138788},
            {"atom1": "CZ", "atom2": "CE2", "distance_nm": 0.139176},
            {"atom1": "CE2", "atom2": "CD2", "distance_nm": 0.138755},
        ],
        "constraints_13": [
            {"atom1": "CG", "atom2": "CE1", "distance_nm": 0.239785},
            {"atom1": "CG", "atom2": "CE2", "distance_nm": 0.240443},
            {"atom1": "CD1", "atom2": "CZ", "distance_nm": 0.239716},
            {"atom1": "CD1", "atom2": "CD2", "distance_nm": 0.240235},
            {"atom1": "CE1", "atom2": "CE2", "distance_nm": 0.240986},
            {"atom1": "CZ", "atom2": "CD2", "distance_nm": 0.240202},
        ],
        "constraints_14": [
            {"atom1": "CG", "atom2": "CZ", "distance_nm": 0.276946},
            {"atom1": "CD1", "atom2": "CE2", "distance_nm": 0.277658},
            {"atom1": "CE1", "atom2": "CD2", "distance_nm": 0.277563},
        ],
        "angles_13": [
            {"atom1": "CG", "atom2": "CD1", "atom3": "CE1", "angle_deg": 120.2035},
            {"atom1": "CG", "atom2": "CD2", "atom3": "CE2", "angle_deg": 120.0078},
            {"atom1": "CD1", "atom2": "CE1", "atom3": "CZ", "angle_deg": 119.8692},
            {"atom1": "CD1", "atom2": "CG", "atom3": "CD2", "angle_deg": 120.0888},
            {"atom1": "CE1", "atom2": "CZ", "atom3": "CE2", "angle_deg": 120.2163},
            {"atom1": "CZ", "atom2": "CE2", "atom3": "CD2", "angle_deg": 119.5941},
        ]
    },
    "TRP": {
        "six_membered": {
            "constraints_12": [
                {"atom1": "CD2", "atom2": "CE2", "distance_nm": 0.138788},
                {"atom1": "CD2", "atom2": "CE3", "distance_nm": 0.138293},
                {"atom1": "CE2", "atom2": "CZ2", "distance_nm": 0.138802},
                {"atom1": "CZ2", "atom2": "CH2", "distance_nm": 0.13877},
                {"atom1": "CH2", "atom2": "CZ3", "distance_nm": 0.139341},
                {"atom1": "CZ3", "atom2": "CE3", "distance_nm": 0.138802},
            ],
            "constraints_13": [
                {"atom1": "CD2", "atom2": "CZ2", "distance_nm": 0.240242},
                {"atom1": "CD2", "atom2": "CZ3", "distance_nm": 0.239627},
                {"atom1": "CE2", "atom2": "CH2", "distance_nm": 0.240501},
                {"atom1": "CE2", "atom2": "CE3", "distance_nm": 0.240452},
                {"atom1": "CZ2", "atom2": "CZ3", "distance_nm": 0.240452},
                {"atom1": "CH2", "atom2": "CE3", "distance_nm": 0.24116},
            ],
            "constraints_14": [
                {"atom1": "CD2", "atom2": "CH2", "distance_nm": 0.277512},
                {"atom1": "CE2", "atom2": "CZ3", "distance_nm": 0.277353},
                {"atom1": "CZ2", "atom2": "CE3", "distance_nm": 0.277923},
            ],
            "angles_13": [
                {"atom1": "CD2", "atom2": "CE2", "atom3": "CZ2", "angle_deg": 119.8697},
                {"atom1": "CD2", "atom2": "CE3", "atom3": "CZ3", "angle_deg": 119.7153},
                {"atom1": "CE2", "atom2": "CZ2", "atom3": "CH2", "angle_deg": 120.097},
                {"atom1": "CE2", "atom2": "CD2", "atom3": "CE3", "angle_deg": 120.4087},
                {"atom1": "CZ2", "atom2": "CH2", "atom3": "CZ3", "angle_deg": 119.6715},
                {"atom1": "CH2", "atom2": "CZ3", "atom3": "CE3", "angle_deg": 120.2316},
            ]
        },
        "five_membered": {
            "constraints_12": [
                {"atom1": "CG", "atom2": "CD1", "distance_nm": 0.135092},
                {"atom1": "CG", "atom2": "CD2", "distance_nm": 0.14521},
                {"atom1": "CD1", "atom2": "NE1", "distance_nm": 0.13731},
                {"atom1": "NE1", "atom2": "CE2", "distance_nm": 0.136675},
                {"atom1": "CE2", "atom2": "CD2", "distance_nm": 0.138788},
            ],
            "constraints_13": [
                {"atom1": "CG", "atom2": "NE1", "distance_nm": 0.223589},
                {"atom1": "CG", "atom2": "CE2", "distance_nm": 0.226363},
                {"atom1": "CD1", "atom2": "CE2", "distance_nm": 0.221572},
                {"atom1": "CD1", "atom2": "CD2", "distance_nm": 0.224958},
                {"atom1": "NE1", "atom2": "CD2", "distance_nm": 0.224753},
            ],
            "constraints_14": [],
            "angles_13": [
                {"atom1": "CG", "atom2": "CD1", "atom3": "NE1", "angle_deg": 110.328},
                {"atom1": "CG", "atom2": "CD2", "atom3": "CE2", "angle_deg": 105.6778},
                {"atom1": "CD1", "atom2": "NE1", "atom3": "CE2", "angle_deg": 107.9382},
                {"atom1": "CD1", "atom2": "CG", "atom3": "CD2", "angle_deg": 106.6936},
                {"atom1": "NE1", "atom2": "CE2", "atom3": "CD2", "angle_deg": 109.3529},
            ]
        }
    },
    "HIS": {
        "constraints_12": [
            {"atom1": "CG", "atom2": "CD2", "distance_nm": 0.135207},
            {"atom1": "CG", "atom2": "ND1", "distance_nm": 0.138148},
            {"atom1": "CD2", "atom2": "NE2", "distance_nm": 0.136697},
            {"atom1": "NE2", "atom2": "CE1", "distance_nm": 0.129008},
            {"atom1": "CE1", "atom2": "ND1", "distance_nm": 0.134503},
        ],
        "constraints_13": [
            {"atom1": "CG", "atom2": "NE2", "distance_nm": 0.22309},
            {"atom1": "CG", "atom2": "CE1", "distance_nm": 0.218874},
            {"atom1": "CD2", "atom2": "CE1", "distance_nm": 0.212106},
            {"atom1": "CD2", "atom2": "ND1", "distance_nm": 0.216947},
            {"atom1": "NE2", "atom2": "ND1", "distance_nm": 0.218486},
        ],
        "constraints_14": [],
        "angles_13": [
            {"atom1": "CG", "atom2": "CD2", "atom3": "NE2", "angle_deg": 110.263},
            {"atom1": "CG", "atom2": "ND1", "atom3": "CE1", "angle_deg": 106.7819},
            {"atom1": "CD2", "atom2": "NE2", "atom3": "CE1", "angle_deg": 105.8958},
            {"atom1": "CD2", "atom2": "CG", "atom3": "ND1", "angle_deg": 105.0496},
            {"atom1": "NE2", "atom2": "CE1", "atom3": "ND1", "angle_deg": 112.003},
        ]
    }
}

def _find_all_instances_in_string(string, substr):
    """Find indices of all instances of substr in string"""
    indices = []
    idx = string.find(substr, 0)
    while idx > -1:
        indices.append(idx)
        idx = string.find(substr, idx + 1)
    return indices

def _replace_defines(line, defines):
    """Replaces defined tokens in a given line"""
    if not defines:
        return line
    for define in reversed(defines):
        value = defines[define]
        indices = _find_all_instances_in_string(line, define)
        if not indices:
            continue
        # Check to see if it's inside of quotes
        inside = ""
        idx = 0
        n_to_skip = 0
        new_line = []
        for i, char in enumerate(line):
            if n_to_skip:
                n_to_skip -= 1
                continue
            if char in ("'\""):
                if not inside:
                    inside = char
                else:
                    if inside == char:
                        inside = ""
            if idx < len(indices) and i == indices[idx]:
                if inside:
                    new_line.append(char)
                    idx += 1
                    continue
                if i == 0 or novarcharre.match(line[i - 1]):
                    endidx = indices[idx] + len(define)
                    if endidx >= len(line) or novarcharre.match(line[endidx]):
                        new_line.extend(list(value))
                        n_to_skip = len(define) - 1
                        idx += 1
                        continue
                idx += 1
            new_line.append(char)
        line = "".join(new_line)

    return line

class PACETopFile(object):
    """Parse a Gromacs Martini top file and constructs a Topology and (optionally) an OpenMM System from it."""

    class _MoleculeType(object):
        """Inner class to store information about a molecule type."""

        def __init__(self):
            self.molecule_name = ""
            self.atoms = []
            self.bonds = []
            self.harmonic_angles = []
            self.dihedrals = []
            self.exclusions = []
            self.pairs = []
            self.constraints = []
            self.cmaps = []

    def __init__(
        self,
        file,
        periodicBoxVectors=None,
        unitCellDimensions=None,
        includeDir=None,
        defines=None,
        epsilon_r=15.0,
    ):
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
        """
        self.epsilon_r = epsilon_r
        if includeDir is None:
            includeDir = _get_default_gromacs_include_dir()
        self._includeDirs = (os.path.dirname(file), includeDir)
        self._defines = OrderedDict()
        self._genpairs = True
        if defines is not None:
            for define, value in defines.items():
                self._defines[define] = value

        self._currentCategory = None
        self._ifStack = []
        self._elseStack = []
        self._moleculeTypes = {}
        self._molecules = []
        self._currentMoleculeType = None
        self._atom_types = {}
        self._bondTypes = {}
        self._angleTypes = {}
        self._dihedralTypes = {}
        self._implicitTypes = {}
        self._pairTypes = {}
        self._cmapTypes = {}
        self._nonbond_types = {}
        #self._all_vsites = []
        self.nb_force = None
        self.es_self_excl_force = None
        self.es_except_force = None
        self.lj_except_force = None
        self.harmonic_angle_force = None
        #self.g96_angle_force = None
        #self.restricted_angle_force = None
        self.harmonic_bond_force = None
        self.periodic_torsion_force = None
        self.harmonic_torsion_force = None
        #self.combined_bending_torsion_force = None
        #self.rb_torsion_force = None
        self._use_sigma_eps = False
        self._use_g96_angles = False
        self._use_harmonic_angles = False
        self._use_restricted_angles = False
        
        self._process_file(file)

        top = Topology()
        self.topology = top
        if periodicBoxVectors is not None:
            if unitCellDimensions is not None:
                raise ValueError(
                    "specify either periodicBoxVectors or unitCellDimensions, but not both"
                )
            top.setPeriodicBoxVectors(periodicBoxVectors)
        else:
            top.setUnitCellDimensions(unitCellDimensions)
        PDBFile._loadNameReplacementTables()            


        for molecule_name, molecule_count in self._molecules:
            if molecule_name not in self._moleculeTypes:
                raise ValueError("Unknown molecule type: " + molecule_name)
            molecule_type = self._moleculeTypes[molecule_name]
            # Create the specified number of molecules of this type.
            for i in range(molecule_count):
                atoms = []
                lastResidue = None
                c = top.addChain()
                for index, fields in enumerate(molecule_type.atoms):
                    resNumber = fields[2]
                    if resNumber != lastResidue:
                        lastResidue = resNumber
                        resName = fields[3]
                        if resName in PDBFile._residueNameReplacements:
                            resName = PDBFile._residueNameReplacements[resName]
                        r = top.addResidue(resName, c)
                        if resName in PDBFile._atomNameReplacements:
                            atomReplacements = PDBFile._atomNameReplacements[resName]
                        else:
                            atomReplacements = {}
                    atomName = fields[4]
                    if atomName in atomReplacements:
                        atomName = atomReplacements[atomName]

                    # Try to guess the element.
                    upper = atomName.upper()
                    if upper.startswith("CL"):
                        element = elem.chlorine
                    elif upper.startswith("NA"):
                        element = elem.sodium
                    #elif upper.startswith("MG"):
                    #    element = elem.magnesium
                    else:
                        try:
                            element = elem.get_by_symbol(atomName[0])
                        except KeyError:
                            element = None
                    atoms.append(top.addAtom(atomName, element, r))

                # Add bonds to the topology

                for fields in molecule_type.bonds:
                    top.addBond(atoms[int(fields[0]) - 1], atoms[int(fields[1]) - 1])


    def _process_file(self, file):
        append = ""
        with open(file) as lines:
            for line in lines:
                if line.strip().endswith("\\"):
                    append = "%s %s" % (append, line[: line.rfind("\\")])
                else:
                    self._process_line(append + " " + line, file)
                    append = ""


    def _process_line(self, line, file):
        """Process one line from a file."""
        if ";" in line:
            line = line[: line.index(";")]
        stripped = line.strip()
        ignore = not all(self._ifStack)
        if stripped.startswith("*") or len(stripped) == 0:
            # A comment or empty line.
            return

        elif stripped.startswith("[") and not ignore:
            # The start of a category.
            if not stripped.endswith("]"):
                raise ValueError("Illegal line in .top file: " + line)
            self._currentCategory = stripped[1:-1].strip()

        elif stripped.startswith("#"):
            # A preprocessor command.
            fields = stripped.split()
            command = fields[0]
            if len(self._ifStack) != len(self._elseStack):
                raise RuntimeError("#if/#else stack out of sync")

            if command == "#include" and not ignore:
                # Locate the file to include
                name = stripped[len(command) :].strip(' \t"<>')
                searchDirs = self._includeDirs + (os.path.dirname(file),)
                for dir in searchDirs:
                    file = os.path.join(dir, name)
                    if os.path.isfile(file):
                        # We found the file, so process it.
                        self._process_file(file)
                        break
                else:
                    raise ValueError("Could not locate #include file: " + name)

            elif command == "#define" and not ignore:
                # Add a value to our list of defines.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                valueStart = stripped.find(name, len(command)) + len(name) + 1
                value = line[valueStart:].strip()
                value = value or "1"  # Default define is 1
                self._defines[name] = value
    

            elif command == "#ifdef":
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                self._ifStack.append(name in self._defines)
                self._elseStack.append(False)
        
            elif command == "#undef":
                # Un-define a variable
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                if fields[1] in self._defines:
                    self._defines.pop(fields[1])
        
            elif command == "#ifndef":
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                self._ifStack.append(name not in self._defines)
                self._elseStack.append(False)
    
            elif command == "#endif":
                # Pop an entry off the if stack.
                if len(self._ifStack) == 0:
                    raise ValueError("Unexpected line in .top file: " + line)
                del self._ifStack[-1]
                del self._elseStack[-1]

            elif command == "#else":
                # Reverse the last entry on the if stack
                if len(self._ifStack) == 0:
                    raise ValueError("Unexpected line in .top file: " + line)
                if self._elseStack[-1]:
                    raise ValueError(
                        "Unexpected line in .top file: "
                        "#else has already been used " + line
                    )
                self._ifStack[-1] = not self._ifStack[-1]
                self._elseStack[-1] = True

        elif not ignore:
            # Gromacs occasionally uses #define's to introduce specific
            # parameters for individual terms (for instance, this is how
            # ff99SB-ILDN is implemented). So make sure we do the appropriate
            # pre-processor replacements necessary
            line = _replace_defines(line, self._defines)
            # A line of data for the current category
            if self._currentCategory is None:
                raise ValueError("Unexpected line in .top file:  " + line)
            if self._currentCategory == "defaults":
                self._process_defaults(line)
            elif self._currentCategory == "moleculetype":
                self._process_molecule_type(line)
            elif self._currentCategory == "molecules":
                self._process_molecule(line)
            elif self._currentCategory == "atoms":
                self._process_atom(line)
            elif self._currentCategory == "bonds":
                self._process_bond(line)
            elif self._currentCategory == "angles":
                self._process_angle(line)
            elif self._currentCategory == "dihedrals":
                self._process_dihedral(line)
            elif self._currentCategory == "exclusions":
                self._process_exclusion(line)
            elif self._currentCategory == "pairs":
                self._process_pair(line)
            elif self._currentCategory == "constraints":
                self._process_constraint(line)
            elif self._currentCategory == "cmap":
                self._process_cmap(line)
            elif self._currentCategory == "atomtypes":
                self._process_atom_type(line)
            elif self._currentCategory == "bondtypes":
                self._process_bond_type(line)
            elif self._currentCategory == "angletypes":
                self._process_angle_type(line)
            elif self._currentCategory == "dihedraltypes":
                self._process_dihedral_type(line)
            elif self._currentCategory == "pairtypes":
                self._process_pair_type(line)
            elif self._currentCategory == "cmaptypes":
                self._process_cmap_type(line)
            elif self._currentCategory == "nonbond_params":
                self._process_nonbond_type(line)
            elif self._currentCategory == "nonbond-params":
                self._process_nonbond_type(line)

            
            elif self._currentCategory == "system":  # ignore the system description
                pass
            elif self._currentCategory == "position_restraints":
                raise ValueError(
                    "[ position_restraints ] is not currently supported.\n"
                    "The same effect can be achieved using a CustomExternalForce in OpenMM.\n"
                    "Please see github.com/maccallumlab/martini_openmm for details."
                )
            else:
                raise ValueError(
                    f"[ {self._currentCategory} ] is not currently supported.\n"
                    "Please file an issue at github.com/maccallumlab/martini_openmm"
                )

    def _process_defaults(self, line):
        """Process the [ defaults ] line."""
        fields = line.split()
        if len(fields) > 5:
            raise ValueError("Too many fields in [ defaults ] line: " + line)
        if fields[0] != "1":
            raise ValueError("Unsupported nonbonded type: " + fields[0])
        if fields[1] == "1":
            self._use_sigma_eps = False
        elif fields[1] == "2":
            self._use_sigma_eps = True
        else:
            raise ValueError("Unsupported combination rule: " + fields[1])
        self._defaults = fields            


    def _process_molecule_type(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            raise ValueError("Too few fields in [ moleculetypes ] line: " + line)
        type = PACETopFile._MoleculeType()
        type.molecule_name = fields[0]
        self._moleculeTypes[fields[0]] = type
        self._currentMoleculeType = type
    
    def _process_molecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            raise ValueError("Too few fields in [ molecules ] line: " + line)
        self._molecules.append((fields[0], int(fields[1])))

    def _process_atom(self, line):
        """Process a line in the [ atoms ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ atoms ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ atoms ] line: " + line)
        self._currentMoleculeType.atoms.append(fields)

    def _process_bond(self, line):
        """Process a line in the [ bonds ] category."""
        """For PACE we only have type 1 bond"""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ bonds ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 3:
            raise ValueError("Too few fields in [ bonds ] line: " + line)
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ bonds ] line: ' + line)

        #if fields[2] != "1" and fields[2] != "6":
        #    raise ValueError("Unsupported function type in [ bonds ] line: " + line)
        self._currentMoleculeType.bonds.append(fields)


    def _process_angle(self, line):
        """Process a line in the [ angles ] category."""
        '''For PACE we only have type 1 angle'''
        if self._currentMoleculeType is None:
            raise ValueError("Found [ angles ] section before [ moleculetype ]")
        fields = line.split()

        if len(fields) < 4:
            raise ValueError("Too few fields in [ angles ] line: " + line)
        if fields[3] == "1":
            self._use_harmonic_angles = True
            self._currentMoleculeType.harmonic_angles.append(fields)        
        else:
            raise ValueError(
                "Unsupported (nonG96) function type in [ angles ] line: " + line
            )

    def _process_dihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        """In PACE we have type 1 and type 2 angle"""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ dihedrals ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ dihedrals ] line: " + line)
        if fields[4] not in ("1", "2"):
            raise ValueError("Unsupported function type in [ dihedrals ] line: " + line)
        self._currentMoleculeType.dihedrals.append(fields)

    def _process_exclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ exclusions ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 2:
            raise ValueError("Too few fields in [ exclusions ] line: " + line)
        self._currentMoleculeType.exclusions.append(fields)

    def _process_pair(self, line):
        """Process a line in the [ pairs ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ pairs ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 3:
            raise ValueError("Too few fields in [ pairs ] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ pairs ] line: " + line)
        self._currentMoleculeType.pairs.append(fields)


    def _process_cmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ cmap ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 6:
            raise ValueError("Too few fields in [ cmap ] line: " + line)
        self._currentMoleculeType.cmaps.append(fields)

    def _process_atom_type(self, line):
        """Process a line in the [ atomtypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError("Too few fields in [ atomtypes ] line: " + line)
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
        self._atom_types[fields[0]] = fields

    def _process_bond_type(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ bondtypes ] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ bondtypes ] line: " + line)
        self._bondTypes[tuple(fields[:2])] = fields

    def _process_angle_type(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError("Too few fields in [ angletypes ] line: " + line)
        if fields[3] not in ("1", "5"):
            raise ValueError(
                "Unsupported function type in [ angletypes ] line: " + line
            )
        self._angleTypes[tuple(fields[:3])] = fields

    def _process_dihedral_type(self, line):
        """Process a line in the [ dihedraltypes ] category."""
        fields = line.split()
        
        if len(fields) < 7 and fields[2] != '1':
            raise ValueError("Too few fields in [ dihedraltypes ] line: " + line)
        #add X for PACE dihedraltypes
        if fields[2] == '1':
            # modify PACE term 
            # CBX SH  1  0.0  3.9  3  to  
            # X CBX SH X 1 0.0 3.9 3
            fields.insert(0, 'X')
            fields.insert(3, 'X')
        if fields[4] not in ("1", "2"):
            raise ValueError(
                "Unsupported function type in [ dihedraltypes ] line: " + line
            )
        key = tuple(fields[:5])
        self._dihedralTypes[key] = [fields]


    def _process_pair_type(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ pairtypes] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ pairtypes ] line: " + line)
        self._pairTypes[tuple(fields[:2])] = fields

    def _process_cmap_type(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8 + int(fields[6]) * int(fields[7]):
            raise ValueError("Too few fields in [ cmaptypes ] line: " + line)
        if fields[5] != "1":
            raise ValueError("Unsupported function type in [ cmaptypes ] line: " + line)
        self._cmapTypes[tuple(fields[:5])] = fields

    def _process_nonbond_type(self, line):
        """Process a line in the [ nonbond_params ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ nonbond_params ] line: " + line)
        if fields[2] != "1":
            raise ValueError(
                "Unsupported function type in [ nonbond_params ] line: " + line
            )
        self._nonbond_types[tuple(sorted(fields[:2]))] = fields

    def _build_dihedral_lookup_table(self):
        dihedral_type_table = {}
        # a dictionary can be iterated over in Python using a for loop.
        # The for loop will iterate over the keys of the dictionary by default. 
        # for PACE we do not have * X X * wildcard types but X * * X wildcard types.
        
        for key in self._dihedralTypes:
            if key[1] != "X" and key[2] != "X":
                if (key[1], key[2]) not in dihedral_type_table:
                    dihedral_type_table[(key[1], key[2])] = []
                dihedral_type_table[(key[1], key[2])].append(key)
                if (key[2], key[1]) not in dihedral_type_table:
                    dihedral_type_table[(key[2], key[1])] = []
                dihedral_type_table[(key[2], key[1])].append(key)
        wildcard_dihedral_types = []
        
        for key in self._dihedralTypes:
            if key[1] == "X" or key[2] == "X":
                wildcard_dihedral_types.append(key)
                for types in dihedral_type_table.values():
                    types.append(key)
            
        return dihedral_type_table, wildcard_dihedral_types

    def _add_atoms_to_system(self, sys, molecule_type):
        for fields in molecule_type.atoms:
            if len(fields) >= 8:
                mass = float(fields[7])
            else:
                mass = float(self._atom_types[fields[1]][3])
            sys.addParticle(mass)

    def _add_bonds_to_system(self, sys, molecule_type, bonded_types, base_atom_index):
        for fields in molecule_type.bonds:
            atoms = [int(x) - 1 for x in fields[:2]]
            types = tuple(bonded_types[i] for i in atoms)
            if len(fields) >= 5:
                params = fields[3:5]
            elif types in self._bondTypes:
                params = self._bondTypes[types][3:5]
            elif types[::-1] in self._bondTypes:
                params = self._bondTypes[types[::-1]][3:5]
            else:
                raise ValueError(
                    "No parameters specified for bond: " + fields[0] + ", " + fields[1]
                )

            length = float(params[0])

            if self.harmonic_bond_force is None:
                self.harmonic_bond_force = mm.HarmonicBondForce()
                sys.addForce(self.harmonic_bond_force)
            self.harmonic_bond_force.addBond(
                base_atom_index + atoms[0],
                base_atom_index + atoms[1],
                length,
                float(params[1]),
            )

    def _add_angles_to_system(self, sys, molecule_type, bonded_types, base_atom_index):
        #PACE Do not have g96 angle or restricted angles.
        self._add_angle_forces(sys)
        self._add_harmonic_angles_to_system(
            molecule_type, bonded_types, base_atom_index
        )


    def _add_angle_forces(self, sys):
        if self.harmonic_angle_force is None:
            self.harmonic_angle_force = mm.HarmonicAngleForce()
            sys.addForce(self.harmonic_angle_force)


    def _process_angle_params(self, fields, bonded_types):
        degToRad = math.pi / 180
        atoms = [int(x) - 1 for x in fields[:3]]
        types = tuple(bonded_types[i] for i in atoms)
        if len(fields) >= 6:
            params = fields[4:]
        elif types in self._angleTypes:
            params = self._angleTypes[types][4:]
        # reverse
        elif types[::-1] in self._angleTypes:
            params = self._angleTypes[types[::-1]][4:]
        else:
            raise ValueError(
                "No parameters specified for angle: "
                + fields[0]
                + ", "
                + fields[1]
                + ", "
                + fields[2]
            )

        theta = float(params[0]) * degToRad
        return atoms, theta, float(params[1])

    def _add_harmonic_angles_to_system(
        self, molecule_type, bonded_types, base_atom_index
    ):
        for fields in molecule_type.harmonic_angles:
            assert fields[3] == "1"
            atoms, theta, k = self._process_angle_params(fields, bonded_types)
            self.harmonic_angle_force.addAngle(
                base_atom_index + atoms[0],
                base_atom_index + atoms[1],
                base_atom_index + atoms[2],
                theta,
                k,
            )



    def _add_torsions_to_system(
        self,
        sys,
        molecule_type,
        bonded_types,
        dihedral_type_table,
        wildcard_dihedral_types,
        base_atom_index,
    ):
        degToRad = math.pi / 180

        for fields in molecule_type.dihedrals:
            atoms = [int(x) - 1 for x in fields[:4]]
            types = tuple(bonded_types[i] for i in atoms)
            dihedralType = fields[4]
            reversedTypes = types[::-1] + (dihedralType,)
            types = types + (dihedralType,)

            if (
                (dihedralType == '1' and len(fields) > 7) or (dihedralType == '2' and len(fields) > 6)
            ):
                paramsList = [fields]
            else:
                # Look for a matching dihedral type.
                # It seems PACE wildcard types can be treated in dihedral_type_table.
                paramsList = None
                if (types[1], types[2]) in dihedral_type_table:
                    dihedralTypes = dihedral_type_table[(types[1], types[2])]
                else:
                    dihedralTypes = wildcard_dihedral_types
                for key in dihedralTypes:
                    if all(a == b or a == "X" for a, b in zip(key, types)) or all(
                        a == b or a == "X" for a, b in zip(key, reversedTypes)
                    ):
                        paramsList = self._dihedralTypes[key]
                        if "X" not in key:
                            break
                if paramsList is None:
                    raise ValueError(
                        "No parameters specified for dihedral: "
                        + fields[0]
                        + ", "
                        + fields[1]
                        + ", "
                        + fields[2]
                        + ", "
                        + fields[3]
                    )
        
            for params in paramsList:
                if dihedralType in ("1", "4", "9"):
                    # Periodic torsion
                    k = float(params[6])
                    if k != 0:
                        if self.periodic_torsion_force is None:
                            self.periodic_torsion_force = mm.PeriodicTorsionForce()
                            sys.addForce(self.periodic_torsion_force)
                        self.periodic_torsion_force.addTorsion(
                            base_atom_index + atoms[0],
                            base_atom_index + atoms[1],
                            base_atom_index + atoms[2],
                            base_atom_index + atoms[3],
                            int(float(params[7])),
                            float(params[5]) * degToRad,
                            k,
                        )

                elif dihedralType == "2":
                    # Harmonic torsion
                    k = float(params[6])
                    phi0 = float(params[5])
                    if k != 0:
                        if self.harmonic_torsion_force is None:
                            self.harmonic_torsion_force = mm.CustomTorsionForce(
                                "0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g"
                                % math.pi
                            )
                            self.harmonic_torsion_force.addPerTorsionParameter("theta0")
                            self.harmonic_torsion_force.addPerTorsionParameter("k")
                            sys.addForce(self.harmonic_torsion_force)
                        # map phi0 into correct space
                            if phi0 > 180:
                                phi0 = phi0 - 360
                            else:
                                phi0 = phi0
                        self.harmonic_torsion_force.addTorsion(
                            base_atom_index + atoms[0],
                            base_atom_index + atoms[1],
                            base_atom_index + atoms[2],
                            base_atom_index + atoms[3],
                            (phi0 * degToRad, k),
                        )

    
    def _add_cmap_to_system(self, sys, molecule_type, bonded_types, base_atom_index):
        cmap = None
        mapIndices = {}
        for fields in molecule_type.cmaps:
            atoms = [int(x) - 1 for x in fields[:5]]
            types = tuple(bonded_types[i] for i in atoms)
            if len(fields) >= 8 and len(fields) >= 8 + int(fields[6]) * int(fields[7]):
                params = fields
            elif types in self._cmapTypes:
                params = self._cmapTypes[types]
            elif types[::-1] in self._cmapTypes:
                params = self._cmapTypes[types[::-1]]
            else:
                raise ValueError(
                    "No parameters specified for cmap: "
                    + fields[0]
                    + ", "
                    + fields[1]
                    + ", "
                    + fields[2]
                    + ", "
                    + fields[3]
                    + ", "
                    + fields[4]
                )
            if cmap is None:
                cmap = mm.CMAPTorsionForce()
                sys.addForce(cmap)
            mapSize = int(params[6])
            if mapSize != int(params[7]):
                raise ValueError("Non-square CMAPs are not supported")
            map = []
            for i in range(mapSize):
                for j in range(mapSize):
                    map.append(
                        float(
                            params[
                                8
                                + mapSize * ((j + mapSize // 2) % mapSize)
                                + ((i + mapSize // 2) % mapSize)
                            ]
                        )
                    )
            map = tuple(map)
            if map not in mapIndices:
                mapIndices[map] = cmap.addMap(mapSize, map)
            cmap.addTorsion(
                mapIndices[map],
                base_atom_index + atoms[0],
                base_atom_index + atoms[1],
                base_atom_index + atoms[2],
                base_atom_index + atoms[3],
                base_atom_index + atoms[1],
                base_atom_index + atoms[2],
                base_atom_index + atoms[3],
                base_atom_index + atoms[4],
            )

    def _add_ring_constraints(self, sys, molecule_type, base_atom_index):
        """Add rigid triangle constraints for aromatic rings (PHE, TYR, TRP, HIS).
        
        Parameters
        ----------
        sys : mm.System
            The OpenMM system to add constraints to
        molecule_type : _MoleculeType
            The molecule type containing atom information
        base_atom_index : int
            The base atom index for this molecule in the system
        """
        # Force constants:
        # - Distance constraints (12, 13, 14): 125000 kJ/mol/nm^2
        # - Angle constraints (13): 1000 kJ/mol/rad^2
        RING_FORCE_CONSTANT = 125000.0  # kJ/mol/nm^2 for all distance constraints
        RING_ANGLE_FORCE_CONSTANT = 1000.0  # kJ/mol/rad^2 for angle constraints
        
        # Build a mapping from (residue_number, atom_name) to atom_index
        # molecule_type.atoms fields: [index, type, resnum, resname, atomname, ...]
        residue_atom_map = {}  # {(resnum, resname): {atom_name: atom_index}}
        
        for atom_idx, fields in enumerate(molecule_type.atoms):
            if len(fields) < 5:
                continue
            resnum = fields[2]
            resname = fields[3]
            atomname = fields[4]
            
            # Normalize residue name (handle replacements)
            if resname in PDBFile._residueNameReplacements:
                resname = PDBFile._residueNameReplacements[resname]
            
            # Normalize atom name (handle replacements)
            if resname in PDBFile._atomNameReplacements:
                atomReplacements = PDBFile._atomNameReplacements[resname]
                if atomname in atomReplacements:
                    atomname = atomReplacements[atomname]
            
            # Store mapping (use uppercase for consistency with constraint data)
            key = (resnum, resname.upper())
            if key not in residue_atom_map:
                residue_atom_map[key] = {}
            residue_atom_map[key][atomname.upper()] = base_atom_index + atom_idx
        
        # Ensure harmonic_bond_force exists
        if self.harmonic_bond_force is None:
            self.harmonic_bond_force = mm.HarmonicBondForce()
            sys.addForce(self.harmonic_bond_force)
        
        # Ensure harmonic_angle_force exists for angle constraints
        if self.harmonic_angle_force is None:
            self.harmonic_angle_force = mm.HarmonicAngleForce()
            sys.addForce(self.harmonic_angle_force)
        
        # Track existing bonds to avoid duplicates
        existing_bonds = set()
        for i in range(self.harmonic_bond_force.getNumBonds()):
            bond = self.harmonic_bond_force.getBondParameters(i)
            atom1, atom2 = bond[0], bond[1]
            if atom1 < atom2:
                existing_bonds.add((atom1, atom2))
            else:
                existing_bonds.add((atom2, atom1))
        
        # Track existing angles to avoid duplicates
        existing_angles = set()
        for i in range(self.harmonic_angle_force.getNumAngles()):
            angle = self.harmonic_angle_force.getAngleParameters(i)
            atom1, atom2, atom3 = angle[0], angle[1], angle[2]
            # Store angles in canonical order (middle atom is always atom2)
            existing_angles.add((atom1, atom2, atom3))
        
        # Helper function to add constraints
        def add_constraint(constraint_list, atom_map, existing_bonds, is_bond=False):
            """Add constraints from a constraint list
            
            Parameters
            ----------
            constraint_list : list
                List of constraint dictionaries
            atom_map : dict
                Mapping from atom names to atom indices
            existing_bonds : set
                Set of existing bond keys to avoid duplicates
            is_bond : bool
                Unused parameter (kept for compatibility), all constraints use same force constant
            """
            force_constant = RING_FORCE_CONSTANT
            for constraint in constraint_list:
                atom1_name = constraint["atom1"].upper()
                atom2_name = constraint["atom2"].upper()
                distance_nm = constraint["distance_nm"]
                
                if atom1_name in atom_map and atom2_name in atom_map:
                    atom1_idx = atom_map[atom1_name]
                    atom2_idx = atom_map[atom2_name]
                    
                    # Avoid duplicates
                    bond_key = (min(atom1_idx, atom2_idx), max(atom1_idx, atom2_idx))
                    if bond_key not in existing_bonds:
                        self.harmonic_bond_force.addBond(
                            atom1_idx,
                            atom2_idx,
                            distance_nm * unit.nanometer,
                            force_constant * unit.kilojoule_per_mole / unit.nanometer**2
                        )
                        existing_bonds.add(bond_key)
        
        # Helper function to add angle constraints
        def add_angle_constraint(angle_list, atom_map, existing_angles):
            """Add angle constraints from an angle list
            
            Parameters
            ----------
            angle_list : list
                List of angle dictionaries with atom1, atom2, atom3, angle_deg
            atom_map : dict
                Mapping from atom names to atom indices
            existing_angles : set
                Set of existing angle tuples to avoid duplicates
            """
            degToRad = math.pi / 180.0
            for angle_constraint in angle_list:
                atom1_name = angle_constraint["atom1"].upper()
                atom2_name = angle_constraint["atom2"].upper()
                atom3_name = angle_constraint["atom3"].upper()
                angle_deg = angle_constraint["angle_deg"]
                
                if atom1_name in atom_map and atom2_name in atom_map and atom3_name in atom_map:
                    atom1_idx = atom_map[atom1_name]
                    atom2_idx = atom_map[atom2_name]
                    atom3_idx = atom_map[atom3_name]
                    
                    # Avoid duplicates
                    angle_key = (atom1_idx, atom2_idx, atom3_idx)
                    if angle_key not in existing_angles:
                        angle_rad = angle_deg * degToRad
                        self.harmonic_angle_force.addAngle(
                            atom1_idx,
                            atom2_idx,
                            atom3_idx,
                            angle_rad,
                            RING_ANGLE_FORCE_CONSTANT * unit.kilojoule_per_mole / unit.radian**2
                        )
                        existing_angles.add(angle_key)
        
        # Process each residue
        for (resnum, resname), atom_map in residue_atom_map.items():
            # Only process aromatic residues (resname is already uppercase)
            if resname not in RING_CONSTRAINTS:
                continue
            
            constraints = RING_CONSTRAINTS[resname]
            
            # Handle TRP special structure (six_membered and five_membered rings)
            if resname == "TRP":
                # Process six-membered ring
                if "six_membered" in constraints:
                    six_membered = constraints["six_membered"]
                    add_constraint(six_membered.get("constraints_12", []), atom_map, existing_bonds, is_bond=True)
                    add_constraint(six_membered.get("constraints_13", []), atom_map, existing_bonds, is_bond=False)
                    add_constraint(six_membered.get("constraints_14", []), atom_map, existing_bonds, is_bond=False)
                    add_angle_constraint(six_membered.get("angles_13", []), atom_map, existing_angles)
                
                # Process five-membered ring
                if "five_membered" in constraints:
                    five_membered = constraints["five_membered"]
                    add_constraint(five_membered.get("constraints_12", []), atom_map, existing_bonds, is_bond=True)
                    add_constraint(five_membered.get("constraints_13", []), atom_map, existing_bonds, is_bond=False)
                    add_constraint(five_membered.get("constraints_14", []), atom_map, existing_bonds, is_bond=False)
                    add_angle_constraint(five_membered.get("angles_13", []), atom_map, existing_angles)
            else:
                # Handle PHE, TYR, HIS
                # Add 1-2 constraints (ring bonds)
                add_constraint(constraints.get("constraints_12", []), atom_map, existing_bonds, is_bond=True)
                
                # Add 1-3 constraints
                add_constraint(constraints.get("constraints_13", []), atom_map, existing_bonds, is_bond=False)
                
                # Add 1-4 constraints
                add_constraint(constraints.get("constraints_14", []), atom_map, existing_bonds, is_bond=False)
                
                # Add angle constraints
                add_angle_constraint(constraints.get("angles_13", []), atom_map, existing_angles)

    def _set_nonbonded_parameters(
        self,
        molecule_type,
        base_atom_index,
        atom_types,
        atom_type_map,
    ):
        charges = []
        for i, fields in enumerate(molecule_type.atoms):
            params = self._atom_types[fields[1]]
            if len(fields) > 6:
                q = float(fields[6])
            else:
                q = float(params[4])

            charges.append(q)
            index = base_atom_index + i
            atomType = atom_type_map[params[0]]
            self.nb_force.addParticle([atomType, q])
            # Add in the self term for the reaction field correction
            # PACE do not use reaction field correction
            #self.es_self_excl_force.addBond(index, index, [0.5 * q * q])

        pairExceptions = self._gen_pair_exceptions(
            molecule_type, base_atom_index, atom_types
        )
        exclusion_exceptions = self._gen_exclusion_exceptions(
            molecule_type, base_atom_index, atom_types
        )
        bondedExceptions = self._gen_bonded_exceptions(molecule_type, base_atom_index)

        angleExceptions = self._gen_angle_exceptions(molecule_type, base_atom_index)

        dihedralExceptions = self._gen_dihedral_exceptions(molecule_type, base_atom_index)

        '''
        constraint_exceptions = self._gen_constraint_exceptions(
            molecule_type, base_atom_index
        )
        '''
        # order matters
        '''
        exceptions = (
            pairExceptions
            + exclusion_exceptions
            + bondedExceptions
            #+ constraint_exceptions
        )
        '''
        #print (pairExceptions)
        exceptions = (
            bondedExceptions
            + angleExceptions
            + dihedralExceptions
            + exclusion_exceptions
            + pairExceptions
        )
        
        return exceptions, charges

    def _gen_pair_exceptions(self, molecule_type, base_atom_index, atom_types):
        pairExceptions = []
        for fields in molecule_type.pairs:
            atoms = [int(x) - 1 for x in fields[:2]]
            types = tuple(atom_types[i] for i in atoms)
            if len(fields) >= 5:
                params = fields[3:5]
            elif types in self._pairTypes:
                params = self._pairTypes[types][3:5]
            elif types[::-1] in self._pairTypes:
                params = self._pairTypes[types[::-1]][3:5]
            else:
                # We'll use the automatically generated parameters
                continue

            type1_, q1 = self.nb_force.getParticleParameters(base_atom_index + atoms[0])
            type2_, q2 = self.nb_force.getParticleParameters(base_atom_index + atoms[1])
            pairExceptions.append(
                (
                    base_atom_index + atoms[0],
                    base_atom_index + atoms[1],
                    q1 * q2,
                    params[0],
                    params[1],
                    0
                )
            )
        return pairExceptions

    def _gen_exclusion_exceptions(self, molecule_type, base_atom_index, atom_types):
        exclusion_exceptions = []
        for fields in molecule_type.exclusions:
            atoms = [int(x) - 1 for x in fields]
            for atom in atoms[1:]:
                if atom > atoms[0]:
                    type1_, q1 = self.nb_force.getParticleParameters(base_atom_index + atoms[0])
                    type2_, q2 = self.nb_force.getParticleParameters(base_atom_index + atom)
                    exclusion_exceptions.append(
                        (base_atom_index + atoms[0], base_atom_index + atom, 0, 0, 0, 1)
                    )
        return exclusion_exceptions

    def _gen_bonded_exceptions(self, molecule_type, base_atom_index):
        bond_indices = []
        for fields in molecule_type.bonds:
            atoms = [int(x) - 1 for x in fields[:2]]
            bond_type = int(fields[2])
            # Type 6 bonds do not generate exclusions
            if bond_type == 1:
                bond_indices.append(
                    (base_atom_index + atoms[0], base_atom_index + atoms[1])
                )

        bondedExceptions = [(i, j, 0, 0, 0, 2) for i, j in bond_indices]
        return bondedExceptions

    # Add angel exceptions for PACE

    def _gen_angle_exceptions(self, molecule_type, base_atom_index):
        bond_indices = []
        for fields in molecule_type.harmonic_angles:
            atoms = [int(x) - 1 for x in fields[:3]]
            bond_type = int(fields[3])
            #In PACE angle type is always type '1'
            if bond_type == 1:
                bond_indices.append(
                    (base_atom_index + atoms[0], base_atom_index + atoms[2])
                )
        
        angleExceptions = [(i, j, 0, 0, 0, 3) for i, j in bond_indices]
        return angleExceptions

    def _gen_dihedral_exceptions(self, molecule_type, base_atom_index):
        bond_indices = []
        for fields in molecule_type.dihedrals:
            atoms = [int(x) - 1 for x in fields[:4]]
            bond_type = int(fields[4])
            #In PACE angle type is always type '1'
            if bond_type == 1:
                bond_indices.append(
                    (base_atom_index + atoms[0], base_atom_index + atoms[3])
                )
        
        dihedralExceptions = [(i, j, 0, 0, 0, 4) for i, j in bond_indices]
        return dihedralExceptions
    '''
    # We do not use it as PACE do not have constrains
    def _gen_constraint_exceptions(self, molecule_type, base_atom_index):
        constraint_indices = []
        for fields in molecule_type.constraints:
            atoms = [int(x) - 1 for x in fields[:2]]
            constraint_indices.append(
                (base_atom_index + atoms[0], base_atom_index + atoms[1])
            )

        constraint_exceptions = [(i, j, 0, 0, 0) for i, j in constraint_indices]
        return constraint_exceptions
    '''

    def create_system(
        self,
        nonbonded_cutoff=1.1 * unit.nanometer,
        remove_com_motion=True,
        add_nonbonded_force=True,
        nonbonded_type="standard",
        add_ring_constraints=False,
    ):
        """Construct an OpenMM System representing the topology described by this
        top file.

        Parameters
        ----------
        nonbonded_cutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        remove_com_motion : boolean=True
            If true, a CMMotionRemover will be added to the System
        add_nonbonded_force : boolean=True
            If true, nonbonded forces (LJ and electrostatic) will be added to the System
        nonbonded_type : string="standard"
            Type of nonbonded force to use:
            - "standard": Original Martini nonbonded force with LJ and electrostatics
            - "gaussian": Gaussian "Bulldozer" Force for initial untangling
            - "softcore": Gapsys Soft-Core Force for smooth potential energy landscapes
        add_ring_constraints : boolean=False
            If true, add rigid triangle constraints for aromatic rings (PHE, TYR, TRP).
            Default is False. Typically enabled only in final optimization steps.
        Returns
        -------
        System
             the newly created System
        """
        sys = mm.System()
        box_vectors = self.topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            sys.setDefaultPeriodicBoxVectors(*box_vectors)
        else:
            raise ValueError("periodicBoxVectors must be set")

        # Only add nonbonded forces if requested

        # Custom force to handle martini nonbonded terms
        # nbfix-like terms
        
        # Create the appropriate type of nonbonded force based on nonbonded_type
        if nonbonded_type == "gaussian":
            # Stage 1: Gaussian "Bulldozer" Force
            # Used for initial untangling. Ignores chemistry, only geometric repulsion.
            # Formula: E = Height * exp( - (r / Width)^2 )
            self.nb_force = mm.CustomNonbondedForce(
                "ga_k * ga_h * exp(-(r/ga_w)^2);"
            )
            # Global parameters (can be tuned in Context)
            self.nb_force.addGlobalParameter("ga_k", 1.0)     # ON/OFF switch
            self.nb_force.addGlobalParameter("ga_h", 800.0)   # Height (kJ/mol)
            self.nb_force.addGlobalParameter("ga_w", 0.14)    # Width (nm) - 1.4 Å
            
            # Dummy parameters to match the addParticle signature in _set_nonbonded_parameters
            # _set_nonbonded_parameters calls addParticle([atomType, q])
            self.nb_force.addPerParticleParameter("dummy_type") 
            self.nb_force.addPerParticleParameter("dummy_q")
            
            self.nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            self.nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
            sys.addForce(self.nb_force)
            
            # Note: In Gaussian mode, we do NOT add exception forces (es/lj_except_force)
            # because we want pure geometric repulsion to untangle knots.
            self.es_except_force = None
            self.lj_except_force = None
            
        elif nonbonded_type == "softcore":
            # Stage 2: Gapsys Soft-Core Force (C6/C12 Version)
            # Uses Tabulated Functions for C6/C12 but modifies the expression for softness.
            # Pre-compute constants to avoid repeated method calls
            rcut_val = nonbonded_cutoff.value_in_unit(unit.nanometers)
            epsilon_r_val = self.epsilon_r
            self.nb_force = mm.CustomNonbondedForce(
                # 1. 主公式 (包含所有分量)
                "step(rcut - r) * (LJ_soft - corrLJ + ES_soft);"

                # 2. 静电项 (依赖 crf 和 inv_r_soft)
                "ES_soft = f * q1 * q2 / epsilon_r * (inv_r_soft + (1 / (2 * rcut^3)) * r^2 - crf);"
                "crf = 1 / rcut + (1 / (2 * rcut^3)) * rcut^2;"
                "inv_r_soft = 1 / sqrt(r^2 + soft_alpha_coul);"

                # 3. LJ项 (依赖 reff_sq, reff, corrLJ)
                "LJ_soft = (C12(type1, type2) / reff_sq) - (C6(type1, type2) / reff);"
                "corrLJ = step(r - rswitch) * (C12(type1, type2)/r^12 - C6(type1, type2)/r^6) * switching_func;"
                "switching_func = (10 * x_val^3 - 15 * x_val^4 + 6 * x_val^5);"
                "x_val = (r - rswitch) / (rcut - rswitch);"

                # 4. 软核中间变量 (依赖 sigma6_proxy)
                "reff_sq = reff * reff;"
                "reff = r^6 + soft_alpha * sigma6_proxy;"
                "sigma6_proxy = select(step(1e-8 - abs(C6(type1, type2))), 1.0, C12(type1, type2) / C6(type1, type2));"

                # 5. 常数
                "f = 138.935458;"
                "rswitch = 0.9;"
                f"epsilon_r = {epsilon_r_val};"
                f"rcut = {rcut_val};"
            )
                    
            # Parameters
            self.nb_force.addGlobalParameter("soft_alpha", 0.5)      # VDW Softness
            self.nb_force.addGlobalParameter("soft_alpha_coul", 0.1) # Electrostatic Softness
            
            self.nb_force.addPerParticleParameter("type")
            self.nb_force.addPerParticleParameter("q")
            
            self.nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            self.nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
            sys.addForce(self.nb_force)
            
            # For softcore mode, we still need exception forces
            self.es_except_force = mm.CustomBondForce(
                f"step(rcut-r) * (ES -corrES);"
                f"ES = f*qprod/epsilon_r * (1/r);"
                f"corrES = f*qprod/epsilon_r * (1/rcut);" 
                #f"corrES = 0;" 
                f"epsilon_r = {self.epsilon_r};"
                f"f = 138.935458;"
                f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
            )
            self.es_except_force.addPerBondParameter("qprod")
            sys.addForce(self.es_except_force)
            
            # Softened exception force for LJ interactions
            self.lj_except_force = mm.CustomBondForce(
                "step(rcut-r) * (energy - corr);"
                "energy = C12 / (r^6 + soft_alpha_exc)^2 - C6 / (r^6 + soft_alpha_exc);"
                "corr = C12 / (rcut^6 + soft_alpha_exc)^2 - C6 / (rcut^6 + soft_alpha_exc);"
                #"corr = 0;"
                f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
            )
            self.lj_except_force.addGlobalParameter("soft_alpha_exc", 0.5)  # Exception force softness
            self.lj_except_force.addPerBondParameter("C12")
            self.lj_except_force.addPerBondParameter("C6")
            sys.addForce(self.lj_except_force)
            
        else:  # "standard" (Original code logic)
            self.nb_force = mm.CustomNonbondedForce(
                "step(rcut-r)*(LJ - corrLJ + ES);"
                "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
                "corrLJ = step(r-rswitch) * (C12(type1, type2) / r^12 - C6(type1, type2) / r^6) * (10 *x^3 - 15 * x^4 + 6 * x^5);"
                "x = (r - rswitch)/(rcut - rswitch);"
                "rswitch = 0.9;"
                "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
                "crf = 1 / rcut + (1 / (2 * rcut^3))* rcut^2;"
                "krf = 1 / (2 * rcut^3);"
                f"epsilon_r = {self.epsilon_r};"
                "f = 138.935458;"
                f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
            )
            #self.nb_force.addComputedValue('k', 'step(0.5 - soft) * 1 + step(soft - 0.5) * (1 + cos(3.14159265359 * r / rcut)) ')
            self.nb_force.addGlobalParameter("soft", 1)
            self.nb_force.addPerParticleParameter("type")
            self.nb_force.addPerParticleParameter("q")
            self.nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            self.nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
            sys.addForce(self.nb_force)

            # custom non-bonded force for the ES exceptions
            
            self.es_except_force = mm.CustomBondForce(
                f"step(rcut-r) * (ES -corrES);"
                f"ES = f*qprod/epsilon_r * (1/r);"
                f"corrES = f*qprod/epsilon_r * (1/rcut);" 
                #f"corrES = 0;" 
                f"epsilon_r = {self.epsilon_r};"
                f"f = 138.935458;"
                f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
            )
            self.es_except_force.addPerBondParameter("qprod")
            sys.addForce(self.es_except_force)
            # custom bonded force to handle exceptions

            
            self.lj_except_force = mm.CustomBondForce(
                "step(rcut-r) * (energy - corr);"
                "energy = (C12/r^12 - C6/r^6);"
                "corr = (C12/rcut^12 - C6/rcut^6);"
                #"corr = 0;"
                f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
            )
            self.lj_except_force.addPerBondParameter("C12")
            self.lj_except_force.addPerBondParameter("C6")
            sys.addForce(self.lj_except_force)
        

        all_exceptions = []
        all_charges = []

        (
            dihedral_type_table,
            wildcard_dihedral_types,
        ) = self._build_dihedral_lookup_table()

        # build a lookup table mapping atom types into integer indices
        # that will later be used to lookup LJ combinations
        used_atom_types = set()
        for molecule_name, _ in self._molecules:
            molecule_type = self._moleculeTypes[molecule_name]
            for atom in molecule_type.atoms:
                used_atom_types.add(atom[1])
        atom_type_map = {k: i for i, k in enumerate(sorted(used_atom_types))}

        # now we need to setup our tables of C6 and C12

        n_types = len(atom_type_map)
        C6 = []
        C12 = []
        for type_i in atom_type_map:
            for type_j in atom_type_map:
                type_i_sorted, type_j_sorted = sorted([type_i, type_j])
                if (type_i_sorted, type_j_sorted) in self._nonbond_types:
                    # Parameters were specified for this pair of types,
                    # so use them.
                    params = self._nonbond_types[(type_i_sorted, type_j_sorted)]
                    if self._use_sigma_eps:
                        sigma = float(params[3])
                        eps = float(params[4])
                        c6 = 4 * eps * sigma ** 6
                        c12 = 4 * eps * sigma ** 12
                    else:
                        c6 = float(params[3])
                        c12 = float(params[4])
                else:
                    # We do not use combination rules in PACE
                    # Parameters were not specified for this pair type,
                    # so calculate using combination rules.
                    c6 = 0
                    c12 = 0
                C6.append(c6)
                C12.append(c12)

        if add_nonbonded_force:
            self.nb_force.addTabulatedFunction(
                "C6", mm.Discrete2DFunction(n_types, n_types, C6)
            )
            self.nb_force.addTabulatedFunction(
                "C12", mm.Discrete2DFunction(n_types, n_types, C12)
            )
        else: 
            C6_0 = [0] * len(C6)
            C12_0 = [0] * len(C12)
            self.nb_force.addTabulatedFunction(
                "C6", mm.Discrete2DFunction(n_types, n_types, C6_0)
            )
            self.nb_force.addTabulatedFunction(
                "C12", mm.Discrete2DFunction(n_types, n_types, C12_0)
            )            

        # Loop over molecules and create the specified number of each type.
        for molecule_name, molecule_count in self._molecules:
            molecule_type = self._moleculeTypes[molecule_name]
            for i in range(molecule_count):
                base_atom_index = sys.getNumParticles()
                atom_types = [atom[1] for atom in molecule_type.atoms]
                try:
                    bonded_types = [self._atom_types[t][1] for t in atom_types]
                except KeyError as e:
                    raise ValueError("Unknown atom type: " + e.message)
                bonded_types = [
                    b if b is not None else a for a, b in zip(atom_types, bonded_types)
                ]

                # add bonded parameters
                self._add_atoms_to_system(sys, molecule_type)
                self._add_bonds_to_system(sys, molecule_type, bonded_types, base_atom_index)
                
                self._add_angles_to_system(
                    sys, molecule_type, bonded_types, base_atom_index
                )
                
                self._add_torsions_to_system(
                    sys,
                    molecule_type,
                    bonded_types,
                    dihedral_type_table,
                    wildcard_dihedral_types,
                    base_atom_index,
                )
                
                self._add_cmap_to_system(
                    sys, molecule_type, bonded_types, base_atom_index
                )
                #self._add_vsites_to_system(sys, molecule_type, base_atom_index)

                exceptions, charges = self._set_nonbonded_parameters(
                    molecule_type,
                    base_atom_index,
                    atom_types,
                    atom_type_map,
                )

                all_exceptions.extend(exceptions)
                all_charges.extend(charges)

                # Add explicitly specified constraints.
                for fields in molecule_type.constraints:
                    atoms = [int(x) - 1 for x in fields[:2]]
                    constraint_type = int(fields[2])
                    assert constraint_type == 1
                    length = float(fields[3])
                    sys.addConstraint(
                        base_atom_index + atoms[0], base_atom_index + atoms[1], length
                    )
                
                # Add ring constraints for aromatic residues if enabled
                if add_ring_constraints:
                    self._add_ring_constraints(sys, molecule_type, base_atom_index)

        if all_exceptions:
            # build a map of unique exceptions
            # process in order, so that later entries trump earlier ones
            except_map = defaultdict(list)
            #print (except_map)
            for exception in all_exceptions:
               
                # exp_type 
                # 0 pairs 1 exclusion
                # 2 bond 3 angle 4 dihedral 
                # 
                i, j, q, c6, c12, exc_type = exception
                if exc_type == 4:
                    dihedral_flag = 1
                else:
                    dihedral_flag = 0
                if i < j:
                    if len(except_map[(i, j)]) != 0:
                        origin_dihedral_flag = except_map[(i, j)][-1]
                    else: 
                        origin_dihedral_flag = 0
                else:
                    if len(except_map[(j, i)]) != 0:
                        origin_dihedral_flag = except_map[(j, i)][-1]                    
                    else:
                        origin_dihedral_flag = 0
                if origin_dihedral_flag == 1:
                    dihedral_flag = 1
                
                if i < j:
                    except_map[(i, j)] = [q, c6, c12, exc_type, dihedral_flag]
                else:
                    except_map[(j, i)] = [q, c6, c12, exc_type, dihedral_flag]
     
            # add in all of the exclusions
            for i, j in except_map:

                c6 = float(except_map[(i, j)][1])
                c12 = float(except_map[(i, j)][2])
                q = float(except_map[(i, j)][0])
                exc_type = except_map[(i, j)][3]
                dihedral_flag = except_map[(i, j)][4]
                #print (i, j , c6, c12 , q )
                # Remove i,j from nonbonded interactions for all exceptions / exclusions
                
                self.nb_force.addExclusion(i, j)

                # Handle electrostatic exceptions / exclusions.
                # We're going to assume that q==0 means that this was an
                # exclusion.
                # Note: In gaussian mode, es_except_force is None, so skip these

                if self.es_except_force is not None:
                    if exc_type in [0,4] and dihedral_flag == 1:  
                        if q == 0:
                            # In this case, we still need to add in the reaction field correction
                            # term.
                            qprod = all_charges[i] * all_charges[j]
                            # Don't bother adding the interaction if either particle has zero charge.
                            #if qprod != 0.0:
                            if qprod != 0:
                                self.es_except_force.addBond(i, j, [qprod])
                        # If q !=0, then this is an exception.
                        else:
                        # We don't bother adding the interaction if either the charge product is zero
                            self.es_except_force.addBond(i, j, [q])
                
                # Now we'll add the LJ exceptions. We don't bother
                # adding the interaction if the combined LJ parameters are zero: In gaussian mode.
                # Note, lj_except_force is None, so skip these
                if self.lj_except_force is not None:
                    if exc_type in [0]  :   
                        #if c6 != 0 or c12 != 0:
                            self.lj_except_force.addBond(i, j, [c12, c6])
                        
                    
                    
        # When add_nonbonded_force is False, exclude all nonbonded interactions
        # (both LJ and electrostatic) to allow bond-only minimization
        if not add_nonbonded_force:
            # Set all particle charges to 0 to exclude electrostatic interactions
            # nb_force has per-particle parameters: type (index 0) and q (index 1)
            for particle_index in range(self.nb_force.getNumParticles()):
                params = list(self.nb_force.getParticleParameters(particle_index))
                params[1] = 0  # Set charge to 0
                self.nb_force.setParticleParameters(particle_index, params)
            print(f"Excluded electrostatic interactions for {self.nb_force.getNumParticles()} particles")

        if remove_com_motion:
            sys.addForce(mm.CMMotionRemover())
        return sys




def _get_default_gromacs_include_dir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) for the gromacs binary
    'pdb2gmx' or 'gmx' in the PATH, or (3) just using the default gromacs
    install location, /usr/local/gromacs/share/gromacs/top"""
    if "GMXDATA" in os.environ:
        return os.path.join(os.environ["GMXDATA"], "top")
    if "GMXBIN" in os.environ:
        return os.path.abspath(
            os.path.join(os.environ["GMXBIN"], "..", "share", "gromacs", "top")
        )

    pdb2gmx_path = shutil.which("pdb2gmx")
    if pdb2gmx_path is not None:
        return os.path.abspath(
            os.path.join(os.path.dirname(pdb2gmx_path), "..", "share", "gromacs", "top")
        )
    else:
        gmx_path = shutil.which("gmx")
        if gmx_path is not None:
            return os.path.abspath(
                os.path.join(os.path.dirname(gmx_path), "..", "share", "gromacs", "top")
            )

    return "/usr/local/gromacs/share/gromacs/top"