ring_distance_dict = {'HIS': {'ring': {'CD2-CE1': 2.1841,
                'CD2-ND1': 2.1878,
                'CE1-CG': 2.191,
                'CG-NE2': 2.1907,
                'ND1-NE2': 2.1715}},
'PHE': {'ring': {'CD1-CD2': 2.407,
                'CD1-CE2': 2.7805,
                'CD1-CZ': 2.408,
                'CD2-CE1': 2.7794,
                'CD2-CZ': 2.4073,
                'CE1-CE2': 2.4079,
                'CE1-CG': 2.4075,
                'CE2-CG': 2.4072,
                'CG-CZ': 2.78}},
'TRP': {'five_membered': {'CD1-CD2': 2.2398,
                        'CD1-CE2': 2.258,
                        'CD2-NE1': 2.2095,
                        'CE2-CG': 2.2965,
                        'CG-NE1': 2.2316},
        'six_membered': {'CD2-CH2': 2.7998,
                        'CD2-CZ2': 2.425,
                        'CD2-CZ3': 2.425,
                        'CE2-CE3': 2.3957,
                        'CE2-CH2': 2.4244,
                        'CE2-CZ3': 2.75,
                        'CE3-CH2': 2.3675,
                        'CE3-CZ2': 2.7504,
                        'CZ2-CZ3': 2.3393}},
'TYR': {'ring': {'CD1-CD2': 2.3852,
                'CD1-CE2': 2.742,
                'CD1-CZ': 2.4075,
                'CD2-CE1': 2.7418,
                'CD2-CZ': 2.3642,
                'CE1-CE2': 2.3416,
                'CE1-CG': 2.4078,
                'CE2-CG': 2.4077,
                'CG-CZ': 2.7802}}}

# 键角定义：atom1-atom2-atom3 的理想角度（度）
# 六元环内角：120°，五元环内角：108°
ring_angle_dict = {'HIS': {'ring': {
    'ND1-CE1-NE2': 108.0,
    'CE1-NE2-CD2': 108.0,
    'NE2-CD2-CG': 108.0,
    'CD2-CG-ND1': 108.0,
    'CG-ND1-CE1': 108.0,
    'ND1-CD2-CG': 108.0,
    'CD2-ND1-CE1': 108.0,
    'CD2-CE1-NE2': 108.0
}},
'PHE': {'ring': {
    'CG-CD1-CE1': 120.0,
    'CD1-CE1-CZ': 120.0,
    'CE1-CZ-CE2': 120.0,
    'CZ-CE2-CD2': 120.0,
    'CE2-CD2-CG': 120.0,
    'CD2-CG-CD1': 120.0,
    'CD1-CD2-CE2': 120.0,
    'CD1-CZ-CE2': 120.0,
    'CE1-CD1-CD2': 120.0,
    'CE1-CZ-CD2': 120.0
}},
'TYR': {'ring': {
    'CG-CD1-CE1': 120.0,
    'CD1-CE1-CZ': 120.0,
    'CE1-CZ-CE2': 120.0,
    'CZ-CE2-CD2': 120.0,
    'CE2-CD2-CG': 120.0,
    'CD2-CG-CD1': 120.0,
    'CD1-CD2-CE2': 120.0,
    'CD1-CZ-CE2': 120.0,
    'CE1-CD1-CD2': 120.0,
    'CE1-CZ-CD2': 120.0
}},
'TRP': {'five_membered': {
    'CG-CD1-NE1': 108.0,
    'CD1-NE1-CE2': 108.0,
    'NE1-CE2-CD2': 108.0,
    'CE2-CD2-CG': 108.0,
    'CD2-CG-CD1': 108.0,
    'CD1-CE2-NE1': 108.0,
    'CD1-CD2-CE2': 108.0,
    'NE1-CD2-CG': 108.0
},
'six_membered': {
    'CG-CD2-CE3': 120.0,
    'CD2-CE3-CZ3': 120.0,
    'CE3-CZ3-CH2': 120.0,
    'CZ3-CH2-CZ2': 120.0,
    'CH2-CZ2-CE2': 120.0,
    'CZ2-CE2-CD2': 120.0,
    'CD2-CZ3-CE3': 120.0,
    'CD2-CH2-CZ2': 120.0,
    'CD2-CE3-CZ2': 120.0,
    'CE3-CD2-CZ2': 120.0,
    'CE3-CH2-CZ3': 120.0
}}}

ring_dihedral_dict = {'HIS': {'CD2-CG-ND1-CE1': 0,
            'CE1-NE2-CD2-CG': 0,
            'CG-ND1-CE1-NE2': 0,
            'ND1-CE1-NE2-CD2': 0,
            'NE2-CD2-CG-ND1': 0},
    'PHE': {'CD1-CE1-CZ-CE2': 0,
            'CD2-CG-CD1-CE1': 0,
            'CE1-CZ-CE2-CD2': 0,
            'CE2-CD2-CG-CD1': 0,
            'CG-CD1-CE1-CZ': 0,
            'CZ-CE2-CD2-CG': 0},
    'TRP': {'CD1-CG-CD2-CE3': 180,
            'CD1-NE1-CE2-CD2': 0,
            'CD1-NE1-CE2-CZ2': 180,
            'CD2-CE2-CZ2-CH2': 0,
            'CD2-CG-CD1-NE1': 0,
            'CE2-CD2-CG-CD1': 0,
            'CE2-CZ2-CH2-CZ3': 0,
            'CE3-CD2-CE2-CZ2': 0,
            'CG-CD1-NE1-CE2': 0,
            'CG-CD2-CE3-CZ3': 180,
            'CH2-CZ3-CE3-CD2': 0,
            'CZ2-CE2-CD2-CG': 180,
            'CZ2-CH2-CZ3-CE3': 0,
            'CZ3-CE3-CD2-CE2': 0,
            'NE1-CE2-CD2-CE3': 180,
            'NE1-CE2-CD2-CG': 0,
            'NE1-CE2-CZ2-CH2': 180},
    'TYR': {'CD1-CE1-CZ-CE2': 0,
            'CD2-CG-CD1-CE1': 0,
            'CE1-CZ-CE2-CD2': 0,
            'CE2-CD2-CG-CD1': 0,
            'CG-CD1-CE1-CZ': 0,
            'CZ-CE2-CD2-CG': 0}}

# 对角线距离限制字典：用于United-Atom模型中维持芳香环平面结构
# 格式: {残基名: {环类型: {(原子对): (low, up1, up2, kdr)}}}
# Type 10 (Distance Restraints): low与up1之间为零电势，up2后施加惩罚力常数kdr
ring_diagonal_dict = {
    'PHE': {
        'ring': {
            ('CG', 'CZ'): (0.270, 0.290, 0.310, 125000),
            ('CD1', 'CE2'): (0.270, 0.290, 0.310, 125000),
            ('CD2', 'CE1'): (0.270, 0.290, 0.310, 125000)
        }
    },
    'TYR': {
        'ring': {
            ('CG', 'CZ'): (0.270, 0.290, 0.310, 125000),
            ('CD1', 'CE2'): (0.270, 0.290, 0.310, 125000),
            ('CD2', 'CE1'): (0.270, 0.290, 0.310, 125000)
        }
    },
    'HIS': {
        'ring': {
            ('CG', 'CE1'): (0.215, 0.235, 0.255, 125000),
            ('CG', 'NE2'): (0.215, 0.235, 0.255, 125000),
            ('ND1', 'NE2'): (0.215, 0.235, 0.255, 125000),
            ('ND1', 'CD2'): (0.215, 0.235, 0.255, 125000),
            ('CE1', 'CD2'): (0.215, 0.235, 0.255, 125000)
        }
    },
    'TRP': {
        'six_membered': {
            ('CD2', 'CH2'): (0.270, 0.290, 0.310, 125000),
            ('CE2', 'CZ3'): (0.270, 0.290, 0.310, 125000),
            ('CZ2', 'CE3'): (0.270, 0.290, 0.310, 125000)
        },
        'five_membered': {
            ('CG', 'NE1'): (0.220, 0.240, 0.260, 125000),
            ('CG', 'CE2'): (0.220, 0.240, 0.260, 125000),
            ('CD1', 'CE2'): (0.220, 0.240, 0.260, 125000),
            ('CD1', 'CD2'): (0.220, 0.240, 0.260, 125000),
            ('NE1', 'CD2'): (0.220, 0.240, 0.260, 125000)
        }
    }
}

import sys 
input_topology = sys.argv[1]

class GromacsAtom:
    def __init__(self, input_string):
        parts = input_string.split()
        if len(parts) >= 8:
            self.nr = int(parts[0])               
            self.atom_type = parts[1]             
            self.resnr = int(parts[2])            
            self.residue = parts[3]               
            self.atom_name = parts[4]             
            self.charge = float(parts[6])         
            self.mass = float(parts[7])           
        else:
            raise ValueError("Invalid input string format")

    def __str__(self):
        return f"Atom {self.nr}: {self.atom_name} ({self.atom_type}), " \
               f"Residue {self.resnr} ({self.residue}), " \
               f"Charge: {self.charge}, Mass: {self.mass}"


def atoms_in_residue(atoms, nres):
    residues = [[] for i in range(nres)]
    residueNames = []
    pre_resnr = 0
    for atom in atoms:
        residues[atom.resnr - 1].append(atom)
        if pre_resnr != atom.resnr:
            pre_resnr = atom.resnr
            residueNames.append(atom.residue)
    return residues, residueNames

f = open(input_topology, 'r')
atomFlag = False
atomlist = []
for rl in f:
    srl = rl.split()
    if rl.startswith(';') or len(srl) == 0:
        continue
    if rl.startswith('[ atoms'):
        atomFlag = True
        continue
    if rl.startswith('[ bonds'):
        atomFlag = False
    if atomFlag:    
        atomlist.append(GromacsAtom(rl))
f.close()


# example for bond: 
# atom1.nr atom2.nr 1 r_bond 125000 
# example for dihedral: 
# atom1.nr atom2.nr atom3.nr atom4.nr 1 theta_dihedral force_constant 1 

def find_atom_by_name(atoms, residue_name, atom_name):
    """根据残基名和原子名查找原子"""
    for atom in atoms:
        if atom.residue == residue_name and atom.atom_name == atom_name:
            return atom
    return None

def generate_ring_angles(atoms, residues, residue_names, force_constant):
    """生成环状氨基酸的键角约束"""
    angles = []
    
    for i, residue_name in enumerate(residue_names):
        if residue_name in ring_angle_dict:
            residue_atoms = residues[i]
            angle_data = ring_angle_dict[residue_name]
            
            for angle_name, ideal_angle in angle_data.items():
                atom_names = angle_name.split('-')
                if len(atom_names) == 3:
                    atoms_found = []
                    for atom_name in atom_names:
                        atom = find_atom_by_name(residue_atoms, residue_name, atom_name)
                        if atom:
                            atoms_found.append(atom)
                        else:
                            break
                    
                    if len(atoms_found) == 3:
                        # 使用GROMACS的键角约束格式: atom1 atom2 atom3 1 theta force_constant
                        # 力常数通常设为400 kJ/mol/rad^2
                        angles.append(f"{atoms_found[0].nr:5d} {atoms_found[1].nr:5d} {atoms_found[2].nr:5d} 1 {ideal_angle:8.2f} {force_constant:8.0f}")
    
    return angles

def generate_ring_dihedrals(atoms, residues, residue_names, force_constant):
    """生成环状氨基酸的二面角约束"""
    dihedrals = []
    
    for i, residue_name in enumerate(residue_names):
        if residue_name in ring_dihedral_dict:
            residue_atoms = residues[i]
            dihedral_data = ring_dihedral_dict[residue_name]
            
            for dihedral_name, angle in dihedral_data.items():
                atom_names = dihedral_name.split('-')
                if len(atom_names) == 4:
                    atoms_found = []
                    for atom_name in atom_names:
                        atom = find_atom_by_name(residue_atoms, residue_name, atom_name)
                        if atom:
                            atoms_found.append(atom)
                        else:
                            break
                    
                    if len(atoms_found) == 4:
                        # 使用GROMACS的二面角约束格式: atom1 atom2 atom3 atom4 1 angle force_constant 1
                        # 力常数通常设为1000 kJ/mol
                        dihedrals.append(f"{atoms_found[0].nr:5d} {atoms_found[1].nr:5d} {atoms_found[2].nr:5d} {atoms_found[3].nr:5d} 1 {angle:8.1f} {force_constant:8.0f} 1")
    
    return dihedrals

def generate_ring_diagonals(atoms, residues, residue_names):
    """生成环状氨基酸的对角线距离限制 (Type 10 Distance Restraints)
    
    用于United-Atom模型中维持芳香环平面结构:
    - low与up1之间为零电势，允许自然热涨落
    - 超过up1后施加力常数kdr的惩罚
    """
    diagonals = []
    
    for i, residue_name in enumerate(residue_names):
        if residue_name in ring_diagonal_dict:
            residue_atoms = residues[i]
            diagonal_data = ring_diagonal_dict[residue_name]
            
            for ring_type, diagonal_pairs in diagonal_data.items():
                for (atom1_name, atom2_name), (low, up1, up2, kdr) in diagonal_pairs.items():
                    atom1 = find_atom_by_name(residue_atoms, residue_name, atom1_name)
                    atom2 = find_atom_by_name(residue_atoms, residue_name, atom2_name)
                    
                    if atom1 and atom2:
                        # 使用GROMACS的Type 10距离限制格式: atom1 atom2 10 low up1 up2 kdr
                        # funct=10 表示Distance Restraints
                        diagonals.append(f"{atom1.nr:5d} {atom2.nr:5d} 10 {low:6.3f} {up1:6.3f} {up2:6.3f} {kdr:6.0f}")
    
    return diagonals

def get_itp_lines(angles, dihedrals, diagonals=None):
    """将键角约束、二面角约束和对角线距离限制写入itp文件
    
    diagonals: 对角线距离限制列表 (Type 10 Distance Restraints)
    
    控制逻辑:
    - -DRINGREFINE: 控制键角、二面角约束
    - -DRINGREFINE_HEAVY: 控制对角线距离限制 (Type 10 Distance Restraints)
    """
    if diagonals is None:
        diagonals = []
    
    alllines = []
    
    # 对角线距离限制 - 由 RINGREFINE_HEAVY 控制
    if diagonals:
        alllines.append("#ifdef RINGREFINE_HEAVY")
        alllines.append("; Diagonal Distance Restraints for aromatic rings (Type 10)")
        alllines.append("; Generated by ring_refine.py")
        alllines.append("; Purpose: Maintain planar structure of aromatic rings in UA model")
        alllines.append("; Format: atom1 atom2 funct low up1 up2 kdr")
        alllines.append("; low-up1: zero potential (natural thermal fluctuations)")
        alllines.append("; >up1: penalty with force constant kdr")
        alllines.append("[ bonds ]")
        alllines.append("; ai   aj    funct  low (nm)  up1 (nm)  up2 (nm)  kdr (kJ/mol/nm^2)")
        for diagonal in diagonals:
            alllines.append(diagonal)
        alllines.append("#endif")
        alllines.append("")  # 空行分隔
    
    # 键角、二面角约束 - 由 RINGREFINE 控制
    alllines.append("#ifdef RINGREFINE")
    alllines.append("; Ring angle and dihedral constraints for aromatic amino acids")
    alllines.append("; Generated by ring_refine.py")

    if angles:
        alllines.append("; Ring angles")
        alllines.append("[ angles ]")
        alllines.append("; atom1 atom2 atom3 type angle force_constant")
        for angle in angles:
            alllines.append(angle)
    
    if dihedrals:
        alllines.append("; Ring dihedrals")
        alllines.append("[ dihedrals ]")
        alllines.append("; atom1 atom2 atom3 atom4 type angle force_constant multiplicity")
        for dihedral in dihedrals:
            alllines.append(dihedral)
        
    alllines.append("#endif")
    return alllines
# 主程序执行
if __name__ == "__main__":
    # 解析原子信息
    residues, residue_names = atoms_in_residue(atomlist, len(set(atom.resnr for atom in atomlist)))
    
    # 生成环约束（键角、二面角、对角线距离限制）
    ring_angles = generate_ring_angles(atomlist, residues, residue_names, force_constant=400)
    ring_dihedrals = generate_ring_dihedrals(atomlist, residues, residue_names, force_constant=-50)
    ring_diagonals = generate_ring_diagonals(atomlist, residues, residue_names)

    alllines = get_itp_lines(ring_angles, ring_dihedrals, ring_diagonals)

    f = open(input_topology, 'r')
    for rl in f: 
        if rl.startswith('; Include Position'):
            break 
        if rl.endswith('\n'):
            rl = rl[:-1]
        print (rl)
    for line in alllines:
        print (line)
    for rl in f:
        if rl.endswith('\n'):
            rl = rl[:-1]
        print (rl)
