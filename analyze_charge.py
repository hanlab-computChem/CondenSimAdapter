#!/usr/bin/env python3
"""
计算蛋白质序列的电荷数
正电荷氨基酸: K (赖氨酸), R (精氨酸), H (组氨酸)
负电荷氨基酸: D (天冬氨酸), E (谷氨酸)
"""

from collections import OrderedDict

def extract_sequence_from_pdb(pdb_file):
    """从 PDB 文件中提取氨基酸序列"""
    residues = OrderedDict()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # PDB 格式: 
                # COLUMNS        DATA  TYPE    FIELD        DEFINITION
                # 18-20          Atom   resName   Residue name
                # 22-27          Atom   resSeq    Residue sequence number
                res_name = line[17:20].strip()
                res_seq = int(line[22:26].strip())
                
                # 只保留 CA (alpha carbon) 原子来代表残基
                atom_name = line[12:16].strip()
                if atom_name == 'CA' and res_seq not in residues:
                    residues[res_seq] = res_name
    
    # 将三字母代码转换为一字母代码
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XAA': 'X',
        'UNK': 'X'
    }
    
    sequence = ''.join(three_to_one.get(res, 'X') for res in residues.values())
    return sequence

def calculate_charge(sequence):
    """计算序列的正负电荷数"""
    # 带正电的氨基酸 (pH 7)
    positive = {'K': 1, 'R': 1, 'H': 0.5}  # His 在 pH 7 只带部分正电
    # 带负电的氨基酸 (pH 7)
    negative = {'D': -1, 'E': -1}
    
    pos_charge = sum(positive.get(aa, 0) for aa in sequence)
    neg_charge = sum(negative.get(aa, 0) for aa in sequence)
    
    return pos_charge, neg_charge

def main():
    # 文件路径
    h1_pdb = 'Tutorials/Tutorial_2_H1_ProTalpha_droplet/input/H1.pdb'
    prota_fasta = 'Tutorials/Tutorial_2_H1_ProTalpha_droplet/input/prota.fasta'
    
    print("=" * 60)
    print("蛋白质电荷分析")
    print("=" * 60)
    
    # 提取 H1 序列
    print("\n正在从 PDB 文件提取 H1 序列...")
    h1_sequence = extract_sequence_from_pdb(h1_pdb)
    print(f"H1 序列长度: {len(h1_sequence)} 氨基酸")
    print(f"H1 序列: {h1_sequence}")
    
    # 读取 prota 序列
    print("\n正在读取 ProtAlpha 序列...")
    with open(prota_fasta, 'r') as f:
        lines = f.readlines()
        prota_sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    print(f"ProtAlpha 序列长度: {len(prota_sequence)} 氨基酸")
    print(f"ProtAlpha 序列: {prota_sequence}")
    
    # 计算电荷
    print("\n" + "=" * 60)
    print("电荷计算结果")
    print("=" * 60)
    
    # H1 电荷
    pos_h1, neg_h1 = calculate_charge(h1_sequence)
    net_h1 = pos_h1 + neg_h1
    
    print(f"\n【H1 蛋白】")
    print(f"  正电荷 (K, R, H): {pos_h1:.1f}")
    print(f"  负电荷 (D, E):    {neg_h1:.1f}")
    print(f"  净电荷:           {net_h1:.1f}")
    
    # 统计各氨基酸数量
    aa_counts = {aa: h1_sequence.count(aa) for aa in 'KREDCH'}
    print(f"  各带电氨基酸数量:")
    print(f"    K (赖氨酸): {aa_counts['K']}")
    print(f"    R (精氨酸): {aa_counts['R']}")
    print(f"    H (组氨酸): {aa_counts['H']}")
    print(f"    D (天冬氨酸): {aa_counts['D']}")
    print(f"    E (谷氨酸): {aa_counts['E']}")
    
    # ProtAlpha 电荷
    pos_prota, neg_prota = calculate_charge(prota_sequence)
    net_prota = pos_prota + neg_prota
    
    print(f"\n【ProtAlpha 蛋白】")
    print(f"  正电荷 (K, R, H): {pos_prota:.1f}")
    print(f"  负电荷 (D, E):    {neg_prota:.1f}")
    print(f"  净电荷:           {net_prota:.1f}")
    
    # 统计各氨基酸数量
    aa_counts_prota = {aa: prota_sequence.count(aa) for aa in 'KREDCH'}
    print(f"  各带电氨基酸数量:")
    print(f"    K (赖氨酸): {aa_counts_prota['K']}")
    print(f"    R (精氨酸): {aa_counts_prota['R']}")
    print(f"    H (组氨酸): {aa_counts_prota['H']}")
    print(f"    D (天冬氨酸): {aa_counts_prota['D']}")
    print(f"    E (谷氨酸): {aa_counts_prota['E']}")
    
    # 总结
    print("\n" + "=" * 60)
    print("总结")
    print("=" * 60)
    print(f"H1 蛋白带有 {pos_h1:.0f} 个正电荷")
    print(f"ProtAlpha 蛋白带有 {abs(neg_prota):.0f} 个负电荷")
    print(f"\n电荷比: H1(+) : ProtAlpha(-) = {pos_h1:.0f} : {abs(neg_prota):.0f}")
    if pos_h1 > 0 and neg_prota < 0:
        ratio = abs(neg_prota) / pos_h1 if pos_h1 > 0 else float('inf')
        print(f"大约 {ratio:.2f} 个 ProtAlpha 负电荷 对应 1 个 H1 正电荷")

if __name__ == '__main__':
    main()
