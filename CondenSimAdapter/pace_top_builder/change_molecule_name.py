import os
import sys

topol = sys.argv[1]
system_name = sys.argv[2]

lines = []

with open(topol, 'r') as f:
    mole_flag = False
    moleculetype_flag = False
    for rl in f:
        
        if rl.startswith('[ moleculetype'):
            moleculetype_flag = True
        if moleculetype_flag:
            srl = rl.split()
            if len(srl) == 2 and srl[1] =='3':
                moleculetype = srl[0]
                print (moleculetype)
                lines.append("{}      3 \n".format(system_name))
                moleculetype_flag = False
                continue
        if rl.startswith('[ molecules'):
            mole_flag = True
        if not mole_flag:
            lines.append(rl)
            continue
        if rl.startswith(moleculetype):
            srl = rl.split()
            monomer_number = srl[1]
            lines.append('{}        {}\n'.format(system_name, monomer_number))
            continue
        lines.append(rl)

with open(topol, 'w') as f:
    for rl in lines:
        f.write(rl)
    