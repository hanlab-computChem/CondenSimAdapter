import os
import sys

top_name = sys.argv[1]
posre_name = sys.argv[2]

'''#ifdef POSRES'''
alllines = []
with open(top_name, 'r') as f:
    for line in f:
        if line.startswith('#ifdef POSRES') and not line.startswith('#ifdef POSRES_WATER'):
            alllines.append(line)
            with open(posre_name, 'r') as posre_f:
                for line in posre_f:
                    alllines.append(line)
        elif line.startswith('#include "posre.itp"'):
            continue
        else:
            alllines.append(line)

with open(top_name, 'w') as f:
    for line in alllines:
        f.write(line)