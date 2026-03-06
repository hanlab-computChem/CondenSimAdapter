#!/usr/bin/env python3
import sys

f = open(sys.argv[1],'r')
last = sys.argv[2]

dih = []
cmap = []
ATM = ["N","CA","C","O"]

print(f"DEBUG: last = {last}")

for rl in f:
    if rl.startswith("ATOM"):
        sprl = rl.split()
        resnum = sprl[5]
        
        if resnum == '%d' % (int(last) - 1):
            if sprl[2] == "C":
                print(f"DEBUG: Found C in res {resnum}: atom {sprl[1]}")
                cmap.append(sprl[1])
        if resnum == last:
            if sprl[2] in ATM:
                print(f"DEBUG: Found {sprl[2]} in res {last}: atom {sprl[1]}")
                cmap.append(sprl[1])

print(f"\nDEBUG: cmap = {cmap}")
print(f"DEBUG: length = {len(cmap)}")

f.close()

# Now run the actual C-N-ter.py logic
f = open(sys.argv[3],'r')
change_N = False
change_C = False
for rl in f:
    if change_N and (sys.argv[4]=='Nter' or sys.argv[4]=='both'):
        sprl = rl.split()
        if set(sprl[:4]) == set(dih):
            if sprl[-1] == "1":
                sprl[6] = "4.0"
                print( "\t".join(sprl))
                continue

    if change_C and (sys.argv[4]=='Cter' or sys.argv[4]=='both'):
        sprl = rl.split()
        if len(sprl) == 0:
            print ("")
            print ("[ cmap ]")
            print ("\t".join(cmap),"\t1")

    if ("dihedral" in rl) and change_N :
        change_C = True
    
    if "; Include Position restraint file" in rl:
        change_C = False

    if "dihedral" in rl:
        if change_N == False:
            change_N = True
        else:
            change_N = False

    print( rl[:-1])
f.close()
