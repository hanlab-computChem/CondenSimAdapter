import sys

f = open(sys.argv[1],'r')
last = sys.argv[2]

dih = []
cmap = []
ATM = ["N","CA","C","O1"]
for rl in f:
    if rl.startswith("ATOM"):
        sprl = rl.split()
        if rl[22:26].strip() == "1":
        #if sprl[5] == "1":
            if sprl[2] in ATM:
                dih.append(sprl[1])
        if rl[22:26].strip() == "2":
        #if sprl[5] == "2":
            if sprl[2] == "N":
                dih.append(sprl[1])
        if rl[22:26].strip() == "%d"%(int(last)-1):
        #if sprl[5] == "%d"%(int(last)-1):
            if sprl[2] == "C":
                cmap.append(sprl[1])
        if rl[22:26].strip() == last:
        #if sprl[5] == last:
            if sprl[2] in ATM:
                cmap.append(sprl[1])
f.close()

with open('dih_cmap.txt','w') as f:
    f.write("dihedral: %s\n"%(" ".join(dih)))
    f.write("cmap: %s\n"%(" ".join(cmap)))

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
