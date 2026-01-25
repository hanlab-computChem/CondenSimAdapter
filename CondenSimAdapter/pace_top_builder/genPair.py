import sys

gromacs_top = sys.argv[1]
#gromacs_top = 'draft.top'

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

def find_atom(residue, atn):
    # find an atom by name in a certain residue
    for atom in residue:
        if atom.atom_name == atn:
            return atom.nr 
    return 0

def find_O_atom(residue, atn):
    # find an atom by name in a certain residue
    for atom in residue:
        if atom.atom_name == atn:
            return atom.nr 
        elif atom.atom_name == 'OT1' and atom.atom_type == 'OCT':
            return atom.nr
    return 0


def find_CA_atom(residue, atn):
    # find an atom by name in a certain residue
    for atom in residue:
        if atom.atom_name == atn:
            return atom.nr 
        elif atom.residue == 'FMO' and atom.atom_name == 'OM':
            return atom.nr
    return 0


# initialize atom information from Gromacs Topology
f = open(gromacs_top, 'r')
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

def recount_residue(atomlist):
    residue_number = 1 
    last_residue_number = atomlist[0].resnr
    for atom in atomlist:
        if atom.resnr != last_residue_number:
            residue_number += 1
            last_residue_number = atom.resnr
        
        atom.resnr = residue_number

    return atomlist

atom_list = recount_residue(atomlist)
for atom in atom_list:
    print (atom.resnr, atom.atom_name,  atom.residue)


nres = atomlist[-1].resnr 

atoms_in_each_residues, residues = atoms_in_residue(atomlist, nres=nres) 
residues.append('nothing')


#initialize atomName list
CA = [-1 for i in range(nres + 1)]
C = [-1 for i in range(nres + 1)]
O = [-1 for i in range(nres + 1)]
NH = [-1 for i in range(nres + 1)]
H = [-1 for i in range(nres + 1)]
CB = [-1 for i in range(nres + 1)]
CG = [-1 for i in range(nres + 1)]
CG1 = [-1 for i in range(nres + 1)]
CG2 = [-1 for i in range(nres + 1)]
OG = [-1 for i in range(nres + 1)]
OE = [-1 for i in range(nres + 1)]
NE = [-1 for i in range(nres + 1)]
CD1 = [-1 for i in range(nres + 1)]
CD2 = [-1 for i in range(nres + 1)]
AUX = [-1 for i in range(nres + 1)]
OG1 = [-1 for i in range(nres + 1)]
OD1 = [-1 for i in range(nres + 1)]
OD2 = [-1 for i in range(nres + 1)]
HD1 = [-1 for i in range(nres + 1)]
ND1 =  [-1 for i in range(nres + 1)]
ND2 =  [-1 for i in range(nres + 1)]
HD21 = [-1 for i in range(nres + 1)]
HD22 = [-1 for i in range(nres + 1)]
HG = [-1 for i in range(nres + 1)]
HG1 = [-1 for i in range(nres + 1)]

for i in range(0, nres ):
    CA[i] = find_CA_atom(atoms_in_each_residues[i], 'CA')
    C[i] = find_atom(atoms_in_each_residues[i], 'C')
    O[i] = find_O_atom(atoms_in_each_residues[i], 'O')
    NH[i] = find_atom(atoms_in_each_residues[i], 'N')
    H[i] = find_atom(atoms_in_each_residues[i],'H')
    CB[i] = find_atom(atoms_in_each_residues[i],'CB')
    CG[i] = find_atom(atoms_in_each_residues[i],'CG')
    CG1[i] = find_atom(atoms_in_each_residues[i],'CG1')
    CG2[i] = find_atom(atoms_in_each_residues[i],'CG2')
    OG[i] = find_atom(atoms_in_each_residues[i],'OG')
    OE[i] = find_atom(atoms_in_each_residues[i],'OE')
    NE[i] = find_atom(atoms_in_each_residues[i],'NE2')
    CD1[i] = find_atom(atoms_in_each_residues[i],'CD1')
    CD2[i] = find_atom(atoms_in_each_residues[i],'CD2')
    AUX[i] = find_atom(atoms_in_each_residues[i],'AUXL')
    OG1[i] = find_atom(atoms_in_each_residues[i],'OG1')
    OD1[i] = find_atom(atoms_in_each_residues[i],'OD1')   
    OD2[i] = find_atom(atoms_in_each_residues[i],'OD2')
    HD1[i] = find_atom(atoms_in_each_residues[i],'HD1')
    ND1[i] = find_atom(atoms_in_each_residues[i],'ND1')
    ND2[i] = find_atom(atoms_in_each_residues[i],'ND2')
    HD21[i] = find_atom(atoms_in_each_residues[i],'HD21')    
    HD22[i] = find_atom(atoms_in_each_residues[i],'HD22')   
    HG1[i] = find_atom(atoms_in_each_residues[i],'HG1')
    HG[i] = find_atom(atoms_in_each_residues[i],'HG')

    print ("{:s} CA:{:d} C:{:d} O:{:d} N:{:d} H:{:d} CB:{:d} CG1:{:d} CG2:{:d} CG:{:d} OG:{:d} OG1:{:d} OE:{:d} NE:{:d} CD1:{:d} CD2:{:d}OD1:{:d} OD2:{:d} HD1:{:d} ND1:{:d} ND2:{:d} HG:{:d} HG1:{:d}".format( residues[i], CA[i],\
        C[i], O[i],NH[i],H[i],CB[i],CG1[i],CG2[i],CG[i],OG[i],OG1[i],OE[i],NE[i],CD1[i],CD2[i],OD1[i],OD2[i],HD1[i],ND1[i],ND2[i],HG[i],HG1[i]))
    
print ("[ exclusions ]")
# strcmp() == 0 mean "==" strcmp() != 0 mean '!='
nres0 = nres
nres = nres 
for i in range(0, nres):
    # i is res_index
    if (i <= nres - 4) :
        if (CA[i] > 0) and (CB[i + 4]>0) :
            print ('{}    {}'.format(CA[i], CB[i+4]))
        if (CA[i] > 0) and (CA[i + 4]>0) and (residues[i+4] =='GLY'):
            print ("{}    {}   ;  0.671  0.50 GLY nm".format(CA[i], CA[i+4]))
    if (i <= nres - 3):
        if (residues[i+3] != 'PRO') and (residues[i+3] != 'DPR'):   
            if ((O[i] >0) and NH[i+3] > 0):
                print("{}    {}  ;  Oi-Ni+3 ".format(O[i], NH[i+3])) 
    if (i <= nres-4):
        if (residues[i+3] != 'PRO') and (residues[i+3] != 'DPR'):   
            if ((O[i] >0) and NH[i+4] > 0):
                print("{}    {}  ;  Oi-Ni+4 ".format(O[i], NH[i+4]))         

    #correct for gamma particles:
    if (OG[i] > 0) and (H[i] > 0):
        print('{}    {}'.format(OG[i], H[i]))
    if (OG[i] > 0) and (O[i] > 0):
        print('{}    {}'.format(OG[i], O[i]))
    if (OG[i] > 0) and (NH[i+1] > 0):
        print('{}    {}'.format(OG[i], NH[i+1])) 
    if (OG[i] > 0) and (H[i+1] > 0):
        print('{}    {}'.format(OG[i], H[i+1]))                
    if (i >= 1):
        if (OG[i] >0) and (O[i-1] > 0):
            print('{}    {}'.format(OG[i], O[i-1]))
    #HG exclusion
    if (HG[i] > 0) and (H[i] > 0):
        print('{}    {}'.format(HG[i], H[i]))
    if (HG[i] > 0) and (O[i] > 0):
        print('{}    {}'.format(HG[i], O[i]))
    if (HG[i] > 0) and (NH[i+1] > 0):
        print('{}    {}'.format(HG[i], NH[i+1])) 
    if (HG[i] > 0) and (H[i+1] > 0):
        print('{}    {}'.format(HG[i], H[i+1]))                
    if (i >= 1):
        if (HG[i] >0) and (O[i-1] > 0):
            print('{}    {}'.format(HG[i], O[i-1]))

    # og1 
    if ((OG1[i]>0) and (H[i]>0)):
        print('{}    {}'.format(OG1[i], H[i]))
    if ((OG1[i]>0) and (O[i] > 0)):
        print('{}    {}'.format(OG1[i], O[i])) 
    if ((OG1[i]>0) and (NH[i+1]>0)):
        print('{}    {}'.format(OG1[i], NH[i+1]))
    if ((OG1[i]>0) and (H[i+1]>0)):
        print('{}    {}'.format(OG1[i], H[i+1]))
    if i >= 1:
        if ((OG1[i] > 0) and (O[i-1]>0)):
            print('{}    {}'.format(OG1[i], O[i-1]))

    # HG1 exclusion
    if ((HG1[i]>0) and (H[i]>0)):
        print('{}    {}'.format(HG1[i], H[i]))
    if ((HG1[i]>0) and (O[i] > 0)):
        print('{}    {}'.format(HG1[i], O[i])) 
    if ((HG1[i]>0) and (NH[i+1]>0)):
        print('{}    {}'.format(HG1[i], NH[i+1]))
    if ((HG1[i]>0) and (H[i+1]>0)):
        print('{}    {}'.format(HG1[i], H[i+1]))
    if i >= 1:
        if ((HG1[i] > 0) and (O[i-1]>0)):
            print('{}    {}'.format(HG1[i], O[i-1]))
    # ASP and ASN sideChain
    if (residues[i] == 'ASP') or (residues[i] == 'ASPH') or (residues[i] == 'ASN'):
        # As residue is ASP ASPH ASN, CG[i] always > 0
        if ((CG[i]>0) and (H[i]>0)):
            print('{}    {} ;  CG-Hi asx'.format(CG[i], H[i]))
        if ((CG[i]>0) and (H[i+1] >0)):
            print ('{}    {}   ;  CG-Hi+1 asx'.format(CG[i], H[i+1]))
        if ((CG[i]>0) and (NH[i+1] >0)):
            print ('{}    {}    ;  CG-Ni+1 asx'.format(CG[i], NH[i+1]))
        if ((ND2[i]>0) and NH[i] > 0):
            print('{}    {}    ;  ND2-Ni+1 asx '.format(ND2[i], NH[i]))
        if ((ND2[i]>0) and NH[i+1] > 0):
            print('{}    {}    ;  ND2-Ni+1 asx '.format(ND2[i], NH[i+1]))
        if residues[i] == 'ASN':
            if ((HD21[i]>0) and (O[i]>0)):
                print("{}    {}    ;  HD21-Oi asx".format(HD21[i], O[i]))
            if ((OD1[i]>0) and (H[i] >0)):
                print("{}    {}    ;  OD1-Hi asx".format(OD1[i], H[i]))
            if ((OD1[i]>0) and (H[i+1] >0)):
                print("{}    {}    ;  OD1-Hi+1 asx".format(OD1[i], H[i+1]))
        if residues[i] == 'ASP' or residues == 'ASPH':
            if ((OD1[i] >0 ) and (H[i] >0)):
                print("{}    {} ;  OD1-Hi asx".format(OD1[i], H[i]))
            if ((OD1[i] >0 ) and (H[i+1] >0)):
                print("{}    {} ;  OD1-Hi+1 asx".format(OD1[i], H[i+1]))                     
            if ((OD2[i] >0 ) and (H[i] >0)):
                print("{}    {} ;  OD2-Hi asx".format(OD2[i], H[i]))
            if ((OD2[i] >0 ) and (H[i+1] >0)):
                print("{}    {} ;  OD2-Hi+1 asx".format(OD2[i], H[i+1]))
            if ((CG[i] > 0) and (O[i]>0)):
                print ('{}    {}   ;  CG-Oi asx'.format(CG[i], O[i]))
    #PHE, TYR, TRP, HIS sidechain
    # HIS[:3]
    elif residues[i] == 'PHE' or residues[i] == 'TYR' or residues[i] == 'TRP' or residues[i][:3] == 'HIS':
        if ((CG[i] > 0) and (H[i] > 0)):
            print ('{}    {}'.format(CG[i], H[i]))
        if ((CG[i] > 0) and (H[i+1] > 0)):
            print ('{}    {}'.format(CG[i], H[i+1]))
        if ((CG[i] > 0) and (NH[i+1] > 0)):
            print ('{}    {}'.format(CG[i], NH[i+1]))
        if ((CD1[i]>0) and (H[i]>0)):
            print ('{}    {}'.format(CD1[i], H[i]))
        if ((ND1[i]>0) and (H[i]>0)):
            print ('{}    {}'.format(ND1[i], H[i]))
        if ((CD2[i]>0) and (H[i]>0)):
            print ('{}    {}'.format(CD2[i], H[i]))
        if ((CD1[i]>0) and (NH[i]>0)):
            print ('{}    {}'.format(CD1[i], NH[i]))    
        if ((ND1[i]>0) and (NH[i]>0)):
            print ('{}    {}'.format(ND1[i], NH[i]))           
        if ((CD2[i]>0) and (NH[i]>0)):
            print ('{}    {}'.format(CD2[i], NH[i]))          
        if i >= 1:
            if ((CD1[i] > 0) and (O[i-1]>0)):
                print ('{}    {}'.format(CD1[i], O[i-1])) 
            if ((CD2[i] > 0) and (O[i-1]>0)):
                print ('{}    {}'.format(CD2[i], O[i-1]))
        if ((CD1[i]>0) and (NH[i+1]>0)):
            print ('{}    {}'.format(CD1[i], NH[i+1]))
        if ((ND1[i]>0) and (NH[i+1]>0)):  
            print ('{}    {}'.format(ND1[i], NH[i+1])) 
        if ((CD2[i]>0) and (NH[i+1]>0)):
            print ('{}    {}'.format(CD2[i], NH[i+1]))        
    else:
        if((CG[i]>0) and (H[i]>0)):
            print ('{}    {}'.format(CG[i], H[i]))
        if((CG[i]>0) and (H[i+1]>0)):
            print ('{}    {}'.format(CG[i], H[i+1]))
    if((CG1[i]>0) and (H[i]>0)):
        print ('{}    {}'.format(CG1[i], H[i]))
    if((CG1[i]>0) and (H[i+1]>0)):
        print ('{}    {}'.format(CG1[i], H[i+1]))
    if((CG2[i]>0) and (H[i]>0)):
        print ('{}    {}'.format(CG2[i], H[i]))
    if((CG2[i]>0) and (H[i+1]>0)):
        print ('{}    {}'.format(CG2[i], H[i+1]))
    
    if((NH[i]>0) and (C[i+1]>0)):
        print ('{}    {}'.format(NH[i], C[i+1]))
    if((O[i]>0) and (O[i+1]>0)):
        print ('{}    {}'.format(O[i], O[i+1]))

    if((NH[i]>0) and (H[i+1]>0)):
        print ('{}    {}'.format(NH[i], H[i+1]))
    if((NH[i]>0) and (O[i+1]>0)):
        print ('{}    {}'.format(NH[i], O[i+1]))
    if((CA[i]>0) and (O[i+1]>0)):
        print ('{}    {}'.format(CA[i], O[i+1]))
    if((C[i]>0) and (O[i+1]>0)):
        print ('{}    {}'.format(C[i], O[i+1]))
    if((O[i]>0) and (C[i+1]>0)):
        print ('{}    {}'.format(O[i], C[i+1]))
    if((CA[i]>0) and (C[i+1]>0)):
        print ('{}    {}'.format(CA[i], C[i+1]))

    if (i <= nres - 2):
        if ((NH[i]>0) and (C[i+2]>0)):
            print ("{}    {}".format(NH[i], C[i+2]))
        if ((NH[i]>0) and (O[i+2]>0)):
            print ("{}    {}".format(NH[i], O[i+2]))
        if ((CA[i]>0) and (O[i+2]>0)):
            print ("{}    {}".format(CA[i], O[i+2]))
        if ((C[i]>0) and (O[i+2]>0)):
            print ("{}    {}".format(C[i], O[i+2]))
        if ((C[i]>0) and (NH[i+2]>0)):
            print ("{}    {}".format(C[i], NH[i+2]))
        if ((O[i]>0) and (CA[i+2]>0)):
            print ("{}    {}".format(O[i], CA[i+2])) 
        if ((O[i]>0) and (NH[i+2]>0)):
            print ("{}    {}".format(O[i], NH[i+2]))
        if ((CA[i]>0) and (CA[i+2]>0)):
            print ("{}    {}".format(CA[i], CA[i+2]))  
        if ((CA[i]>0) and (C[i+2]>0)):
            print ("{}    {}".format(CA[i], C[i+2]))
        if ((C[i]>0) and (CA[i+2]>0)):
            print ("{}    {}".format(C[i], CA[i+2]))
        if ((C[i]>0) and (C[i+2]>0)):
            print ("{}    {}".format(C[i], C[i+2]))     
        if (i >= 1):
            if ((C[i-1]>0) and  (O[i+2]>0)):
                print ("{}    {}".format(C[i-1], O[i+2])) 

print ("[ pairs ]")
for i in range(0, nres):
    # Cai --> Cbi+4 0.671 0.47
    if (i <= nres - 4):
        if (CA[i] > 0) and (CB[i + 4]>0):
            print ('{}    {}  1  2.89314E-02 3.11858E-04  ;CAi-CBi+4  0.671  0.47 nm'.format(CA[i], CB[i+4]))
        if (CA[i] > 0) and (CA[i + 4]>0) and (residues[i+4] =='GLY'):
            print ("{}    {}  1  4.19375E-02  6.55273E-04  ;CAi CAi+4 GLY  0.671  0.50  nm".format(CA[i], CA[i+4]))    
    # refinforce Oi-Ni+3 HB for turn
    if residues[i] != 'FMO':
        if (i <= nres - 3):
            if (residues[i+3] != 'PRO') and (residues[i+3] != 'DPR'):   
                if ((O[i] >0) and NH[i+3] > 0):
                    print("{}    {}  1  8.427637e-03   1.610547e-06;  Oi-Ni+3 14.7 0.24nm".format(O[i], NH[i+3])) # 1.123685E-02    2.147396E-06
        if (i <= nres-4):
            if (residues[i+3] != 'PRO') and (residues[i+3] != 'DPR'):   
                if ((O[i] >0) and NH[i+4] > 0):
                    print("{}    {}  1  8.427637e-03   1.610547e-06;  Oi-Ni+4  14.7 0.24nm".format(O[i], NH[i+4]))       
    else:
        if (i <= nres - 3):
            if (residues[i+3] != 'PRO') and (residues[i+3] != 'DPR'):   
                if ((O[i] >0) and NH[i+3] > 0):
                    print("{}    {}  1  6.742110e-03   1.288438e-06 ;  Oi-Ni+3 8.8 0.24nm FMO weaker HB".format(O[i], NH[i+3])) 
        if (i <= nres-4):
            if (residues[i+3] != 'PRO') and (residues[i+3] != 'DPR'):    
                if ((O[i] >0) and NH[i+4] > 0):
                    print("{}    {}  1  6.742110e-03   1.288438e-06 ;  Oi-Ni+4  8.8 0.24nm FMO weaker HB".format(O[i], NH[i+4]))           
    # correct for gamma particles
    if (OG[i] > 0) and (H[i] > 0):
        print('{}    {}  1  1.13380E-03  1.28550E-07 ; OG-H 2.5 0.22 nm '.format(OG[i], H[i]))
    if (OG[i] > 0) and (O[i] > 0):
        print('{}    {}  1 7.74841E-04   1.50095E-06 ; OG-Oi  2.5  0.27 nm  rep'.format(OG[i], O[i]))
    if (OG[i] > 0) and (NH[i+1] > 0):
        print('{}    {}  1  3.83970E-03  4.12285E-06 ;   OGi-Ni+1  0.894 0.32nm'.format(OG[i], NH[i+1])) 
    if (OG[i] > 0) and (H[i+1] > 0):
        print('{}    {}  1  1.13380E-03  1.28550E-07  ;  OGi-Hi+1 2.5  0.22nm'.format(OG[i], H[i+1]))                
    if (i >= 1):
        if (OG[i] >0) and (O[i-1] > 0):
            print('{}    {}  1  3.87420E-03   1.50095E-06 ;  OGi-Oi-1 2.5  0.27nm'.format(OG[i], O[i-1]))
    # HG exclusion        
    if (HG[i] > 0) and (H[i] > 0):
        print('{}    {}  1     0     0 ; HG serine exclusion '.format(HG[i], H[i]))
    if (HG[i] > 0) and (O[i] > 0):
        print('{}    {}  1     0     0 ; HG serine exclusion'.format(HG[i], O[i]))
    if (HG[i] > 0) and (NH[i+1] > 0):
        print('{}    {}  1     0     0 ;   HG serine exclusion'.format(HG[i], NH[i+1])) 
    if (HG[i] > 0) and (H[i+1] > 0):
        print('{}    {}  1     0     0  ;  HG serine exclusion'.format(HG[i], H[i+1]))                
    if (i >= 1):
        if (HG[i] >0) and (O[i-1] > 0):
            print('{}    {}  1     0     0  ; HG serine exclusion'.format(HG[i], O[i-1]))
    #og1
    if ((OG1[i]>0) and (H[i]>0)):
        print('{}    {}  1  2.26760E-04  1.28550E-07  ; OG1-Hi  0.1 0.287685nm rep'.format(OG1[i], H[i]))
    if ((OG1[i]>0) and (O[i] > 0)):
        print('{}    {}  1  7.74841E-04   1.50095E-06 ; OG1-Oi 0.1 0.353068nm rep'.format(OG1[i], O[i])) 
    if ((OG1[i]>0) and (NH[i+1]>0)):
        print('{}    {}  1  3.83970E-03  4.12285E-06 ;  OG1-Ni+1  0.894 0.32nm'.format(OG1[i], NH[i+1]))
    if ((OG1[i]>0) and (H[i+1]>0)):
        print('{}    {}  1  1.13380E-03  1.28550E-07  ;  OG1-Hi+1  2.5 0.22nm'.format(OG1[i], H[i+1]))
    if i >= 1:
        if ((OG1[i] > 0) and (O[i-1]>0)):
            print('{}    {}  1  3.87420E-03   1.50095E-06 ;  OG1-Oi-1 2.5 0.27nm '.format(OG1[i], O[i-1]))
    #HG1
    if ((HG1[i]>0) and (H[i]>0)):
        print('{}    {}  1 0 0 ;  HG1 threonine exclusion'.format(HG1[i], H[i]))
    if ((HG1[i]>0) and (O[i] > 0)):
        print('{}    {}  1 0 0 ;  HG1 threonine exclusion'.format(HG1[i], O[i])) 
    if ((HG1[i]>0) and (NH[i+1]>0)):
        print('{}    {}  1 0 0 ;  HG1 threonine exclusion'.format(HG1[i], NH[i+1]))
    if ((HG1[i]>0) and (H[i+1]>0)):
        print('{}    {}  1 0 0 ;  HG1 threonine exclusion'.format(HG1[i], H[i+1]))
    if i >= 1:
        if ((HG1[i] > 0) and (O[i-1]>0)):
            print('{}    {}  1 0 0 ;  HG1 threonine exclusion'.format(HG1[i], O[i-1]))
    # ASP and ASN sideChain
    if (residues[i] == 'ASP') or (residues[i] == 'ASPH') or (residues[i] == 'ASN'):
        # As residue is ASP ASPH ASN, CG[i] always > 0
        if ((CG[i]>0) and (H[i]>0)):
            print('{}    {}  1  0.0     8.5E-07  ;  CG-Hi asx'.format(CG[i], H[i]))
        if ((CG[i]>0) and (H[i+1] >0)):
            print ('{}    {}  1  0.0     8.5E-07  ;  CG-Hi+1 asx'.format(CG[i], H[i+1]))
        if ((CG[i]>0) and (NH[i+1] >0)):
            print ('{}    {}  1  1.91985E-03   2.06142E-06  ;  CG-Ni+1  0.447  0.32nm asx'.format(CG[i], NH[i+1]))
        if ((ND2[i]>0) and NH[i] > 0):
            print('{}    {}  1  1.30345E-03   9.50217E-07  ;  ND2-Ni 0.447 0.30nm asx'.format(ND2[i], NH[i]))
        if ((ND2[i]>0) and NH[i+1] > 0):
            print('{}    {}  1  2.76211E-03   4.26692E-06  ;  ND2-Ni+1 0.447 0.34nm asx'.format(ND2[i], NH[i+1]))
        if residues[i] == 'ASN':
            if ((HD21[i]>0) and (O[i]>0)):
                print("{}    {}  1  1.275068E-03  2.139210E-08  ;  HD21-Oi 19.0 0.16nm asx".format(HD21[i], O[i]))
            if ((OD1[i]>0) and (H[i] >0)):
                print("{}    {}   1  1.275068E-03  2.139210E-08  ;  OD1-Hi 19.0 0.16nm asx".format(OD1[i], H[i]))
            if ((OD1[i]>0) and (H[i+1] >0)):
                print("{}    {}   1  1.275068E-03  2.139210E-08 ;  OD1-Hi+1 19.0 0.16nm asx".format(OD1[i], H[i+1]))
        if residues[i] == 'ASP' or residues == 'ASPH' or residues == 'ASH':
            if ((OD1[i] >0 ) and (H[i] >0)):
                print("{}    {}   1  1.275068E-03  2.139210E-08 ;  OD1-Hi 19.0 0.16nm asx".format(OD1[i], H[i]))
            if ((OD1[i] >0 ) and (H[i+1] >0)):
                print("{}    {}   1  1.275068E-03  2.139210E-08 ;  OD1-Hi+1 19.0 0.16nm asx".format(OD1[i], H[i+1]))                     
            if ((OD2[i] >0 ) and (H[i] >0)):
                print("{}    {}   1  1.275068E-03  2.139210E-08 ;  OD2-Hi 19.0 0.16nm asx".format(OD2[i], H[i]))
            if ((OD2[i] >0 ) and (H[i+1] >0)):
                print("{}    {}   1  1.275068E-03  2.139210E-08 ;  OD2-Hi+1 19.0 0.16nm asx".format(OD2[i], H[i+1]))
            if ((CG[i] > 0) and (O[i]>0)):
                print ('{}    {}   1  7.41398E-03   2.29029E-06 ; 6.0  0.26nm CG-Oi asx'.format(CG[i], O[i]))
    #PHE, TYR, TRP, HIS sidechain
    # HIS[:3]
    elif (residues[i] == 'PHE') or (residues[i] == 'TYR') or (residues[i] == 'TRP') or (residues[i][:3] == 'HIS'):
        
        if ((CG[i] > 0) and (H[i] > 0)):
            print ('{}    {}  1  6.53184E-04  4.76171E-07  ;  CG-Hi 0.224 0.30nm aro'.format(CG[i], H[i]))
        if ((CG[i] > 0) and (H[i+1] > 0)):
            print ('{}    {}  1  6.53184E-04  4.76171E-07  ;  CG-Hi+1 0.224 0.30nm aro'.format(CG[i], H[i+1]))
        if ((CG[i] > 0) and (NH[i+1] > 0)):
            print ('{}    {}  1  1.91985E-03   2.06142E-06  ;  CG-Ni+1 0.447 0.32nm aro '.format(CG[i], NH[i+1]))
        if ((CD1[i]>0) and (H[i]>0)):
            print ('{}    {} 1  6.53184E-04  4.76171E-07  ;  Hi-CD1 0.224 0.30nm aro'.format(CD1[i], H[i]))
        if ((ND1[i]>0) and (H[i]>0)):
            print ('{}    {} 1  6.53184E-04  4.76171E-07  ;  Hi-ND1 0.224 0.30nm aro'.format(ND1[i], H[i]))
        if ((CD2[i]>0) and (H[i]>0)):
            print ('{}    {} 1  6.53184E-04  4.76171E-07  ;  Hi-CD2 0.224 0.30nm aro'.format(CD2[i], H[i]))
        if ((CD1[i]>0) and (NH[i]>0)):
            print ('{}    {}   1  1.30345E-03   9.50217E-07  ;  CD1-Ni 0.447 0.30nm aro'.format(CD1[i], NH[i]))    
        if ((ND1[i]>0) and (NH[i]>0)):
            print ('{}    {}   1  1.30345E-03   9.50217E-07  ;  ND1-Ni 0.447 0.30nm aro'.format(ND1[i], NH[i]))           
        if ((CD2[i]>0) and (NH[i]>0)):
            print ('{}    {}   1  1.30345E-03   9.50217E-07  ;  CD2-Ni 0.447 0.30nm aro'.format(CD2[i], NH[i]))          
        if i >= 1:
            if ((CD1[i] > 0) and (O[i-1]>0)):
                print ('{}    {}   1     1.06354E-03  6.32621E-07  ;   CD1-Oi-1 0.447 0.29nm  aro'.format(CD1[i], O[i-1])) 
            if ((CD2[i] > 0) and (O[i-1]>0)):
                print ('{}    {}   1     1.06354E-03  6.32621E-07  ;   CD2-Oi-1 0.447 0.29nm aro'.format(CD2[i], O[i-1]))
        if ((CD1[i]>0) and (NH[i+1]>0)):
            print ('{}    {}   1  2.76211E-03  4.26692E-06  ;  CD1-Ni+1 0.447 0.34nm aro'.format(CD1[i], NH[i+1]))
        if ((ND1[i]>0) and (NH[i+1]>0)):  
            print ('{}    {}   1  2.76211E-03  4.26692E-06  ;  ND1-Ni+1 0.447 0.34nm aro'.format(ND1[i], NH[i+1])) 
        if ((CD2[i]>0) and (NH[i+1]>0)):
            print ('{}    {}   1  2.76211E-03  4.26692E-06  ;  CD2-Ni+1 0.447 0.34nm aro'.format(CD2[i], NH[i+1])) 
    else:
        if((CG[i]>0) and (H[i]>0)):
            print ('{}    {}  1    1.74675E-03   1.70645E-06      ; Hi-CGi  0.2 0.24nm rep'.format(CG[i], H[i]))
        if((CG[i]>0) and (H[i+1]>0)):
            print ('{}    {}  1    1.74675E-03   1.70645E-06      ; Hi+1-CGi  0.2 0.24nm rep'.format(CG[i], H[i+1]))
    
    if((CG1[i]>0) and (H[i]>0)):
        print ('{}    {}  1  1.06354E-03  6.32621E-07      ; Hi-CGi H 0.2 0.24nm rep'.format(CG1[i], H[i]))
    if((CG1[i]>0) and (H[i+1]>0)):
        print ('{}    {}  1  1.06354E-03  6.32621E-07      ; Hi+1-CGi H 0.2 0.24nm rep'.format(CG1[i], H[i+1]))
    if((CG2[i]>0) and (H[i]>0)):
        print ('{}    {}  1  1.06354E-03  6.32621E-07      ; Hi-CGi H 0.2 0.24nm rep'.format(CG2[i], H[i]))
    if((CG2[i]>0) and (H[i+1]>0)):
        print ('{}    {}  1  1.06354E-03  6.32621E-07      ; Hi+1-CGi H 0.2 0.24nm rep'.format(CG2[i], H[i+1]))

    if((NH[i]>0) and (C[i+1]>0)):
        print ('{}    {}  1   3.37244E-03   5.68668E-06 ; Ni-Ci+1 0.5 0.345nm'.format(NH[i], C[i+1]))
    if((O[i]>0) and (O[i+1]>0) and residues[i] != 'FMO'):
        print ('{}    {}  1   4.13270E-03   5.33725E-06 ; Oi-Oi+1 0.8 0.33nm'.format(O[i], O[i+1]))

    if((NH[i]>0) and (H[i+1]>0)):
        print ('{}    {}  1   4.42E-04   1.64E-08  ;   Ni-Hi+1 2.98 0.1826nm'.format(NH[i], H[i+1]))
    if((NH[i]>0) and (O[i+1]>0)):
        print ('{}    {}  1   2.87870E-03  2.31737E-06 ; Ni-Oi+1 0.894 0.305nm'.format(NH[i], O[i+1]))
    if((CA[i]>0) and (O[i+1]>0) and residues[i]!= 'FMO'):
        print ('{}    {}  1   2.26146E-03   3.19637E-06  ; CAi-Oi+1 0.6  0.335nm'.format(CA[i], O[i+1]))
    if((C[i]>0) and (O[i+1]>0) and residues[i] != 'FMO'):
        print ('{}    {}  1    1.91985E-03   2.06142E-06  ; Ci-Oi+1 0.6  0.32nm'.format(C[i], O[i+1]))
    if((O[i]>0) and (C[i+1]>0) and residues[i] != 'FMO'):
        print ('{}    {}  1    9.29809E-04   3.60227E-07  ; Oi-Ci+1 0.6 0.27nm'.format(O[i], C[i+1]))
    if((CA[i]>0) and (C[i+1]>0) and residues[i] != 'FMO'):
        print ('{}    {}  1    2.49170E-03   6.92920E-06  ; CAi-Ci+1  0.45  0.375nm'.format(CA[i], C[i+1]))
    if((CB[i]>0) and (H[i+1]>0)):
        print ("{}    {}  1    1.74675E-03   1.70645E-06  ; CB-Hi+1  0.2 0.24 nm rep".format( CB[i],  H[i+1]))
    
    if (i <= nres - 2):
        if ((NH[i]>0) and (C[i+2]>0)):
            print ("{}    {}  1    3.37244E-03   5.68668E-06 ; Ni-Ci+2 0.671 0.345nm".format(NH[i], C[i+2]))
        if ((NH[i]>0) and (O[i+2]>0)):
            print ("{}    {}  1    2.87870E-03  2.31737E-06 ; Ni-Oi+2 0.894 0.305nm".format(NH[i], O[i+2]))
        if ((CA[i]>0) and (O[i+2]>0)):
            print ("{}    {}  1    2.26146E-03   3.19637E-06   ; CAi-Oi+2 0.6  0.335nm".format(CA[i], O[i+2]))
        if ((C[i]>0) and (O[i+2]>0)):
            print ("{}    {}  1    1.91985E-03   2.06142E-06   ; Ci-Oi+2 0.6  0.32nm ".format(C[i], O[i+2]))
        if ((C[i]>0) and (NH[i+2]>0)):
            print ("{}    {}  1     3.37244E-03   5.68668E-06 ; Ci-NHi+2 0.671 0.345nm".format(C[i], NH[i+2]))
        if ((O[i]>0) and (CA[i+2]>0)):
            print ("{}    {}  1    2.26146E-03   3.19637E-06  ; Oi-CAi+2 0.6  0.335nm".format(O[i], CA[i+2])) 
        
        if ((O[i]>0) and (NH[i+2]>0)):
            print ("{}    {}  1    1.10468E-03   3.41254E-07 ; Oi-Ni+2 0.894 0.26nm".format(O[i], NH[i+2]))
        if ((CA[i]>0) and (CA[i+2]>0)):
            print ("{}    {}  1    2.81500E-03   9.90525E-06 ; CAi-CAi+2 0.45  0.39nm".format(CA[i], CA[i+2]))  
        if ((CA[i]>0) and (C[i+2]>0)):
            print ("{}    {}  1    2.49170E-03   6.92920E-06  ; CAi-Ci+2 0.45  0.375nm".format(CA[i], C[i+2]))
        if ((C[i]>0) and (CA[i+2]>0)):
            print ("{}    {}  1    2.49170E-03   6.92920E-06  ; Ci-CAi+2 0.45  0.375nm".format(C[i], CA[i+2]))
        if ((C[i]>0) and (C[i+2]>0)):
            print ("{}    {}  1    2.17678E-03   4.73838E-06  ; Ci-Ci+2 0.45  0.36nm".format(C[i], C[i+2]))     
        if (i >= 1):
            if ((C[i-1]>0) and  (O[i+2]>0)):
                print ("{}    {}  1     1.91985E-03   2.06142E-06   ; 0.6  0.32 nm".format(C[i-1], O[i+2]))
    if (residues[i] == 'FMO'):
        if (O[i] > 0 and O[i+1] > 0):
            print ("{}    {}  1    2.10369E-02   1.10649E-05 ; FMO".format(O[i], O[i+1]))
        if (CA[i] > 0 and O[i+1] > 0):
            print ("{}    {}  1    1.03010E-07   1.69147E-15 ; FMO ".format(CA[i], O[i+1]))
        if (C[i] > 0 and O[i+1] > 0):
            print ("{}    {}  1    4.20128E-04   7.36847E-06 ; FMO".format(C[i], O[i+1]))
        if (O[i] > 0 and C[i+1] > 0):
            print ("{}    {}  1     7.83407E-04   8.59933E-07 ; FMO".format(O[i], C[i+1]))
        if (CA[i] > 0 and C[i+1] > 0):
            print ("{}    {}  1    1.42390E-02   2.50539E-04 ; FMO".format(CA[i], C[i+1]))
      

print ("[ impropers ]") 
for i in range(0, nres):
    if ((C[i]>0) and (O[i]>0) and (NH[i+1]>0) and (CA[i+1]>0)):
        print("{}  {}  {}  {}  2  0.0 75.0".format(O[i], C[i], NH[i+1], CA[i+1]))
    if ((C[i]>0) and (NH[i+1]>0)) and ((H[i+1]>0) and (CA[i+1]>0)):
        print("{} {} {} {} 2 0.0 150".format( NH[i+1], H[i+1], CA[i+1],C[i]))
f = open('res.temp','w')
f.write(str(nres))
f.close()