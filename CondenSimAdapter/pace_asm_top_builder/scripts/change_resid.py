import sys
# this code make sure all resid always start with 1
# resid between adjacent residues differ by 1
# also check if chain ID is not defined and provide 'A' if undefined 

f = open(sys.argv[1],'r')
#Read = True
rid_last = -1
rid_cnt = 0
for rl in f:
    if rl.startswith("ATOM"):
        Value = int(rl[22:26])-1
#        Read = False
        if rid_last!= Value:
         rid_cnt+=1
         rid_last = Value
    
    if rl.startswith("ATOM"):
     chainID = rl[21:22]
     if chainID==' ': chainID ='A'
    
    if rl.startswith("ATOM"):
        # Preserve the exact format: columns 0-20, then chain ID (1 char),
        # then residue number (4 chars), then continue from column 26
        # Original format: columns 0-20 (atom+serial+name+resname+chain), 22-26 (resid)
        # After change: columns 0-20, chain ID (1 char), resid (4 chars), rest from col 26
        print( rl[:21] + chainID + "%4d"%rid_cnt + rl[26:])
    else:
        print( rl[:-1])
f.close()
