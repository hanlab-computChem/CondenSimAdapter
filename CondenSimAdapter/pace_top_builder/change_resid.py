import sys
# this code make sure all resid always start with 1
# resid between adjacent residues differ by 1

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
        print( rl[:22]+"%4d"%rid_cnt + rl[26:-1])
    else:
        print( rl[:-1])
f.close()
