import sys

pep_name = sys.argv[2]
old_name = None
read_mol = False

for rl in open(sys.argv[1],'r'):
 if 'moleculetype' in rl:
  read_mol = True
#  print(rl,end='')
  continue

 if read_mol:
  if 'Protein' in rl:
   old_name = rl[:-1].split()[0]

  if '[' in rl:
   read_mol =False

if old_name is None:
 print('************')
 print('* ERROR    *')
 print('************')
 print('Incorrect molecular name in',sys.argv[1])
 sys.exit()

f = open(sys.argv[3],'w')
for rl in open(sys.argv[1],'r'):
 if old_name in rl:
  f.write(rl.replace(old_name, 'Protein_'+pep_name))
 elif 'posre.itp' in rl:
  f.write(rl.replace('posre.itp', 'posre_%s.itp'%pep_name))
 else:
  f.write(rl)



f.close()



