import sys

def _replace_line(_tuple_list, _replace_rule, rl):

  for _ri in _replace_rule:
   _rep_tuple = _ri[0]
   if _rep_tuple in _tuple_list:
    return _ri[1]


  return rl
   


fnm = sys.argv[1]

_read_data = False
_raw_data = []
for rl in open(fnm,'r'):
 if '[ atoms ]' in rl:
  _read_data = True
  continue

# reading atom information data 
 
 if _read_data:
  srl=rl[:-1].split()
  if len(srl)==0: continue   # remove blank line
  if ';' in srl[0]: continue # remove commented line  
  if '[' in srl[0]:
   _read_data = False
   continue

  _raw_data.append(srl)


# all atom info in _raw_data
# _raw_data[i] all text in each column for atom i+1


natom = len(_raw_data)
nres = int(_raw_data[-1][2])
 
print( '; There are',natom,'atoms in',nres,'residues')

# initiate index for all atoms
_CA = [-1]*nres
_CB = [-1]*nres
_H  = [-1]*nres
_N  = [-1]*nres
_C  = [-1]*nres
_O  = [-1]*nres


_fill_list = [ [_CA, 'CA'],[_CB, 'CB'],[_H,'H'],[_N,'N'],\
               [_C,'C'],[_O,'O']]

_atom_in_res=[[] for i in range(nres)]
for aid,i in enumerate(_raw_data):
 
 resid = int(i[2])-1
 atomnm = i[4]
 _atom_in_res[resid].append(aid)

 for _fi in _fill_list:
  if _fi[1]==atomnm:
   _fi[0][resid]= aid+1
   break

# _atom_in_res atom idx in each residue 
# till here all the code can be used for general purpose


# start for patching NH2 at C-terminus
# by replacing NMA to NH2

if _raw_data[-1][3]!='NMA':
 print( '******Error: there is no NMA residue at C-terminus.')
 sys.exit()


_atom_section = False
_excl_section = False
_pair_section = False
_add_pair_yet = False
_add_excl_yet = False

# DB: atom replacement
_a1 = _atom_in_res[-1][0]
_a2 = _atom_in_res[-1][1]
_a3 = _atom_in_res[-1][2]

_replace = [[[_a1+1],\
 '    %d     NT    %d  NH2   N     %s      0      10  ; N in NH2'%(_a1+1, nres, _raw_data[_a1][5])],\
[[_a2+1],\
 '    %d     HA    %d  NH2   H1     %s      0      5  ; N in NH2'%(_a2+1, nres, _raw_data[_a2][5])],\
[[_a3+1],\
 '    %d     HA    %d  NH2   H2     %s      0      5  ; N in NH2'%(_a3+1, nres, _raw_data[_a3][5])]\
]

_ar = _atom_in_res


# pairs to be replaced

if nres>=2:
 _idx=nres-1
 _replace_2=[[ [_CA[_idx-1], _H[_idx]],'%d  %d  1  0.0  0.0  ; NH2 replace'%(_CA[_idx-1], _H[_idx])  ],\
[ [_CA[_idx-1], _CA[_idx]],'%d  %d  1  8.70713E-04  1.89535E-06  ; NH2 replace'%(_CA[_idx-1], _CA[_idx])  ]\
]
else:
 _replace_2=[]

#print(_replace_2)

if nres>=2: 
 _idx = nres-1
 _add = [[[_N[_idx-1],_CA[_idx]],'%d  %d  1  3.89209E-03  8.47223E-06 ; NH2'%(_N[_idx-1],_CA[_idx])],\
[ [_H[_idx-1],_N[_idx]],'%d  %d  1  0.0   0.0 ; NH2'%(_H[_idx-1],_N[_idx])   ],\
[ [_H[_idx-1],_H[_idx]],'%d  %d  1  0.0   1.0E-07 ; NH2'%(_H[_idx-1],_H[_idx])   ],\
[ [_H[_idx-1],_CA[_idx]],'%d  %d  1  0.0   0.0 ;  NH2'%(_H[_idx-1],_CA[_idx])   ]\
]

if nres>=3:
 _idx = nres-1
 _add+= [[[_C[_idx-2],_H[_idx]],'%d  %d  1  0.0  0.0 ;  NH2'%(_C[_idx-2],_H[_idx])],\
[ [_O[_idx-2],_H[_idx]],'%d  %d  1  0.0   0.0 ;  NH2'%(_O[_idx-2],_H[_idx])  ],\
[ [_N[_idx-2],_CA[_idx]],'%d  %d  1  3.89209E-03  8.47223E-06 ;  NH2'%(_N[_idx-2],_CA[_idx])  ],\
[ [_H[_idx-2],_N[_idx]],'%d  %d  1  0.0   0.0 ;  NH2'%(_H[_idx-2],_N[_idx])  ],\
[ [_H[_idx-2],_H[_idx]],'%d  %d  1  0.0   1.0E-07 ;  NH2'%(_H[_idx-2],_H[_idx])  ],\
[ [_H[_idx-2],_CA[_idx]],'%d  %d  1  0.0   0.0 ;  NH2'%(_H[_idx-2],_CA[_idx])  ]\
]


if nres>=4:
 _idx = nres-1
 _add+= [[[_CA[_idx-3],_CA[_idx]],'%d  %d  1  6.33374E-03 2.22868E-05 ;  NH2'%(_CA[_idx-3],_CA[_idx])],\
[ [_C[_idx-3],_N[_idx]],'%d  %d  1  0.0  4.00293E-05 ;  NH2'%(_C[_idx-3],_N[_idx])  ],\
[ [_C[_idx-3],_H[_idx]],'%d  %d  1  0.0  0.0 ;  NH2'%(_C[_idx-3],_H[_idx])  ],\
[ [_C[_idx-3],_CA[_idx]],'%d  %d  1  1.34313E-02  1.73461E-05 ;  NH2'%(_C[_idx-3],_CA[_idx])  ],\
[ [_O[_idx-3],_H[_idx]],'%d  %d  1  0.0  0.0 ;  NH2'%(_O[_idx-3],_H[_idx])  ],\
[ [_O[_idx-3],_CA[_idx]],'%d  %d  1  0.0  6.36867E-06 ;  NH2'%(_O[_idx-3],_CA[_idx])  ]\
]


if nres>=5:
 _idx = nres-1
 _add+= [\
[ [_CA[_idx-4],_CA[_idx]],'%d  %d  1  6.33374E-03 2.22868E-05 ;  NH2'%(_CA[_idx-4],_CA[_idx])],\
[ [_C[_idx-4],_N[_idx]],'%d  %d  1  0.0  4.00293E-05 ;  NH2'%(_C[_idx-4],_N[_idx])  ],\
[ [_C[_idx-4],_H[_idx]],'%d  %d  1  0.0  0.0 ;  NH2'%(_C[_idx-4],_H[_idx])  ],\
[ [_C[_idx-4],_CA[_idx]],'%d  %d  1  1.34313E-02  1.73461E-05 ;  NH2'%(_C[_idx-4],_CA[_idx])  ],\
[ [_O[_idx-4],_H[_idx]],'%d  %d  1  0.0  0.0 ;  NH2'%(_O[_idx-4],_H[_idx])  ],\
[ [_O[_idx-4],_CA[_idx]],'%d  %d  1  0.0  6.36867E-06 ;  NH2'%(_O[_idx-4],_CA[_idx])  ]\
]


for rl in open(fnm, 'r'):

# modify atom section of top

 srl=rl[:-1].split()

 if len(srl)==0:
   print( rl,end='')
   continue

 if ';' in srl[0]:
   print( rl,end='')
   continue

 if '[' in srl[0]:
  _atom_section = False
  _excl_section = False
  _pair_section = False

  if '[ atoms ]' in rl:
   _atom_section = True
   print( rl,end='')
   continue

  if '[ exclusions ]' in rl:
   _excl_section = True
   print( rl,end='')
   continue

  if '[ pairs ]' in rl:
   _pair_section = True
   print( rl,end='')
   continue

  print( rl,end='')
  continue
  

 if _atom_section:


  __a = int(srl[0])
  print( _replace_line([[__a]],_replace,rl[:-1]))

 elif _excl_section:
  if not _add_excl_yet:
   _add_excl_yet = True
   for __i in _add+_replace_2:
    if -1 in __i[0]: continue
    print( __i[0][0],'  ',__i[0][1])

  print( rl,end='')
     
 elif _pair_section:

## add pairs at the beginning of the pair section
  if not _add_pair_yet:
   _add_pair_yet = True
   for __i in _add:
    if -1 in __i[0]: continue
    print( __i[1])


  __a = int(srl[0])
  __b = int(srl[1])
   # take into account the permutation of indices
  print( _replace_line([[__a, __b] , [__b, __a]], _replace_2, rl[:-1]))

 else:
  print( rl,end='')



