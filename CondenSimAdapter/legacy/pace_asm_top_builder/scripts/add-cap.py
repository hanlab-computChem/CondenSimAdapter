import math
import sys


def get_xp(atom):
 return [atom[5],atom[6],atom[7]]

def pdb_line(rl,autoAssign=False, ass_idx=0):
    """
    Read a line and parse it into PDB components
    [0]: atom number
    [1]: atom name
    [2]: residu name
    [3]: chain name
    [4]: resid
    [5-7]: x,y,z
    [8]: segname
    """
    atominfo=[]
#atom indice in PDB will not be used if autoAssign is on
#instead, all atom indices start at 0
    if (rl[0:4]!='ATOM'):
        return 'None'
    if not autoAssign:
     atominfo.append(int(rl[6:11]))
    else:
     atominfo.append(ass_idx)

    atominfo.append(rl[12:17].strip())
    atominfo.append(rl[17:21].strip())
    atominfo.append(rl[21:22])
    atominfo.append(int(rl[22:26]))
    atominfo.append(float(rl[30:38]))
    atominfo.append(float(rl[38:46]))
    atominfo.append(float(rl[46:54]))
    if len(rl)>=75:
     kl=len(rl)
     if kl>76: kl=76
     atominfo.append(rl[72:kl].strip())
    else:
     atominfo.append('')
    return atominfo

def print_pdb_line(aid, anm, rnm, cnm, rid, x, y, z, snm):
 if len(anm)<3:
  anm1=anm+' '
 else:
  anm1=anm
 if len(rnm)<4:
  rnm1=rnm+' '
 else:
  rnm1=rnm
 if aid>99999: aid=99999
 print( "ATOM  %5d%5s %4s%1s"%(aid, anm1,rnm1,cnm)+\
       '%4d    %8.3f%8.3f%8.3f'%(rid,x,y,z)+\
       '%18s%3s'%(' ',snm))

def print_pdb_atom(atom):
 print_pdb_line(atom[0],atom[1],atom[2],atom[3],atom[4],\
                atom[5],atom[6],atom[7],atom[8])


def _fsgn(a):
 if a>0: return 1
 if a<0: return -1
 return 0

def VECADD(c,a,b):
 c[0]=a[0]+b[0]
 c[1]=a[1]+b[1]
 c[2]=a[2]+b[2]

def VECSUB(c,a,b):
 c[0]=a[0]-b[0]
 c[1]=a[1]-b[1]
 c[2]=a[2]-b[2]

def VECLEN(a):
 return (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])**0.5

#define VECDOT(a,b)   ((a)[XX]*(b)[XX]+(a)[YY]*(b)[YY]+(a)[ZZ]*(b)[ZZ])

def VECDOT(a,b):
 return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def VECOUT(c,a,b):
  c[0]=a[1]*b[2]-a[2]*b[1]
  c[1]=a[2]*b[0]-a[0]*b[2]
  c[2]=a[0]*b[1]-a[1]*b[0]

def VECELONG(c,a,v):
  c[0]=a[0]*(v)
  c[1]=a[1]*(v)
  c[2]=a[2]*(v)



def cal_xyz( c1, c2, c3,  b,  a,  d):
    ci=[0.,0.,0.]
    r=[0.,0.,0.]
    h=[0.,0.,0.]
    r21=[0.,0.,0.]
    r23=[0.,0.,0.]
    rb=[0.,0.,0.]
    rb21=[0.,0.,0.]
    rc=[0.,0.,0.]




    PI=3.141592654
    ca=a
    cb=b
    cd=d
    ca=ca/180.*PI
    cd=cd/180.*PI

    VECSUB(r21, c1, c2)
    VECSUB(r23, c3, c2)

    VECOUT(r, r21, r23)

    VECOUT(h, r21, r)



    rv=VECLEN(r)
    hv=VECLEN(h)
    r21v=VECLEN(r21)


    VECELONG(r,r,1./rv)

    VECELONG(h,h,1./hv)

    VECELONG(r21,r21,1./r21v)


    rbv=cb*math.sin(ca)

    VECELONG(rb21, r21, math.cos(PI-ca)*cb)

    VECELONG(rc, h, math.cos(PI-math.fabs(cd)) )
    VECELONG(rb, r, math.sin(math.fabs(cd))*_fsgn(cd) )
    VECADD (rb, rb, rc)
    VECELONG (rb, rb, rbv)


    VECADD( ci, rb, rb21)
    VECADD( ci, ci, c1)

    return ci

#main#

atoms=[]
nm=[]
xp=[]
aid=[]
rid=[]

for rl in open(sys.argv[1],'r'):
 if 'ATOM' in  rl:
  if rl[0:4]!='ATOM': continue
  atom = pdb_line(rl)
  atoms.append(atom)

  if sys.argv[2]=='Nter' or sys.argv[2]=='both':
   atom[4]+=1
   atom[0]+=3
  aid.append(atom[0])
  rid.append(atom[4])
  nm.append(atom[1])

#print nm

max_aid = max(aid)
max_rid = max(rid)
min_aid = min(aid)
min_rid = min(rid)

ft_C=-1
ft_CA=-1
ft_N=-1
_id=0
while ft_C<0 or ft_CA<0 or ft_N<0:
 if ft_C<0:
  if nm[_id]=="C":
   ft_C=_id
 if ft_CA<0:
  if nm[_id]=="CA":
   ft_CA=_id
 if ft_N<0:
  if nm[_id]=="N":
   ft_N=_id

 _id+=1


l_C=-1
l_CA=-1
l_O=-1
_id=len(nm)

while l_C<0 or l_CA<0 or l_O<0:
 _id-=1
 if l_C<0:
  if nm[_id]=="C":
   l_C=_id
 if l_CA<0:
  if nm[_id]=="CA":
   l_CA=_id
 if l_O<0:
  if nm[_id]=="O":
   l_O=_id


xp_f_N = get_xp(atoms[ft_N])
xp_f_C = get_xp(atoms[ft_C])
xp_f_CA = get_xp(atoms[ft_CA])
xp_l_C = get_xp(atoms[l_C])
xp_l_CA = get_xp(atoms[l_CA])
xp_l_O = get_xp(atoms[l_O])

xp_ace_C = cal_xyz(xp_f_N, xp_f_CA, xp_f_C, 1.33, 120.0, -60.0)
xp_ace_O = cal_xyz(xp_ace_C, xp_f_N, xp_f_CA, 1.2, 120.0, 0.0)
xp_ace_CA = cal_xyz(xp_ace_C, xp_f_N, xp_f_CA, 1.5, 120.0, 180.0)

xp_nma_N = cal_xyz(xp_l_C, xp_l_CA, xp_l_O, 1.3, 120.0, 180.0)
xp_nma_CA = cal_xyz(xp_nma_N, xp_l_C, xp_l_O, 1.5, 120.0, 0.0)
xp_nma_H = cal_xyz(xp_nma_N, xp_l_C, xp_l_CA, 1.0, 120.0, 0.0)


#print_pdb_line(aid, anm, rnm, cnm, rid, x, y, z, snm):
if sys.argv[2]=='Nter' or sys.argv[2]=='both':
 print_pdb_line(min_aid-3, 'CA', 'ACE', atoms[0][3], min_rid-1,\
               xp_ace_CA[0], xp_ace_CA[1], xp_ace_CA[2], atoms[0][8])
 print_pdb_line(min_aid-2, 'C', 'ACE', atoms[0][3], min_rid-1,\
               xp_ace_C[0], xp_ace_C[1], xp_ace_C[2], atoms[0][8])
 print_pdb_line(min_aid-1, 'O', 'ACE', atoms[0][3], min_rid-1,\
               xp_ace_O[0], xp_ace_O[1], xp_ace_O[2], atoms[0][8])

for _a in atoms: print_pdb_atom(_a)

if sys.argv[2]=='Cter' or sys.argv[2]=='both':
 print_pdb_line(max_aid+1, 'N', 'NMA', atoms[0][3], max_rid+1,\
              xp_nma_N[0], xp_nma_N[1], xp_nma_N[2], atoms[0][8])
 print_pdb_line(max_aid+2, 'CA', 'NMA', atoms[0][3], max_rid+1,\
              xp_nma_CA[0], xp_nma_CA[1], xp_nma_CA[2], atoms[0][8])
 print_pdb_line(max_aid+3, 'H', 'NMA', atoms[0][3], max_rid+1,\
              xp_nma_H[0], xp_nma_H[1], xp_nma_H[2], atoms[0][8])






