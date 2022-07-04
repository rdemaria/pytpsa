from numpy import matrix, r_, c_
from .dop import *

def jacobian(self):
  vars=self.keys()
  o= [ dop(i)(self) for i in vars ]
  o= [ [ c[j] for j in vars ] for c in o ]
  return matrix(o)


def frommatrix(self,vars):
  pvars=[pol(i) for i in vars]
  m=polmap()
  for ni,i in enumerate(vars):
    m[i]=0
    for nj,j in enumerate(vars):
      m[i]+=self[ni,nj]*pvars[nj]
  return m

def tomatrix(self,vars=None):
  if not vars:
    vars=self[self.keys()[0]].vars
  out=[]
  for vi in vars:
    p=self[vi]
    oi=[]
    for vj in vars:
      oi.append(p.getlcoef(vj))
    out.append(oi)
  return matrix(out)


def linearpart(self):
  vars=self[self.keys()[0]].vars
  m=tomatrix(self)
  return frommatrix(m,vars)


def linv(self):
  vars=self[self.keys()[0]].vars
  m=tomatrix(self,vars)
  return frommatrix(m**-1,vars)

def lin(self):
  vars=self[self.keys()[0]].vars
  m=tomatrix(self,vars)
  return frommatrix(m,vars)

j2=matrix([[0.,1],[-1,0]])
n2=matrix([[0.,0],[0,0]])
i2=matrix([[1,0],[0,1]])
j4=r_[c_[j2,n2],c_[n2,j2]]
j6=r_[c_[j2,n2,n2],c_[n2,j2,n2],c_[n2,n2,j2]]


def normalfactor(eig,i,j):
  # Prod_k lambda[k]**j[k] - lambda[i]
  res=1
  for k,jk in enumerate(j):
    res*=eig[k]**jk
  res-=eig[i]
  return res

