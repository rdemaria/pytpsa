#Copyright (c) 2008, Riccardo De Maria
#All rights reserved.

from polmap import *
from pol import *


def Jmat(vars):
  """ Simplectic matrix or order 2 n
    >>> from lie import *
    >>> print Jmat('x px y py z pz'.split())
    pz=- 1.0*z
    px=- 1.0*x
    py=- 1.0*y
    y=py
    x=px
    z=pz
  """
  d=polmap()
  vars=list(vars)
  while vars:
    q,p=vars.pop(0),vars.pop(0)
    d[q]=pol(p)
    d[p]=-pol(q)
  return d

def parder(p,vars):
  """ Return derivative of the pol p given a vector of var
    >>> from lie import *
    >>> f=pol('l*(p**2+q**2)/2.')
    >>> print parder(f,'qp')
    q=q*l
    p=p*l
  """
  d=polmap()
  for i in vars:
    d[i]=p.der(i)
  return d

def lieop(f,vars):
  """ Return the lie operator for f
    >>> from lie import *
    >>> f=pol('l*(p**2+q**2)/2.')
    >>> Lf=lieop(f,'qp')
    >>> print Lf
    q=p*l
    p=- 1.0*q*l
    >>> print Jmat('qp')*parder(f,'qp')
    q=p*l
    p=- 1.0*q*l
    >>> print Lf*polmap(q='q',p=0)
    q=0
    p=- 1.0*q*l
    >>> print Lf*polmap(p='p',q=0)
    q=p*l
    p=0
  """
  d=polmap()
  vars=list(vars)
  while vars:
    q,p=vars.pop(0),vars.pop(0)
    d[q]=f.der(p)
    d[p]=-f.der(q)
  return d


def pbracket(f,g,vars):
  out=0
  vars=list(vars)
  while vars:
    q,p=vars.pop(0),vars.pop(0)
    out+=f.der(q)*g.der(p) - f.der(p)*g.der(q)
  return out


def pmbracket(f,g,vars):
  out=polmap()
  for i in g:
    out[i]=pbracket(f,g[i],vars)
  return out

def idmap(vars):
  vars=list(vars)
  g=polmap()
  for i in vars:
    g[i]=pol(i)
  return g

def liexp(f,vars):
  vars=list(vars)
  g=idmap(vars)
  out=polmap(g)
  out+=pmbracket(f,g,vars)
  out+=pmbracket(f,pmbracket(f,g,vars),vars)



