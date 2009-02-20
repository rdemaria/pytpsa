#Copyright (c) 2008, Riccardo De Maria
#All rights reserved.

from pol import *

class polmap(dict):
  """
  >>> from polmap import *
  >>> x,xp,l=mkpol('x,xp,l')
  >>> drift=polmap(x=x+l*xp,xp=xp)
  >>> o=drift(x=0,xp=1.)
  >>> print o
  x=0.0 + l
  xp=1.0
  """

  def __init__(self,*s,**vars):
    dict.__init__(self,*s,**vars)
    for k,v in self.items():
      self[k]=pol(v)

  def __setattr__(self,k,v):
    dict.__setitem__(self,k,v)
    dict.__setattr__(self,k,v)

  __setitem__=__setattr__

  def eval(self,*v,**a):
    """Compose self with other using the python interpreter. The argument can be a dictionary, another polmap or a set of named arguments.
    >>> from polmap import *
    >>> p=polmap(x='1+x+y')
    >>> print p({'x':4})
    x=5.0 + y
    >>> print p(x=4)
    x=5.0 + y
    >>> print p(polmap(x=4))
    x=5.0 + y
    >>> print p(p)
    x=2.0 + x + 2.0*y
    """
    loc={}
    try:
      loc.update(v[0])
    except:
      pass
    loc.update(a)
    for n,v in loc.items():
      if isinstance(v,str):
        loc[n]=pol(v)
    out=polmap()
    for n,v in self.items():
      out[n]=pol(v.eval(**loc))
    return out

  def __call__(self,other):
    return other.eval(self)

  def __repr__(self):
    return getattr(polmap,'out')(self)

  def table(self):
    """Table  print
    """
    out=[]
    for n,v in self.items():
      out.append('%s=\n%s' % (n,v.table()))
    return '\n'.join(out)

  def pretty(self):
    """Pretty print
    """
    out=[]
    for n,v in self.items():
      out.append('%s=%s' % (n,pol(v).pretty()))
    return '\n'.join(out)

  def int(self,x):
    new=polmap(self)
    for n,v in self.items():
      new[n]=v.int(x)
    return new

  def float(self):
    new=polmap(self)
    for n,v in self.items():
      new[n]=float(v)
    return new

  def __mul__(self,other):
    """Compose self with other using a fast but memory consuming algorithm.
    >>> from polmap import *
    >>> p=polmap(x='1+x+x**2')
    >>> print p(p)
    x=3.0 + 3.0*x + 4.0*x**2 + 2.0*x**3 + x**4
    >>> print p*p
    x=3.0 + 3.0*x + 4.0*x**2 + 2.0*x**3 + x**4
    """

    return compose(self,other)

  def __pow__(self,other):
    """ Return the power of a map defined as a sequence of composition
    >>> from polmap import *
    >>> p=polmap(x='x+y',y='2*x-y')
    >>> print p**0
    y=y
    x=x
    >>> print p**1
    y=2.0*x - 1.0*y
    x=x + y
    >>> print p**2
    y=3.0*y
    x=3.0*x
    >>> print p**3
    y=6.0*x - 3.0*y
    x=3.0*x + 3.0*y
    """
    if other==0:
      return polmap( [ (i,i) for i  in self.keys() ])
    else:
      new=polmap(self)
    for i in range(1,other):
      new=new*self
    return new

  def __add__(self,other):
    new=polmap(self)
    for i in other.keys():
      new[i]=self[i]+other[i]
    return new

  def __sub__(self,other):
    new=polmap(self)
    for i in other.keys():
      new[i]=self[i]-other[i]
    return new

  def __neg__(self):
    new=polmap(self)
    for i in other.keys():
      new[i]=-self[i]
    return new

  def matrix(self,vars=None):
    """Return the first order matrix:
    >>> from polmap import *
    >>> f=polmap(x='1+x+l*y',y='y')
    >>> f*=polmap(x='x',y='y+k*x**2')
    >>> print f
    y=y + x**2*k
    x=1.0 + x + y*l + x**2*k*l
    >>> print f.matrix()
    [[1.0, x*k], [l, 1.0 + x*k*l]]
    """
    if not vars:
      vars=self.keys()
    return [ [ self[j].divterm(i) for i in vars] for j in vars if j in self]

  def imatrix(self,vars=None):
    """ Return inverse first order matrix: EXPERIMENTAL
    """
    ((a,b),(c,d))=self.matrix(vars=vars)
    det=a*d-b*c
    return [[d/det,-b/det],[-c/det,a/det]]

  def trace(self):
    vars=self.keys()
    return [self[i].divterm(i)  for i in vars]

  def der(self,var):
    """ Return the derivative with respect to var
    >>> from polmap import *
    >>> f=polmap(x='1+x+l*y',y='y')
    >>> print f.der('y')
    y=1.0
    x=l
    """
    new=polmap(self)
    for j in new.keys():
      new[j]=new[j].der(var)
    return new

  def const(self,vars):
    """ Return the terms that do non depend on vars
    >>> from polmap import *
    >>> f=polmap(x='1+x+l*y',y='y')
    >>> print f.const('ly')
    y=0
    x=1.0 + x
    """
    new=polmap(self)
    for j in self.keys():
      new[j]=new[j].const(vars)
    return new

  __str__=__repr__
  out=pretty


def cdot(a,b):
  """ Dot product
    >>> from polmap import *
    >>> x=polmap(x=1,y=2,z=3)
    >>> print cdot(x,x)
    14.0
  """
  vars=list(set(a).union(set(b)))
  out=0
  for i in vars:
    out+=a.get(i,0)*b.get(i,0)
  return out


def compose(a,b):
  temp={} #bookkeeping vars
  tempvalues={} #bookkeeping values
  newm=polmap() #final map
  lm=0 #couter for temp values
  newv=pol(0) # fix for empty a
  for i in a: #store input pol in temp
    tempvalues[i]=b[i]
  for n,i in a.items(): # for each name,pol in map
    new=pol(0) # initialize new pol
    for j,c in i.coef.items(): # for each mon in pol
      m=[] # prod vector
      for k in range(len(j)): # expand powers
        vv=i.vars[k] # mon var
        if vv in b: # keep
          m+=[i.vars[k]]*j[k] # add has list of mul
        else:
          c*=pol(vv)**j[k] # transfer in the coefficient
#      print "coefficient and monomial"
#      print c,m
#      print "reduce multiplications"
      if m: # if mon is not a constant
        while len(m)>1: # reduce multiplication
          p=m.pop(),m.pop() # terms
          if p in temp: # if terms already calculated
            t=temp[p]  # take it
          else: # calculate new term
            lm+=1  # create temp var
            temp[p]=lm # bookkeep var
            t=lm # store for putting in m
            newv=tempvalues[p[0]]*tempvalues[p[1]] # compute
            tempvalues[t]=newv # store
          m.append(t) # put var
#        print m
        new+=c*tempvalues[m[0]] # add terms, newv carries always the result
      else:
        new+=c # case mon is a constant
#      print new
    newm[n]=new # put in the resulting map
  return newm

def idmap(vars):
  return polmap([(i,i) for i in vars])


