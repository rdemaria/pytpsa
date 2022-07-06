#Copyright (c) 2008, Riccardo De Maria
#All rights reserved.

import math
from functools import reduce

pi=math.pi

_mabs=abs

def abs(a):
  if type(a) is pol:
    a0,c=a.separate()
    if c:
      return 1e20
    else:
      return _mabs(a0)
  return _mabs(a)

def normint(n,m):
  """coefficient of int^m x^n dx"""
  c=1
  for i in range(m):
    c*=n+i+1
  return c

def normder(n,m):
  """coefficient of d^m / dx^m x^n"""
  c=1
  for i in range(m):
    c*=n-i
  return c

def mkpol(r):
  return map(pol,r.split(','))

def pinv(c):
  a0,p=c.separate(); p/=a0
  lst=[1/a0]
  for n in range(1,c.order+1):
    lst.append(-lst[-1])
  return phorner(lst,p)

def phorner(lst,p):
  out=lst.pop()
  for i in reversed(lst):
    out=out*p+i
  return out


class mydict(dict):
  """ Dictionary for evaluation"""
  def __getitem__(self,k):
    if not k in self:
      return pol(k)
    else:
      return dict.__getitem__(self,k)

class pol(dict):
  """ Class for multivariate polinomials
    >>> from pol import *
    >>> print pol('x**3+z-5-x**2')
    - 5.0 + z - 1.0*x**2 + x**3
    >>> x=pol('x')
    >>> y=pol('y')
    >>> c=pol(pi)
    >>> p=c+pi*x+y**2
    >>> print p
    3.14159265359 + 3.14159265359*x + y**2
    >>> p2=pol(p)
    >>> print p(x=p2,y=p)
    22.8808014558 + 29.6088132033*x + 9.86960440109*x**2 + 9.42477796077*y**2 + 6.28318530718*y**2*x + y**4
    >>> pol.out='table'
    >>> print p(x=p2,y=p)
                     o   y  x
    22.8808014558    0   0  0
    29.6088132033    1   0  1
     9.86960440109   2   0  2
     9.42477796077   2   2  0
     6.28318530718   3   2  1
     1.0             4   4  0
    >>> pol.out='pretty'
    >>> print pol('1j*x')
    1j*x
    """
  out='pretty'
  def __init__(self,val=None,order=None,eps=1E-18,loc={},m='eval'):
    self.vars=[]
    self.order=10
    self.eps=eps
    if val!=None:
      if isinstance(val,self.__class__):
        self.update(val)
        self.vars=val.vars[:]
        self.order=val.order
        self.eps=val.eps
      elif isinstance(val,str):
        if m=='eval':
          c=compile(val,'eval','eval')
          l=dict( (i,pol(i,m='name')) for i in c.co_names if i not in globals())
          l.update(loc)
          pol.__init__(self,eval(c,globals(),l))
        elif m=='name':
          self.vars=[val]
          self[(1,)]=1.
      else:
        self.vars=[]
        self[()]=val
    if order is not None:
      self.order=order

  def zero(self):
    """Extract zero order
    """
    return self.get((0,)*len(self.vars),0.)

  def linear(self):
    a=[1]+[0]*(len(self.vars)-1)
    out=[]
    for i in self.vars:
      out.append( self.get(tuple(a),0) )
      a.insert(0,a.pop())
    return out

  def getlind(self,v):
    ind=[0]*len(self.vars)
    ind[ self.vars.index(v)]=1
    return tuple(ind)

  def getlcoef(self,v):
    return self.get(self.getlind(v),0.)

  def setlcoef(self,v,val):
    self[self.getlind(v)]=val

  def separate(self):
    """Return a copy in couple of zero and high order
    """
    c=self.__class__(self)
    a0=c.pop((0,)*len(c.vars),0.)
    return a0,c

  def truncate(self,order=None,eps=None):
    """Truncate to the order indicate in pol.order and 
    damp the elements smaller than self.eps"""
    if not order: order=self.order
    if not eps: eps=self.eps
    for k in list(self.keys()):
      if sum(k)>order or abs(self[k])<eps:
        del self[k]
    return self

  def const(self,vars=None):
    """Extract terms that do not depends on vars
    >>> from pol import *
    >>> print pol('x**3+y+l').const()
    0.0
    >>> print pol('x**3+y+l').const('yl')
    x**3
    """
    if not vars:
      return self.get((0,)*len(self.vars),0.)
    else:
      c=self.__class__(self)
      for exp in self:
        expd=dict(zip(self.vars,exp))
        if sum(expd.get(j,0) for j in vars) >0 :
          del c[exp]
      return c

  def dropneg(self):
    """ Delete terms with negative exponent"""
    for k in self.keys():
      for j in k:
        if j<0:
          del self[k]
          break
    return self

  def addcoef(self,other):
    """Add a number to pol"""
    new=self.__class__(self)
    i=(0,)*len(new.vars)
    new[i]=new.get(i,0.)+other
    return new.truncate()

  def mulcoef(self,other):
    """Mul coef to pol"""
    new=self.__class__(self)
    for i in new:
      new[i]*=other
    return new.truncate()

  def reorder(self,vars):
    new=self.__class__(order=self.order)
    new.vars=vars
    for exp in self:
      expd=dict(zip(self.vars,exp))
      newexp=tuple([expd.get(j,0) for j in new.vars])
      new[newexp]=self[exp]
    return new

  def addpol(self,other):
    """Add pol to pol """
    new=self.__class__(order=min(self.order,other.order))
    new.vars=list(set(other.vars+self.vars))
    for exp in self:
      expd=dict(zip(self.vars,exp))
      newexp=tuple([expd.get(j,0) for j in new.vars])
      new[newexp]=new.get(newexp,0.)+self[exp]
    for exp in other:
      expd=dict(zip(other.vars,exp))
      newexp=tuple([expd.get(j,0) for j in new.vars])
      new[newexp]=new.get(newexp,0.)+other[exp]
    return new.truncate()

  def mulpol(self,other):
    """Mul pol to pol """
    if other.vars==self.vars:
      return self.fmulpol(other)
    else:
      new=pol(order=min(self.order,other.order))
      new.vars=list(set(other.vars+self.vars))
      for i in self:
        for j in other:
          c=self[i]*other[j]
          expi=dict(zip(self.vars,i))
          expj=dict(zip(other.vars,j))
          newexp=tuple([expi.get(k,0)+expj.get(k,0) for k in new.vars])
          new[newexp]=new.get(newexp,0.)+c
      return new.truncate()

  def fmulpol(self,other):
    """fast mul pol to pol """
    new=self.__class__(order=min(self.order,other.order))
    new.vars=self.vars[:]
    for i in self:
      for j in other:
        c=self[i]*other[j]
        newexp=tuple([l+m for l,m in zip(i,j)])
        new[newexp]=new.get(newexp,0.)+c
    return new.truncate()

  def divterm(self,other):
    """Extract a term from a pol EXPERIMENTAL:
    >>> from pol import *
    >>> r=pol('x**3+z-5').divterm('x**3')
    >>> print r
    1.0
    """
    other=pol(other)
    new=pol(order=min(self.order,other.order))
    other=pol(other)
    new.vars=list(set(other.vars+self.vars))
    for i in self:
      for j in other:
        c=self[i]*other[j]
        expi=dict(zip(self.vars,i))
        expj=dict(zip(other.vars,j))
        newexp=tuple([expi.get(k,0)-expj.get(k,0) for k in new.vars])
        new[newexp]=new.get(newexp,0.)+c
    return new.truncate().dropneg()

  def __pow__(self,n):
    new=self.__class__(self)
    if n==0:
      return 1
    elif n <0:
      return pinv(self)**n
    else:
      for i in range(n-1):
        new*=self
    return new

  def __add__(self,other):
    """Addition
    >>> from pol import *
    >>> x=pol('x')
    >>> 1-x+x
    1.0
    """
    if isinstance(other,self.__class__):
      return self.addpol(other)
    else:
      return self.addcoef(other)

  def __radd__(self,other):
    return self.addcoef(other)

  def __mul__(self,other):
    """Addition
    >>> from pol import *
    >>> x=pol('x')
    >>> (1+1.0*x)/-(-x-1)
    1.0
    """
    if isinstance(other,self.__class__):
      return self.mulpol(other)
    else:
      return self.mulcoef(other)

  def __rmul__(self,other):
    return self.mulcoef(other)

  def __sub__(self,other):
    if isinstance(other,pol):
      return self.addpol(-other)
    else:
      return self.addcoef(-other)

  def __rsub__(self,other):
    return (-self).addcoef(other)

  def __truediv__(self,other):
    if isinstance(other,pol):
      return self.mulpol(pinv(other))
    else:
      return self.mulcoef(1/other)

  def __rtruediv__(self,other):
    return pinv(self).mulcoef(other)

  def __neg__(self):
    return self.mulcoef(-1)

  def __pos__(self):
    return self

  def eval(self,*args,**loc):
    if len(args)>0:
      loc.update(args[0])
    for i in self.vars:
      loc.setdefault(i,pol(i))
    loc=mydict(loc)
    return eval(self.pretty(),{},loc)

  __call__=eval

  def _pexp(self,i):
    out=[]
    for j in range(len(i)):
      if i[j]==1:
        out.append( '%s' % self.vars[j] )
      elif i[j]!=0.:
        out.append( '%s**%d' %(self.vars[j],i[j]) )
    return '*'.join(out)


  def _pcoeff(self,c,i):
    if isinstance(c,complex):
      if abs(c.imag)<self.eps:
        return self._pcoeff(c.real,i)
      else:
        c='+ '+str(c)
    elif isinstance(c,float):
      sign=c<0 and '-' or '+'
      if abs(c-1.0)<self.eps and i:
        c='%s %s' % (sign,i)
        i=''
      else:
        c='%s %s' % (sign,abs(c))
    else:
      c='+ (%s)' % c
    if i:
      return '%s*%s' % (c,i)
    else:
      return c


  def pretty(self):
    lst=sorted([ (sum(i),i,c) for i,c in self.items()])
    m=[]
    for o,i,c in lst:
      i=self._pexp(i)
      c=self._pcoeff(c,i)
      m.append( c )
    if m:
      m=' '.join(m)
      if m.startswith('+ '):
        return '  '+m[2:]
      else:
        return m
    else:
      return '0'

  def table(self):
    fvar=lambda x,y:'%3s%3s' %(x,y)
    out=[['',0,0,reduce(fvar,self.vars),' o']]
    lst=sorted([ (sum(i),i) for i in self])
    rmax,cmax=0,0
    for order,exp in lst:
      coef=str(self[exp])
      r=[coef,len(coef),coef.find('.'),reduce(fvar,exp),str(order)]
      rmax=r[1]>rmax and r[1] or rmax
      cmax=r[2]>cmax and r[2] or cmax
      out.append(r)
    nout=[]
    for c,l,p,e,o in out:
      c=('%%-%ds'%(rmax+cmax)) % (' '*(cmax-p)+ c)
      nout.append('%(c)s %(o)2s %(e)s' % locals() )
    return '\n'.join(nout)

  def __repr__(self):
    return getattr(pol,pol.out)(self)

#if __name__=='__main__':
#  import doctest
#  doctest.testmod()
#  import profile
#  pol.order=9
#  profile.run('pol("sqrt(1+x+y+z+px+py+pz)")',sort='time')

