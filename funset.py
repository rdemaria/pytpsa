#Copyright (c) 2008, Riccardo De Maria
#All rights reserved.

import cmath as math
from pol import pol, phorner

pi=math.pi

""" The idea is that you evaluate a+x being x a small quantity,
    then you rewrite f(a+x) so that you have a combination of
    g_i(a) and g_k(x h_l(a) ).
    g_i(a) and h_l(a) can be evaluated exactly by usual math libraries
    g_k(x) because are proportianal to x are well approximated
    by a taylor expansion.
"""


def sqrt(c):
  """Square root of a polynomial
    >>> from pol import *
    >>> print sqrt(pol('1+x'))**2
    1.0 + x
    >>> print sqrt(pol('1j+x'))**2
    1j -1j*x
  """
  if not isinstance(c,pol): return math.sqrt(c)
  a0,p=c.separate(); p/=a0
  lst=[math.sqrt(a0)]
  for n in range(1,c.order+1):
    lst.append(-lst[-1]/2/n*(2*n-3))
  return phorner(lst,p)

def exp(c):
  """Exponential of a polynomial
    exp(a+x) = exp(a)exp(x)
    >>> from pol import *
    >>> print log(exp(pol('1+x')))
    1.0 + x
  """
  if not isinstance(c,pol): return math.exp(c)
  a0,p=c.separate();
  lst=[exp(a0)]
  for n in range(1,c.order+1):
    lst.append(lst[-1]/n)
  return phorner(lst,p)

def log(c):
  """Logarithm of a polynomial
    log(a+x)=log(a)-log(1+x/a)
    >>> from pol import *
    >>> print log(exp(pol('1+x')))
    1.0 + x
  """
  if not isinstance(c,pol): return math.log(c)
  a0,p=c.separate(); p/=a0
  lst=[log(a0),1.]
  for n in range(2,c.order+1):
    lst.append( -lst[-1]/n*(n-1)   )
  return phorner(lst,p)

def sin(c):
  """
  sin(a+x)= sin(a) cos(x) + cos(a) sin(x)
  """
  if not isinstance(c,pol): return math.sin(c)
  a0,p=c.separate();
  lst=[math.sin(a0),math.cos(a0)]
  for n in range(2,c.order+1):
    lst.append( -lst[-2]/n/(n-1))
  return phorner(lst,p)

def cos(c):
  """
  cos(a+x)= cos(a) sin(x) - sin(a) cos(x)
  """
  if not isinstance(c,pol): return math.cos(c)
  a0,p=c.separate();
  lst=[math.cos(a0),-math.sin(a0)]
  for n in range(2,c.order+1):
    lst.append( -lst[-2]/n/(n-1))
  return phorner(lst,p)


def tan(y):
  """Compute Tan
  """
  return sin(y)/cos(y)



def sinh(c):
  """Compute Sinh using a Taylor expansion
  """
  if not isinstance(c,pol): return math.sinh(c)
  a0,p=c.separate();
  lst=[math.sinh(a0),math.cosh(a0)]
  for n in range(2,c.order+1):
    lst.append( lst[-2]/n/(n-1))
  return phorner(lst,p)

def cosh(c):
  """Compute Cosh using a Taylor expansion
  """
  if not isinstance(c,pol): return math.cosh(c)
  a0,p=c.separate();
  lst=[math.cosh(a0),math.sinh(a0)]
  for n in range(2,c.order+1):
    lst.append( lst[-2]/n/(n-1))
  return phorner(lst,p)

def isqrt(c):
  if not isinstance(c,pol): return 1/math.sqrt(c)
  a0,p=c.separate(); p/=a0
  lst=[1/math.sqrt(a0)]
  for n in range(1,c.order+1):
    lst.append(-lst[-1]/a0/2/n*(2*n-1))
  return phorner(lst,p)


#def asin_t(c):
#  """Compute ArcSin using a Taylor expansion
#    >>> from pol import *
#    >>> x,y=mkpol('x,y')
#    >>> asin(sin(.7+x+y))
#    0.7 + x + y
#  """
#  p=c.copy()
#  x=pol('x')
#  return pisqrt(1-x**2).int(x)(x=p)

def asin(y):
  """Compute ArcSin
    >>> from pol import *
    >>> x,y=mkpol('x,y')
    >>> asin(sin(.4+x+y))
    0.4 + x + y
  """
  x0=pol(math.asin(y.zero()))
  for i in range(y.order):
    x0=x0 + (y-sin(x0))/cos(x0)
  return x0


def acos(y):
  """Compute ArcCos
    >>> from pol import *
    >>> x,y=mkpol('x,y')
    >>> acos(cos(.4+x+y))
    0.4 + x + y
  """
  x0=pol(math.acos(y.zero()))
  for i in range(y.order):
    x0=x0 + -(y-cos(x0))/sin(x0)
  return x0


def atan(y):
  """Compute ArcTan using Newton method
    x=f(y); y=g(x)
    x0=f(y0); x0=x0 +(y-g(x0))/g'(x)
    >>> from pol import *
    >>> p=pol('.4+x+y')
    >>> print tan(atan(p))
    0.4 + x + y
    >>> print atan(tan(p))
    0.4 + x + y
  """
  x0=pol(math.atan(y.zero()))
  for i in range(y.order):
    x0=x0 + (y-tan(x0))*cos(x0)**2
  return x0



def asinh(y):
  """Compute ArcSinh
    >>> from pol import *
    >>> x,y=mkpol('x,y')
    >>> asinh(sinh(.4+x+y))
    0.4 + x + y
  """
  x0=pol(log(y+sqrt(y**2+1)))
  for i in range(y.order):
    x0=x0 + (y-sinh(x0))/cosh(x0)
  return x0


def acosh(y):
  """Compute ArcCosh
    >>> from pol import *
    >>> x,y=mkpol('x,y')
    >>> acosh(cosh(1.4+x+y))
    1.4 + x + y
  """
  x0=pol(log(y+sqrt(y**2-1)))
  for i in range(y.order):
    x0=x0 + (y-cosh(x0))/sinh(x0)
  return x0

def newton(f,y,x0):
  """Inverse polynomial using newton method EXPERIMENTAL
  >>> f=pol('1+x**2')
  >>> f(x=newton(f,pol('2+x'),1))
  2.0 + x
  """
  for i in range(y.order):
    x0=x0+ ( y-f(x=x0) ) / f.der('x')(x=x0)
  return x0


#if __name__=='__main__':
#  import doctest
#  doctest.testmod()
#  import profile
#  pol.order=9
#  profile.run('pol("sqrt(1+x+y+z+px+py+pz)")',sort='time')


