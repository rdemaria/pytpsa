#Copyright (c) 2008, Riccardo De Maria
#All rights reserved.


from pol import pol, normder
from polmap import polmap


class dop(pol):
  def __init__(self,val=None,order=None,eps=1E-13,loc={},m='eval'):
    self.coef={}
    self.vars=[]
    self.order=1000
    self.eps=eps
    if val!=None:
      if isinstance(val,dop):
        self.coef.update(val.coef)
        self.vars=val.vars[:]
        self.order=val.order
        self.eps=val.eps
      elif isinstance(val,str):
        if m=='eval':
          c=compile(val,'eval','eval')
          l=dict( (i,dop(i,m='name')) for i in c.co_names)
          l.update(loc)
          dop.__init__(self,eval(c,{},l))
        elif m=='name':
          self.vars=[val]
          self.coef[(1,)]=pol(1.)
    if order is not None:
      self.order=order

  def __call__(self,other):
    if type(other) is pol:
      return other.der(self)
    elif type(other) is polmap:
      new=polmap(other)
      for j in new.keys():
        new[j]=new[j].der(self)
      return new

  def __call__(self,other):
    new=pol(order=min(self.order,other.order))
    new.vars=list(set(other.vars+self.vars))
    for i in self.coef:
      for j in other.coef:
        c=self.coef[i]*other.coef[j]
        expi=dict(zip(self.vars,i))
        expj=dict(zip(other.vars,j))
        newexp=tuple([-expi.get(k,0)+expj.get(k,0) for k in new.vars])
        for k in self.vars:
          c*=normder(expj.get(k,0), expi.get(k,0))
        new.coef[newexp]=new.coef.get(newexp,0.)+c
    return new.truncate().dropneg()

#  def truncate(self):
#    return self


def pop(h,vars):
  h=pol(h)
  vars=list(vars)
  out=0
  while vars:
    x=vars.pop(0)
    p=vars.pop(0)
    out+=dop(x)*(-dop(p)(h))+dop(p)*dop(x)(h)
  return out

  
