class pol(list):
  def __init__(self,s):
    if isinstance(s,str):
      self.extend([0,1,s])
    elif isinstance(s,pol):
      self.extend(s)

  def copy(self):
    n=pol(self)
    if isinstance(n[0],pol):
      n[0]=n[0].copy()
    if isinstance(n[1],pol):
      n[1]=n[1].copy
    return n

  def __add__(a,b):
    a=a.copy()
    if isinstance(b,float):
      a[0]+=b
    return a

  def __mul__(a,b):
    a=a.copy()
    if isinstance(b,float):
      a[0]*=b
      a[1]*=b
    return a


