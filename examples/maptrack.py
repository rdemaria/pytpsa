from pytpsa import *

class Drift(object):
  def __init__(self,l):
    self.l=l
  def track(self,m):
    m['x']+=self.l*m['px']
    return m

class Kick(object):
  def __init__(self,k,order=1):
    self.k=k
    self.order=order
  def track(self,m):
    m['px']+=self.k*m['x']**self.order
    return m

line=[Drift(5),Kick(0.03,1),Kick(.0003,5),Drift(5),Kick(-.03)]
p=polmap(x='0.01+x',px='px')
p['x'].order=1
p['px'].order=1

for elem in line:
  p=elem.track(p)

print p

