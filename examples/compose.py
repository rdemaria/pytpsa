from pytpsa import *

pol.order=10

p=polmap(x='x+px+y**2',y='y+py+x**2',px='px+y',py='py+x')
n=5

def slow():
  a=p
  for i in range(n-1):
    a=a(p)
#  print a

def fast():
  a=p**n
#  print a
