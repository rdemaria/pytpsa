from pytpsa import *
d=pol('d1*px+d2*px**2+d3*px**3')
k=pol('a1*x+a2*x**2+a3*x**3')
vars='x px'.split()
pbracket(d,k,vars)


d=pol('d2*px**2+d3*px**3')
k=pol('a2*x**2+a3*x**3')
vars='x px'.split()
pbracket(d,k,vars)


def pb(lst):
  x=lst.pop(0)
  if len(lst)==0:
    return x
  else:
    return pbracket(x,pb(lst),vars)
#    return '[%s,%s]' % (x,pb(lst))

d+k+1./2*pb([d,k])+1./12*pb([d,d,k])-1./12*pb([k,d,k])
