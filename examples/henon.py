from pytpsa import *
setorder(4)

drift=polmap(x='x+l*px',px='px')
quad1=polmap(x='x',px='px+k1*x')
quad2=polmap(x='x',px='px+k2*x')
sext=polmap(x='x',px='px+k3*x**2')

oneturn=drift*quad1*drift*quad2
oneturn.trace()

drift=polmap(x='x+l*px',px='px',y='y+l*py',py='py')
quad1=polmap(x='x',px='px+k1*x',y='y',py='py-k1*y')
quad2=polmap(x='x',px='px+k2*x',y='y',py='py-k2*y')
sext=polmap( x='x',px='px+k3*(x**2-y**2)',y='y',py='py-k3*y*x')

oneturn=drift*quad1*drift*quad2
oneturn.trace()
[ sum([oneturn[i].divterm(i) for i in 'x px'.split()]) for j in ['x px','y py'] ]

oneturn=oneturn(l=1,k1=0.1,k2=-0.12,k3=0.001)


