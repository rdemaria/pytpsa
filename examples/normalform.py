nu=.27
s=sin(2*pi*nu)
c=cos(2*pi*nu)
eps=.001
m=polmap(x='c*x+s*p+eps*s*x**2',p='-s*x+c*p+eps*c*x**2')
m=m.eval(c=c,s=s,eps=eps)
linv(m)

m.linear()

ml=polmap(m).truncate()

jacobian(ml)

mlm=array([ [ml.x.coef[(1, 0)],ml.x.coef[(0, 1)]], [ml.p.coef[(1, 0)],ml.p.coef[(0, 1)]]])

mlinv=polmap(
    x='- 0.125333233564*x -0.992114701314*p',
    p='- 0.125333233564*p +0.992114701314*x',)

mlinv(m)

## TOD0
# make sure that parameraters are separated
