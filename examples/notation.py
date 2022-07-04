from pytpsa import *

# definition
x=pol('x')
p=pol('p')
nu=.27
s=sin(2*pi*nu)
c=cos(2*pi*nu)
eps=.001



#map definition
h1=polmap(x=c*x+s*p,p=-s*x+c*p)
h2=polmap(x=x,p=p+eps*x**2)

#observable definition
o=pol('x**2')

print o(h1)
print h1 * o
print h1(h2)
print h1*h2

m=h1*h2
r=tomatrix(m,'xp')
print r
print frommatrix(r,'xp')

mi=linv(m)

print mi*m
print m*mi

print h1*h2*x
print x(h1*h2)
print x(h1(h2))

print h2*h1*x
print x(h2*h1)
print x(h2(h1))




# matrix multiplication
a=polmap(x='a11*x+a12*p',p='a21*x+a22*p')
b=polmap(x='b11*x+b12*p',p='b21*x+b22*p')

print a*b


x=pol('x')
p=pol('p')
nu=.27
s=sin(2*pi*nu)
c=cos(2*pi*nu)
eps=.001
#map definition
h1=polmap(x=c*x+s*p,p=-s*x+c*p)
h2=polmap(x=x,p=p+eps*x**2)

m=h1(h2)
r=tomatrix(m,'xp')
r=lin(m)
ri=linv(m)

rb=polmap(x='.5*(hp+hm)',p='.5j*(hp-hm)')



