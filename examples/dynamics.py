from pytpsa import *

f=polmap(x='x+x**2')
g=polmap(x='x**3')

print (g*f)(x=2)
print g(f(x=2))

