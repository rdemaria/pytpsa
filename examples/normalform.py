from pytpsa import *

# build map
nu=.27
s=sin(2*pi*nu)
c=cos(2*pi*nu)
eps=.001
m=polmap(x='c*x+s*p+eps*s*x**2',p='-s*x+c*p+eps*c*x**2')
vars='xp'
M=m.eval(c=c,s=s,eps=eps).reorder(vars)

# linearpart
R=linearpart(M)
keep=[]
partM=[M.copy()]

#diagonalization
up=pol('x+1j*p')
um=pol('x-1j*p')
lp=up(R).getlcoef('x')
lm=um(R).getlcoef('x')
l=[lp,lm]

D=polmap(x='.5*(x+p)',p='-.5j*(x-p)').reorder(vars)
Dinv=polmap(x='x+1j*p',p='x-1j*p').reorder(vars)
D(Dinv)
E=idmap('xp')
polmap.out=polmap.table

# restart here

#separation
G=M-R

#go to resonance basis

R_R = Dinv * R * D # Dinv(R(D))
G_R = Dinv * G * D
M_R = R_R + G_R

#solve homological equation
T_R=G_R.copy()
for i in vars:
  for j in G_R[i]:
    nf=normalfactor(l,vars.find(i),j)
    print i,j,G_R[i][j],abs(nf)
    if abs(nf) > 1e-10 :
      T_R[i][j]=G_R[i][j]/nf
    else:
      T_R[i][j]=0
      keep.append(j)

R_R * T_R - T_R * R_R + G_R # show that I get rid of G_R

# solve in resonance basis
A_R    = E + T_R + T_R*(1+M+M**2)
A_R_inv= E - T_R
M2_R= A_R_inv * M_R * A_R

# solve in physical basis
T= D * T_R * Dinv
R * T - T * R + G
A=E+T
Ainv=E-T
M2=Ainv * M * A
partM.append(M2)

# iterate
M=M2



Aj=jacobian(A)
Aj[0,0]*Aj[1,1]-Aj[1,0]*Aj[1,0]



from pytpsa import *
m=polmap(x='c*x+s*p+eps*s*x**2',p='-s*x+c*p+eps*c*x**2')
rinv=polmap(x='c*x-s*p',p='s*x+c*p')


