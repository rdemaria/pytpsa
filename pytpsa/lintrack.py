from pytpsa import *
from numpy import sign

quad=polmap(x='x',px='px+kn1l*x',y='y',py='py-kn1l*y')
bend=polmap(x='x+l*px+l*angle*dpp/2.',px='px+angle*dpp',y='y+l*py',py='py')
drift=bend(a=0)

def tmatrix(m):
  return m.matrix(vars='x px y py dpp'.split())

def pbeta(r11,r12,r21,r22):
  cmu=(r11+r22)/2.
  try:
    smu=sign(r12.zero()) * sqrt(-r12*r21 - (r11-r22)**2 / 4.)
    mu=.5/pi * atan(smu/cmu)
    ismu=1/smu
  except ZeroDivisionError:
    mu=pol('mu')
    smu=pol('smu')
    ismu=pol('ismu')
  bet=r12*ismu
  alf=.5*(r11-r22)*ismu
  return bet,alf,cmu,mu

def pdisp(r11,r12,r21,r22,r16,r26):
  det==r11*r22-r12*r21
  d=(-r16*r22+r26*r12) / det
  dp=(-r26*r11+r16*r21) / det
  return d,dp

def twiss(m):
  r=tmatrix(m)
  betx,alfx,cmux,mux=pbeta(r[0][0],r[0][1],r[1][0],r[1][1])
  bety,alfy,cmuy,muy=pbeta(r[2][2],r[2][3],r[3][2],r[3][3])
  return betx,alfx,cmux,mux,bety,alfy,cmuy,muy

def tfquad(l,kn0l,ks0l,kn1l):
  kx=(+kn1l/l+(kn0l/l)**2)
  skx=sqrt(kx)
  cx=cos(skx*l)
  sx=sin(skx*l)/skx
  dx =(1-cx)/kx
  fx =(l-sx)/kx
  ky=(-kn1l/l+(ks0l/l)**2)
  sky=sqrt(ky)
  cy=cos(sky*l)
  sy=sin(sky*l)/sky
  dy =(1-cy)/ky
  m=polmap()
  m.x =pol('x+cx*x+sx*px+dx*kn0l/l*dpp')
  m.px=pol('px-kx*sx*x+cx*px+sx*kn0l/l*dpp')
  m.y =pol('y+cy*y+sy*py')
  m.px=pol('py-ky*sy*y+cy*px')


