def mkexp1(n,k):
  out=[]
  for i in range(n+1):
    nl=(i,)
    if k==1:
      out.append(nl)
    else:
      for j in  mkexp1(n-i,k-1):
        out.append(nl+j)
  return out

n,k=6,3
for iexp,exp in enumerate(mkexp1(n,k)):
   print exp,iexp


print (n**2+3*n+2)/2


def mkexpn(n,k):
  if k==1:
    return n+1
  else:
    out=0
    for i in range(n+1):
        out+=mkexpn(n-i,k-1)
    return out

mkexpn(n,k)

def mkexpq(jl,n,o):
  j=jl[0]
  if n==1:
    return j
  else:
    out=0
    for i in range(j):
      out+=mkexpn(n-1,o-i)
    out+=mkexpq(jl[1:],n-1,o-j)
    return out

i2e=mkexp1(3,5)
i2e[mkexpq([2,1,1],3,5)]



def mkexp2(*alg):
  alg=list(alg)
  n,o=alg.pop()
  out=[]
  if len(alg)==0:
    for i in mkexp1(n,o):
      out.append(i)
  else:
    for i in mkexp1(n,o):
      for j in mkexp2(*alg):
        out.append(i+j)
  return out

mkexp2((3,5),(5,2))

def mkalg(*alg):
  i2e=mkexp2(*alg)
  e2i={}
  for i,e in enumerate(i2e):
    e2i[e]=i
  i2e=array(i2e)
  return i2e,e2i

i2e,e2i=mkalg((6,2),(5,1))

x=random.rand(len(i2e))
y=random.rand(len(i2e))

def mul(x,y):
  r=zeros(len(x),dtype=x.dtype)
  for ai in range(len(x)):
    for bi in range(len(y)):
      ce=tuple( i2e[ai]+i2e[bi] )
      #ce=i2e[ai]+i2e[bi]
      #ce=tuple(map(sum,i2e[ai],i2e[bi]))
      ci=e2i.get(ce)
      if ci is not None:
        r[ci]+=x[ai]*y[bi]
  return r

%timeit mul(x,y)

def ncomp(n, k):
   """Number of composition of n in k
       (n + k - 1)! / n! (k-1)!
       (n+1)  (n + k -1) / 1 ... k
   """
   result = 1
   for i in range(1, k):
       result = result * (n+i) / i
   return result

def nallcomp(n,k):
   """Sum of the number of composition of i in k for i<=n"""
   return ncomp(n,k+1)

def comp(n,k):
  if k>1:
    for j in range(n+1):
       for rest in comp(j,k-1):
         yield [n-j]+rest
  else:
    yield [n]

def allcomp(n,k):
  for nn in range(n+1):
    for exp in comp(nn,k):
       yield exp


def myhash(exp,mm):
  return sum([ (n+1)**i*exp[i] for i in range(k)])%mm

def myhash(exp,mm):
  res=5432
  for i in exp:
     res=res*33+(i+53)
  return res%



n,k=5,3
for iexp,exp in enumerate(allcomp(n,k)):
   print exp,sum(exp),iexp






def countrep(l):
  count={}
  for h in l:
    count[h]=count.get(h,0)+1
  return max(count.values())

for n in range(1,10):
  for k in range(1,10):
    nmin=nallcomp(n,k)
    hashes=[ myhash(exp,nmin) for exp in allcomp(n,k)]
    print n,k,countrep(hashes)



len(list(allcomp(n,k)))
nallcomp(n,k)

# new trial 2012/08/27

def ncomp(n, k):
   """Number of composition of n in k
       (n + k - 1)! / n! (k-1)!
       (n+1)  (n + k -1) / 1 ... k
   """
   result = 1
   for i in range(1, k):
       result = result * (n+i) / i
   return result

def nallcomp(n,k):
   """Sum of the number of composition of i in k for i<=n"""
   return ncomp(n,k+1)

def comp(n,k):
  if k>1:
    for j in range(n+1):
       for rest in comp(j,k-1):
         yield [n-j]+rest
  else:
    yield [n]

def allcomp(n,k):
  for nn in range(n+1):
    for exp in comp(nn,k):
       yield exp

n,k,ii=5,3,35
exp2il=list(allcomp(n,k))
#for iexp,exp in enumerate(exp2il):
#     print exp,sum(exp),iexp

ee=exp2il[ii]
print ii,exp2i(ee,n,k)


for ii in range(len(exp2il)):
  ee=exp2il[ii]
  print ii,exp2i(ee,n,k)


n,k=5,3
exp2i=dict([tuple(ex),i] for i,ex in  enumerate(allcomp(n,k)))

exp2it=[]
for iex,ex in  enumerate(allcomp(n,k)):
  curr=exp2it
  print iex,ex
  for iee,ee in enumerate(ex):
    print iee,ee,
    if len(curr)<iee+1:
      if iee==k-1:
        ncurr=iex
      else:
        ncurr=[]
      curr.append(ncurr)
      curr=ncurr
    else:
      curr=curr[iee]
    print exp2it


  print ex,i





