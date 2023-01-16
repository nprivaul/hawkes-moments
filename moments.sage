%colors linux
from sage.all import *
nu,a,b,t,z,y = var('nu,a,b,t,z,y')

@CachedFunction
def kz(z,n,a,b,t):
  var('yy,u')
  if n == 1:
    if a == b:
      return 1 + a*(t - z)
    else:
      return b/(b - a) + a*exp((a - b)*(t - z))/(a - b)
  tmp = 0
  z1 = []
  for k in range(1,n):
    z1.append(kz(u, k, a, b, t))
    tmp += bell_polynomial(n, n-k+1)(z1)
  assume(z-t>0)
  mmm = a*integrate(exp((a - b)*yy)*tmp.subs({u:yy+z}), (yy, 0, t - z))
  forget(z-t>0)
  return mmm

def c(a, b, t, n):
  var('temp,ccc,z,k,kkk')
  temp = 0
  ccc = []
  for k in range(n,0,-1):
    ccc.append(kz(z, n - k + 1, a, b, t))
    temp += bell_polynomial(n, k)(ccc)
  if hasattr(a,'assume') and hasattr(b,'assume'): assume(a<b)
  if hasattr(t,'assume'): assume(t>0)
  kkk = integrate(temp, (z, 0, t))
  if hasattr(a,'assume') and hasattr(b,'assume'): forget(a<b)
  if hasattr(t,'assume'):forget(t>0)
  return expand(kkk)
  
def m(nu, a, b, t, n):
  tmp = 0
  if n == 0:
    return 1
  ddd = []
  for k in range(n,0,-1):
    ddd.append(nu*c(a, b, t, n - k + 1))
    tmp += bell_polynomial(n, k)(ddd)
  return expand(tmp)

c(1, 2, t, 3)

%time m(nu, a, b, t, 5)