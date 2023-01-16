%colors linux
from sage.all import *
nu,x,y,z,a,b,z,T = var('nu,x,y,z,a,b,z,T')
def kz(zz,a,b,t):
    var ('yy,u')
    n = len(t)
    if n == 1:
        if a == b:
            return 1+a*(t[0]-zz)
        else:
            return b/(b-a)+a*exp((a-b)*(t[0]-zz))/(a-b)
    tmp = 0
    pp = SetPartitions(n).list()
    for p in pp:
        c = 1
        if len(p) >= 2:
            for i in range(len(p)):
                c = c*kz(u,a,b,[t[j-1] for j in p[i]])
            tmp += c
    assume(t[0]-zz>0)
    kk = a * integrate(exp((a-b)*yy)*tmp.subs({u:yy+zz}),(yy,0,t[0]-zz))
    forget(t[0]-zz>0)
    return kk

def c(a,b,t):
    var ('yyy')
    n = len(t)
    tmp = 0
    pp = SetPartitions(n).list()
    for p in pp:
        ee = 1
        for i in range(len(p)):
            ee = ee*kz(yyy,a,b,[t[j-1] for j in p[i]])
        tmp += flatten([ee])[0]
    if hasattr(t[0],'assume'):assume(t[0]>0)
    if hasattr(a,'assume') and hasattr(b,'assume'):
        assume(a<b)
    if isinstance(a+b,type(1)) or isinstance (a+b,type(1.0)):tt = integrate(tmp,(yyy,0,t[0]))
    else:tt = integrate(simplify(tmp),(yyy,0,t[0]))
    if hasattr(a,'assume') and hasattr(b,'assume'):
        forget(a<b)
    if hasattr(t[0],'assume'):forget(t[0]>0)
    return expand(tt)

def m(nu,a,b,t):
    tmp = 0;
    n = len(t)
    if n == 0:
        return 1
    pp = SetPartitions(n).list()
    for p in pp:
        e = 1
        for i in range(len(p)):
            e = e*nu*c(a,b,[t[j-1] for j in p[i]])
        tmp += e
    return flatten([tmp])[0]

t=var('t')
c(1,2,[t,t])

t1,t2 = var('t1,t2')

m(nu,a,b,[t1,t2])

expand(c(1,1,[t1,t2]))
