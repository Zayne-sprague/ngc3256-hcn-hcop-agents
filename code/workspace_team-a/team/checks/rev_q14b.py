import numpy as np
def load(f):
    L=open(f).read().split('\n')
    n=int(L[5]); lv={}
    for l in L[7:7+n]:
        p=l.split(); lv[int(p[0])]=(float(p[1]),float(p[2]),p[3])
    j=[k for k,l in enumerate(L) if 'NUMBER OF RADIATIVE' in l][0]
    m=int(L[j+1]); tr=[]
    for l in L[j+3:j+3+m]:
        p=l.split(); tr.append((int(p[1]),int(p[2]),float(p[3]),float(p[4]),float(p[5])))
    return lv,tr
def W(f,sel,T):
    lv,tr=load(f)
    Q=sum(g*np.exp(-1.4388*e/T) for e,g,_ in lv.values())
    s=0
    for u,l,A,nu,Eu in tr:
        if sel(nu): s+=A*lv[u][1]*np.exp(-Eu/T)/nu**2
    return s/Q
for T in [8,10,15,20,30,50]:
    e=W('/opt/lamda/e-ch3oh.dat',lambda nu:abs(nu-254.015)<0.01,T)
    h=W('/opt/lamda/e-ch3oh.dat',lambda nu:96.73<nu<96.76,T)+W('/opt/lamda/ch3oh_a.dat',lambda nu:96.73<nu<96.76,T)
    print(T,h/e)
