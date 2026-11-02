import numpy as np
def load(f):
    L=open(f).read().split('\n')
    i=[k for k,l in enumerate(L) if 'NUMBER OF ENERGY LEVELS' in l][0]
    n=int(L[i+1]); lv={}
    for l in L[i+3:i+3+n]:
        p=l.split(); lv[int(p[0])]=(float(p[1]),float(p[2]))
    j=[k for k,l in enumerate(L) if 'NUMBER OF RADIATIVE' in l][0]
    m=int(L[j+1]); tr=[]
    for l in L[j+3:j+3+m]:
        p=l.split(); tr.append((int(p[1]),int(p[2]),float(p[3]),float(p[4]),float(p[5])))
    return lv,tr
c=2.998e10;h=6.626e-27;k=1.381e-16
for f in ['/opt/lamda/e-ch3oh.dat','/opt/lamda/ch3oh_a.dat']:
    lv,tr=load(f)
    sel=[t for t in tr if (abs(t[4]-96.739)<0.03) or abs(t[4]-254.015)<0.01]
    print(f)
    for t in sel: print(t, lv[t[0]], lv[t[1]])
