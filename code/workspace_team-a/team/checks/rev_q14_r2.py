import numpy as np
h=6.626e-27;k=1.381e-16;c=2.998e10
def load(f):
    L=open(f).read().split('\n')
    i=[j for j,l in enumerate(L) if 'NUMBER OF ENERGY' in l][0]
    n=int(L[i+1]);lev=[L[i+3+j].split() for j in range(n)]
    E=np.array([float(x[1]) for x in lev]);g=np.array([float(x[2]) for x in lev])
    j=[j for j,l in enumerate(L) if 'NUMBER OF RADIATIVE' in l][0];m=int(L[j+1])
    tr=[L[j+3+q].split() for q in range(m)]
    return E,g,[(int(t[1])-1,int(t[2])-1,float(t[3]),float(t[4]),float(t[5])) for t in tr]
def W(f,sel,T):
    E,g,tr=load(f);Ek=E*1.4388
    Q=(g*np.exp(-Ek/T)).sum();tot=0
    Tb=2.73
    for u,l,A,nu,Eu in tr:
        if not sel(nu):continue
        nu_=nu*1e9;gu=g[u]
        Nu=gu*np.exp(-Ek[u]/T)/Q
        # thin: int Tb dv [K cm/s]/N = h c^3 A Nu /(8 pi k nu^2)... in K km/s per N
        x=h*nu_/k
        cont=(1/(np.exp(x/T)-1)-1/(np.exp(x/Tb)-1))*x
        tau_over_N=A*c**3/(8*np.pi*nu_**3)*(np.exp(x/T)-1)*Nu  # per cm/s per N
        # Tb = (x)*(1/(e^{x/T}-1)-1/(e^{x/Tb}-1))*(1-e^-tau)
        tot+=cont*tau_over_N*1e-5 # K km/s per N (thin)
    return tot
for T in [6,10,20,30]:
    a=W('/opt/lamda/e-ch3oh.dat',lambda nu:abs(nu-96.739)<0.02 or abs(nu-96.7445)<.01 or abs(nu-96.7555)<.01,T)/2*1
    a+=W('/opt/lamda/ch3oh_a.dat',lambda nu:abs(nu-96.7414)<0.005,T)/2
    b=W('/opt/lamda/e-ch3oh.dat',lambda nu:abs(nu-254.0154)<0.003,T)/2
    print(T,a/b)
