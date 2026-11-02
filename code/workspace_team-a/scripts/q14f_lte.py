# LTE helper for CH3OH (LAMDA e-ch3oh / ch3oh_a): K km/s per cm-2 TOTAL column (N_A=N_E=N/2, each ladder on its own partition sum), CMB contrast on. Importable.
import numpy as np
K1=1.4388; hP=6.62607e-27; cc=2.998e10; kB=1.380649e-16
def _lam(p):
    L=open(p).read().split('\n'); i=[k for k,l in enumerate(L) if 'NUMBER OF ENERGY' in l][0]; n=int(L[i+1])
    lv=np.array([[float(a.split()[1])*K1,float(a.split()[2])] for a in L[i+3:i+3+n]]); j=[k for k,l in enumerate(L) if 'NUMBER OF RADIATIVE' in l][0]
    return lv,np.array([[float(x) for x in l.split()[:6]] for l in L[j+3:j+3+int(L[j+1])]])
_E=_lam('/opt/lamda/e-ch3oh.dat'); _A=_lam('/opt/lamda/ch3oh_a.dat')
def _W(M,sel,T,cmb=True):
    lv,tr=M; Q=np.sum(lv[:,1]*np.exp(-lv[:,0]/T)); s=0.
    for r in tr:
        if not sel(r[4]): continue
        up=int(r[1])-1; Nu=lv[up,1]*np.exp(-lv[up,0]/T)/Q; nu=r[4]; w=hP*cc**3*r[3]*Nu/(8*np.pi*kB*(nu*1e9)**2)/1e5
        if cmb: x=4.799e-11*nu*1e9/T; xb=4.799e-11*nu*1e9/2.73; w*=1-(np.exp(x)-1)/(np.exp(xb)-1)
        s+=w
    return s
S254=lambda nu:abs(nu-254.015)<0.01
def W254(T,cmb=True): return 0.5*_W(_E,S254,T,cmb)
def Qpair(T): return (np.sum(_E[0][:,1]*np.exp(-_E[0][:,0]/T)),np.sum(_A[0][:,1]*np.exp(-_A[0][:,0]/T)))
