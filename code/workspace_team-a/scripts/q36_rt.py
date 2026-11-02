# Q36: minimal RADEX-equivalent (uniform-sphere escape probability, LAMDA, first collision partner = H2) that allows a
# FREQUENCY-DEPENDENT background mean intensity J_bg(nu) [K, Planck units] instead of a single blackbody Tbg.
# Validated against the team RADEX grids (out/q4c/grid4_*.npz, Tbg=2.73/5/10/20) by q36_validate.py.
import numpy as np
H=6.62607e-27; KB=1.380649e-16; C=2.99792458e10
def parse(fn):
    L=[l.rstrip('\n') for l in open(fn)]
    i=[k for k,l in enumerate(L) if l.upper().startswith('!NUMBER OF ENERGY')][0]; nl=int(L[i+1])
    lev=[L[i+3+k].split() for k in range(nl)]; E=np.array([float(x[1]) for x in lev]); g=np.array([float(x[2]) for x in lev])
    i=[k for k,l in enumerate(L) if l.upper().startswith('!NUMBER OF RADIATIVE')][0]; nt=int(L[i+1])
    rad=[L[i+3+k].split() for k in range(nt)]
    up=np.array([int(r[1])-1 for r in rad]); lo=np.array([int(r[2])-1 for r in rad]); A=np.array([float(r[3]) for r in rad]); nu=np.array([float(r[4]) for r in rad])*1e9
    j=[k for k,l in enumerate(L) if l.upper().startswith('!NUMBER OF COLL TRANS')][0]; nct=int(L[j+1])
    nT=int(L[j+3]); T=np.array(L[j+5].split(),float)
    rows=[L[j+7+k].split() for k in range(nct)]
    cu=np.array([int(r[1])-1 for r in rows]); cl=np.array([int(r[2])-1 for r in rows]); K=np.array([[float(x) for x in r[3:3+nT]] for r in rows])
    return dict(nl=nl,E=E,g=g,up=up,lo=lo,A=A,nu=nu,cT=T,cu=cu,cl=cl,K=K)
def beta(tau):
    tau=np.asarray(tau,float); b=np.ones_like(tau); m=np.abs(tau)>1e-6; t=tau[m]
    if m.any():
        # RADEX uniform sphere (Osterbrock 1989)
        b[m]=1.5/t*(1-2/t**2+(2/t+2/t**2)*np.exp(-np.minimum(t,700.)))
        sm=np.abs(t)<0.01; b[m]=np.where(sm,1-0.375*t,b[m])
    return b
def Jp(nu,T): x=H*nu/(KB*T); return (H*nu/KB)/np.expm1(x)       # Planck J in K
def solve(M,Tk,nH2,Ncol,dv_kms,Jbg,niter=1500):
    """Jbg: callable(nu_array)-> K  or float Tbg.  Returns dict(TR,tau,Tex) per radiative transition (RADEX output convention)."""
    nl=M['nl']; up,lo,A,nu=M['up'],M['lo'],M['A'],M['nu']; g=M['g']
    JB=(lambda f: Jp(f,Jbg)) if np.isscalar(Jbg) else Jbg
    jb=JB(nu)
    # collision rates at Tk (linear interp in T, clamp)
    kd=np.array([np.interp(Tk,M['cT'],r) for r in M['K']]); Cdn=np.zeros((nl,nl)); Cup=np.zeros((nl,nl))
    dE=(M['E'][M['cu']]-M['E'][M['cl']])*H*C
    Cdn[M['cu'],M['cl']]=kd*nH2; Cup[M['cl'],M['cu']]=kd*nH2*g[M['cu']]/g[M['cl']]*np.exp(-dE/(KB*Tk))
    Ncd=Ncol/(dv_kms*1e5*1.0645)   # RADEX: N/(1.0645 dv)
    B21=A*C**2/(2*H*nu**3); B12=B21*g[up]/g[lo]
    # start from LTE; rate eqs use A*beta*(1+n_bg) (RADEX mean-intensity form)
    x=g*np.exp(-M['E']*H*C/(KB*Tk)); x/=x.sum()
    for it in range(niter):
        tau=(C**3/(8*np.pi*nu**3))*A*Ncd*(x[lo]*g[up]/g[lo]-x[up])
        bn=beta(tau)
        nbg=jb/(H*nu/KB)       # photon occupation of background
        W=np.zeros((nl,nl))
        # collisional: W[up,lo]=Cup[lo,up]; W[lo,up]=Cdn[up,lo]
        W+=Cup.T+Cdn.T
        # radiative: spontaneous (escaping fraction beta) + induced by background (beta*nbg)
        for k in range(len(A)):
            u,l=up[k],lo[k]
            W[l,u]+=A[k]*bn[k]*(1+nbg[k]); W[u,l]+=A[k]*bn[k]*nbg[k]*g[u]/g[l]
        Mx=W-np.diag(W.sum(0)); Mx[-1,:]=1.; b=np.zeros(nl); b[-1]=1.
        xn=np.linalg.solve(Mx,b); xn=np.maximum(xn,0)
        if np.max(np.abs(xn-x)/np.maximum(x,1e-30)*(x>1e-8))<1e-6: x=xn; break
        x=0.7*x+0.3*xn
    tau=(C**3/(8*np.pi*nu**3))*A*Ncd*(x[lo]*g[up]/g[lo]-x[up])
    with np.errstate(all='ignore'):
        Tex=(H*nu/KB)/np.log(x[lo]*g[up]/(g[lo]*x[up]))
        TR=(Jp(nu,Tex)-jb)*(1-np.exp(-tau))
    return dict(TR=TR,tau=tau,Tex=Tex,x=x)
