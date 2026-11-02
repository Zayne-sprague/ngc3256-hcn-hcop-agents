# Q3 round 2 library: fine (interpolated) RADEX grid, fast Gaussian-likelihood engine with analytic ln(eta) marginalisation (eta<=cap), recovery nuisance model.
# Grid: Tk 31 pts (0.05 dex, 10-300 K), log n 41 pts (0.125 dex, 3-8), log N(HCO+) 0.1 dex steps (11.5-19), log A = log N(HCN)/N(HCO+) 0.1 dex steps (-1..1.5); ln T_R interpolated
# linearly in (log Tk, log n, log N) from the 11x21x31 RADEX LVG grid (out/q3/grid_lvg.npz, s52). Cells where any neighbour RADEX model failed (T_R<=0 or tau<=0) are excluded.
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
from scipy.interpolate import RegularGridInterpolator
from scipy.special import log_ndtr
import s53_ladder_fit as L
OUT='/workspace/out/q3b/'
LT=np.linspace(1,np.log10(300),31); LN=np.linspace(3,8,41); NO=np.round(np.arange(11.5,19.0001,0.1),2); AA=np.round(np.arange(-1,1.5001,0.1),2)
class Fine:
    def __init__(s,geom='lvg'):
        d=np.load(f'/workspace/out/q3/grid_{geom}.npz'); T0=np.log10(d['TK']); N0=d['LN']; C0=d['LNC']
        def lnI(a):
            TR=a[...,0]; ok=(TR>0)&(a[...,2]>0)&np.isfinite(TR); return np.where(ok,np.log(1.0645*150*np.where(ok,TR,1)),np.nan),np.where(ok,a[...,2],np.nan)
        out={m:lnI(d[key]) for m,key in (('H','hcn'),('O','hcop'))}
        # interpolate each J separately
        TT,NN,CC=np.meshgrid(LT,LN,np.r_[NO],indexing='ij'); pts=np.stack([TT.ravel(),NN.ravel(),CC.ravel()],1)
        s.lnI={};s.tau={}
        for m in 'HO':
            l,t=out[m]; li=[];ti=[]
            for j in range(3):
                li.append(RegularGridInterpolator((T0,N0,C0),l[...,j],bounds_error=False,fill_value=np.nan)(pts).reshape(len(LT),len(LN),len(NO)))
                ti.append(RegularGridInterpolator((T0,N0,C0),np.log10(t[...,j]),bounds_error=False,fill_value=np.nan)(pts).reshape(len(LT),len(LN),len(NO)))
            s.lnI[m]=np.stack(li,-1); s.tau[m]=np.stack(ti,-1)       # (T,n,N,3)
        # pairs (iO,iA) with N_H=N_O+A inside the grid
        iO,iA=np.meshgrid(np.arange(len(NO)),np.arange(len(AA)),indexing='ij'); iH=iO+iA-10; ok=(iH>=0)&(iH<len(NO))
        s.iO,s.iA,s.iH=iO[ok],iA[ok],iH[ok]
        M=np.concatenate([s.lnI['H'][:,:,s.iH,:],s.lnI['O'][:,:,s.iO,:]],-1)   # (T,n,P,6): HCN2,3,4,HCO2,3,4
        s.shape=M.shape[:3]; M=M.reshape(-1,6); s.good=np.isfinite(M).all(1); s.M=M[s.good].astype(np.float64)
        ii=np.indices(s.shape); s.iT,s.iN,s.iP=[a.reshape(-1)[s.good] for a in ii]
        s.lT=LT[s.iT]; s.ln=LN[s.iN]; s.lNO=NO[s.iO[s.iP]]; s.lA=AA[s.iA[s.iP]]; s.lNH=NO[s.iH[s.iP]]
        s.tauH=s.tau['H'][s.iT,s.iN,s.iH[s.iP],:]; s.tauO=s.tau['O'][s.iT,s.iN,s.iO[s.iP],:]
class Engine:
    """Gaussian ln-flux likelihood, covariance Cm (6x6), common ln eta (flat, <=cap) marginalised analytically. Precomputes the model-dependent terms."""
    def __init__(s,F,Cm,cap=0.0):
        s.F=F; s.Ci=np.linalg.inv(Cm); s.a=float(np.ones(6)@s.Ci@np.ones(6)); s.u=s.Ci@np.ones(6)
        M=F.M; s.CiM=M@s.Ci; s.qMM=np.einsum('ij,ij->i',s.CiM,M); s.bM=M@s.u; s.cap=cap
    def lnlike(s,y):
        """y: ln fluxes (6,). returns lnL per model (marginalised over ln eta<=cap) and ln eta_hat"""
        yCy=y@s.Ci@y; yCiM=s.CiM@y; q=yCy-2*yCiM+s.qMM
        b=float(y@s.u)-s.bM; mu=b/s.a; chi=q-b*b/s.a
        # marginal over x=ln eta<=cap of exp(-0.5*(chi + a (x-mu)^2))
        lnL=-0.5*chi+log_ndtr(np.sqrt(s.a)*(s.cap-mu))
        return lnL,mu,chi
class Post:
    def __init__(s,F,lnw,mu=None):
        s.F=F; m=lnw.max(); s.p=np.exp(lnw-m); s.lnZ=m+np.log(s.p.sum())-np.log(len(lnw)); s.p/=s.p.sum(); s.mu=mu
    def wq(s,x,q=(2.5,16,50,84,97.5)):
        o=np.argsort(x); cs=np.cumsum(s.p[o]); xs=x[o]; return [float(xs[min(np.searchsorted(cs,qq/100),len(xs)-1)]) for qq in q]
    def summary(s):
        F=s.F; d={'log n':s.wq(F.ln),'Tk':s.wq(10**F.lT),'log N(HCO+)':s.wq(F.lNO),'log A':s.wq(F.lA),'log N(HCO+)/dv':s.wq(F.lNO-np.log10(150)),'log P/k':s.wq(F.ln+F.lT)}
        d['tau(HCN3-2)']=s.wq(F.tauH[:,1]); d['tau(HCO+3-2)']=s.wq(F.tauO[:,1])
        if s.mu is not None: d['log10 eta']=s.wq(s.mu/np.log(10))
        n=F.ln;T=10**F.lT;A=10**F.lA
        d['P']={'n>1e5':float(s.p[n>5].sum()),'n>1e6':float(s.p[n>6].sum()),'n<1e4':float(s.p[n<4].sum()),'Tk>50':float(s.p[T>50].sum()),'Tk<30':float(s.p[T<30].sum()),'Tk<20':float(s.p[T<20].sum()),'A>1':float(s.p[A>1].sum()),'A>3':float(s.p[A>3].sum()),
                'tauHCN3-2>1':float(s.p[10**F.tauH[:,1]>1].sum())}
        # lower/upper limit helpers: edges of prior
        d['edge']={'n_lo':float(s.p[F.ln<=3.0].sum()),'n_hi':float(s.p[F.ln>=8.0].sum()),'T_lo':float(s.p[F.lT<=1.0001].sum()),'T_hi':float(s.p[F.lT>=2.4771].sum()),'A_lo':float(s.p[F.lA<=-0.99].sum()),'A_hi':float(s.p[F.lA>=1.49].sum())}
        return d
    def map2d(s):
        m=np.zeros((len(LT),len(LN))); np.add.at(m,(s.F.iT,s.F.iN),s.p); return m
def covmat(rel,sig_b=(0.15,0.10,0.15),s_ind=0.07,extra=None):
    """6x6 ln-flux covariance: stat (rel, per datum) + s_ind per line + band-scale nuisance shared by the two molecules of a band; extra: additional per-datum sigma (6,)"""
    d=np.asarray(rel)**2+s_ind**2
    if extra is not None: d=d+np.asarray(extra)**2
    Cm=np.diag(d); band=np.array([0,1,2,0,1,2])
    for b,sb in enumerate(sig_b): v=(band==b).astype(float); Cm+=sb**2*np.outer(v,v)
    return Cm
def recovery_vec(Lf,delta,k47,k55):
    """multipliers f (6,) = [HCN 2,3,4, HCO+ 2,3,4]: band loss L_b=Lf*k_b (Xb55 k55, Xb4f 1, Xb47 k47); HCN recovers an extra (1-delta*k_b) relative to HCO+ (HCN more extended/resolved out)"""
    Lb=np.minimum(np.array([Lf*k55,Lf,Lf*k47]),0.8); db=np.array([delta*k55,delta,delta*k47])
    fO=1-Lb; fH=(1-Lb)*(1-db); return np.r_[fH,fO]
