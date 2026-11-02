# Q3 round 4 library. Extended RADEX grid (Tk 10-500 K, J=1-0..4-3; s62_radex_grid.py), likelihood engine with analytic ln(eta) marginalisation (eta<=1),
# and a Monte-Carlo marginalisation over *shared* nuisances (band scales, optional shared per-line errors, recovery factors) so that the same calibration is applied to several positions.
# Model grid: LT 35 pts (10-500 K, 0.05 dex), log n 3-8 (0.125), log N(HCO+) 11.5-19 (0.2), log A=log N(HCN)/N(HCO+) -1..1.4 (0.2); ln T_R interpolated linearly in (logT, log n, logN) from the RADEX grid.
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
from scipy.interpolate import RegularGridInterpolator
from scipy.special import log_ndtr,logsumexp
from scipy.stats import qmc
LT=np.linspace(1,np.log10(500),35); LN=np.linspace(3,8,41); NO=np.round(np.arange(11.5,19.0001,0.2),2); AA=np.round(np.arange(-1,1.4001,0.2),2)
class Grid:
    def __init__(s,step=0.1):
        global NO,AA
        NO=np.round(np.arange(11.5,19.0001,step),2); AA=np.round(np.arange(-1,1.4001,step),2); s.step=step
        d=np.load('/workspace/out/q3d/grid_lvg500.npz'); T0=np.log10(d['TK']); N0=d['LN']; C0=d['LNC']
        TT,NN,CC=np.meshgrid(LT,LN,NO,indexing='ij'); pts=np.stack([TT.ravel(),NN.ravel(),CC.ravel()],1)
        lnI={};tau={}
        for m,key in (('H','hcn'),('O','hcop')):
            a=d[key]; TR=a[...,0]; ok=(TR>0)&(a[...,2]>0)&np.isfinite(TR)
            l=np.where(ok,np.log(1.0645*150*np.where(ok,TR,1)),np.nan); t=np.where(ok,a[...,2],np.nan)
            li=[];ti=[]
            for j in range(4):
                li.append(RegularGridInterpolator((T0,N0,C0),l[...,j],bounds_error=False,fill_value=np.nan)(pts).reshape(len(LT),len(LN),len(NO)))
                ti.append(RegularGridInterpolator((T0,N0,C0),np.log10(t[...,j]),bounds_error=False,fill_value=np.nan)(pts).reshape(len(LT),len(LN),len(NO)))
            lnI[m]=np.stack(li,-1); tau[m]=np.stack(ti,-1)          # (T,n,N,4) J=1..4
        iO,iA=np.meshgrid(np.arange(len(NO)),np.arange(len(AA)),indexing='ij')
        # N_H = N_O + A  -> NO-grid index shift = (AA/0.2) since NO step 0.2 and AA step 0.2
        iH=iO+iA-int(round(1/s.step)); ok=(iH>=0)&(iH<len(NO)); s.iO,s.iA,s.iH=iO[ok],iA[ok],iH[ok]
        # 6 ladder columns (HCN 2,3,4, HCO+ 2,3,4) + J=1 columns
        M=np.concatenate([lnI['H'][:,:,s.iH,1:],lnI['O'][:,:,s.iO,1:]],-1); M1=np.stack([lnI['H'][:,:,s.iH,0],lnI['O'][:,:,s.iO,0]],-1)
        s.shape=M.shape[:3]; M=M.reshape(-1,6); M1=M1.reshape(-1,2); s.good=np.isfinite(M).all(1); s.good1=np.isfinite(M1).all(1)
        s.M=M[s.good].astype(np.float64); s.M1=M1[s.good]
        ii=np.indices(s.shape); s.iT,s.iN,s.iP=[a.reshape(-1)[s.good] for a in ii]
        s.lT=LT[s.iT]; s.ln=LN[s.iN]; s.lNO=NO[s.iO[s.iP]]; s.lA=AA[s.iA[s.iP]]; s.lNH=NO[s.iH[s.iP]]
        s.tauH=tau['H'][s.iT,s.iN,s.iH[s.iP],:]; s.tauO=tau['O'][s.iT,s.iN,s.iO[s.iP],:]
        s.nm=len(s.lT)
class Engine:
    """Gaussian ln-flux likelihood with diagonal/general covariance Cm (6x6), common ln eta (flat, <=cap) marginalised analytically."""
    def __init__(s,F,Cm,cap=0.0):
        s.F=F; s.Ci=np.linalg.inv(Cm); s.a=float(np.ones(6)@s.Ci@np.ones(6)); s.u=s.Ci@np.ones(6)
        M=F.M; s.CiM=M@s.Ci; s.qMM=np.einsum('ij,ij->i',s.CiM,M); s.bM=M@s.u; s.cap=cap
    def lnlike(s,y):
        yCy=y@s.Ci@y; q=yCy-2*(s.CiM@y)+s.qMM; b=float(y@s.u)-s.bM; mu=b/s.a; chi=q-b*b/s.a
        return -0.5*chi+log_ndtr(np.sqrt(s.a)*(s.cap-mu)),mu,chi
def diagcov(rel,s_ind=0.07,extra=None):
    d=np.asarray(rel)**2+s_ind**2
    if extra is not None: d=d+np.asarray(extra)**2
    return np.diag(d)
# ------------------------------------------------------------ priors (weights per grid model)
def prior_w(F,kind):
    T=10**F.lT
    if kind=='logT300': return (T<=300.01).astype(float)
    if kind=='logT500': return np.ones(F.nm)
    if kind=='linT500': return T/T.max()
    if kind=='linT300': return np.where(T<=300.01,T/300.,0.)
    if kind=='linn6': return np.where(T<=300.01,np.where(F.ln<=6.001,10**(F.ln-6),0.),0.)
    if kind=='linn7': return np.where(T<=300.01,np.where(F.ln<=7.001,10**(F.ln-7),0.),0.)
    if kind=='linn': return np.where(T<=300.01,10**(F.ln-8),0.)
    raise KeyError(kind)
# ------------------------------------------------------------ Monte-Carlo marginalisation over shared nuisances
def draws(nd,sig_b=(0.1,0.1,0.1),shared_line=0.0,rec_ranges=None,seed=1):
    """returns offsets (nd,6) [ln obs/true applied identically to all positions], and per-class recovery ln f (dict class->(nd,6))"""
    rng=qmc.Sobol(d=3+6+2*len(rec_ranges or {}),scramble=True,seed=seed); u=rng.random(nd)
    from scipy.stats import norm
    z=norm.ppf(np.clip(u,1e-9,1-1e-9)); band=z[:,:3]*np.array(sig_b); off=np.stack([band[:,[0,1,2,0,1,2]][:,k] for k in range(6)],1)
    off=off+z[:,3:9]*shared_line
    rec={}
    for k,(c,(r2,r4)) in enumerate((rec_ranges or {}).items()):
        f2=np.log(r2[0]+(r2[1]-r2[0])*u[:,9+2*k]); f4=np.log(r4[0]+(r4[1]-r4[0])*u[:,10+2*k]); rec[c]=np.stack([f2,0*f2,f4,f2,0*f2,f4],1)
    return off,rec
def mcfit(F,positions,nd=300,sig_b=(0.1,0.1,0.1),shared_line=0.0,s_ind=0.07,rec_ranges=None,common_theta=True,K=300,seed=1,keep_samples=True):
    """positions: list of dict(y=ln I(6), rel=stat(6), cls=recovery class or None, extra=optional extra sigma(6)).
    common_theta=True: one theta for all positions, per-position eta (exact truncated), shared nuisance draws. Returns lnpost (per model, ln mean_d prod_pos L), samples.
    common_theta=False: returns list of per-position (lnZ_d, lnL_d marginals handled separately) via mcfit_sep."""
    s_i=0.0 if shared_line>0 else s_ind
    En=[Engine(F,diagcov(p['rel'],s_i,p.get('extra'))) for p in positions]
    off,rec=draws(nd,sig_b,shared_line,rec_ranges,seed); rng=np.random.default_rng(seed)
    tot=None; S_idx=[];S_eta=[];S_w=[]
    for d in range(nd):
        lnL=0; mus=[]; 
        for En_,p in zip(En,positions):
            y=p['y']-off[d]+(rec[p['cls']][d] if (p.get('cls') in rec) else 0)
            l,mu,chi=En_.lnlike(y); lnL=lnL+l; mus.append((mu,En_.a))
        tot=lnL if tot is None else np.logaddexp(tot,lnL)
        if keep_samples:
            lz=logsumexp(lnL); p_=np.exp(lnL-lz); p_/=p_.sum(); idx=rng.choice(len(p_),K,p=p_)
            S_idx.append(idx); S_w.append(np.full(K,lz))
            e=np.empty((K,len(positions)))
            for j,(mu,a) in enumerate(mus):
                m=mu[idx]; sd=1/np.sqrt(a); ub=(0.0-m)/sd
                from scipy.stats import truncnorm
                e[:,j]=m+sd*truncnorm.rvs(-10,ub,size=K,random_state=rng.integers(1<<30)) if False else m+sd*_tn(ub,rng)
            S_eta.append(e)
    lnpost=tot-np.log(nd)
    out=dict(lnpost=lnpost)
    if keep_samples:
        w=np.concatenate(S_w); w=np.exp(w-w.max()); out.update(sidx=np.concatenate(S_idx),seta=np.concatenate(S_eta),sw=w/w.sum())
    return out
def _tn(ub,rng):
    """standard normal truncated to z<=ub (vector) by inverse cdf"""
    from scipy.stats import norm
    u=rng.random(len(ub)); return norm.ppf(np.clip(u*norm.cdf(ub),1e-300,1))
class Post:
    """posterior over grid models from lnpost with a prior-weight vector w (default flat log T<=300 K, flat log n)"""
    def __init__(s,F,lnpost,w,out=None):
        s.F=F; s.w=w; l=lnpost+np.log(np.where(w>0,w,1e-300)); l=np.where(w>0,l,-np.inf); m=l.max(); s.p=np.exp(l-m); s.lnZ=m+np.log(s.p.sum())-np.log(max((w>0).sum(),1)); s.p/=s.p.sum(); s.out=out
        s.pr=w/w.sum()
    def wq(s,x,q=(2.5,16,50,84,97.5),p=None):
        p=s.p if p is None else p; o=np.argsort(x); cs=np.cumsum(p[o]); xs=x[o]; return [float(xs[min(np.searchsorted(cs,qq/100),len(xs)-1)]) for qq in q]
    def prob(s,mask,prior=False): return float((s.pr if prior else s.p)[mask].sum())
    def sample_eta(s):
        o=s.out; 
        if o is None: return None
        # importance weights: prior reweight of samples drawn under the base (flat, all-T) prior
        wt=o['sw']*s.w[o['sidx']]; wt=wt/wt.sum(); return o['seta'],wt
    def summary(s):
        F=s.F; T=10**F.lT; d={'log n':s.wq(F.ln),'Tk':s.wq(T),'log A':s.wq(F.lA),'log N(HCO+)':s.wq(F.lNO),'log N(HCN)':s.wq(F.lNH),'log P/k':s.wq(F.ln+F.lT),
           'tau(HCN3-2)':s.wq(F.tauH[:,1]),'tau(HCO+3-2)':s.wq(F.tauO[:,1])}
        pr=lambda m:(s.prob(m),s.prob(m,True))
        d['P']={k:pr(m) for k,m in {'n>1e5':F.ln>5,'n>1e6':F.ln>6,'n<1e4':F.ln<4,'Tk>30':T>30,'Tk>50':T>50,'Tk>100':T>100,'Tk<30':T<30,'A>1':F.lA>0,'A>3':F.lA>np.log10(3),'tauHCN>1':F.tauH[:,1]>0,'tauHCN>5':F.tauH[:,1]>np.log10(5),'tauHCO>1':F.tauO[:,1]>0}.items()}
        e=s.sample_eta()
        if e is not None:
            seta,wt=e
            # eta ranges: sample index with weights (single-position eta = column 0)
            for j in range(seta.shape[1]):
                x=seta[:,j]/np.log(10); o=np.argsort(x); cs=np.cumsum(wt[o]); d[f'log10 eta[{j}]']=[float(x[o][min(np.searchsorted(cs,q/100),len(x)-1)]) for q in (2.5,16,50,84,97.5)]
        return d
    def marg1d(s,x,bins):
        return np.histogram(x,bins=bins,weights=s.p)[0],np.histogram(x,bins=bins,weights=s.pr)[0]
def logmean(a): return logsumexp(a)-np.log(len(a))

# ------------------------------------------------------------ analytic joint engine (Gaussian nuisances marginalised analytically; per-position ln eta flat, <=cap)
REC_CLASSES={'centre':((.85,.95),(1.03,1.12)),'E1':((.95,1.05),(.95,1.05)),'E2':((.6,1.),(1.,1.5)),'W':((.6,1.),(1.,1.5)),'none':((1.,1.),(1.,1.))}
def rec_stats(cls):
    """mean and sd of ln f for the uniform recovery ranges of class cls: vectors (6,) for [HCN2,3,4,HCO2,3,4]"""
    (a2,b2),(a4,b4)=REC_CLASSES[cls]; u=np.linspace(0,1,2001)
    l2=np.log(a2+(b2-a2)*u); l4=np.log(a4+(b4-a4)*u)
    return np.array([l2.mean(),0,l4.mean(),l2.mean(),0,l4.mean()]),np.array([l2.std(),0,l4.std(),l2.std(),0,l4.std()])
class GEngine:
    """positions k=0..n-1 each with 6 data. Covariance: diag(stat^2 + s_ind^2 + extra^2); shared band scales (sig_b) over all positions; shared per-line errors s_line (cross-position);
    recovery (mean shift + sd, correlated over positions of the same class and over the two molecules for the same J); per-position ln eta flat <=cap marginalised:
    Gaussian conditional with product-of-marginals truncation (exact for n=1)."""
    def __init__(s,F,positions,sig_b=(.1,.1,.1),s_ind=0.07,s_line=0.0,use_rec=True,cap=0.0):
        s.F=F; n=len(positions); s.n=n; N6=6*n; C=np.zeros((N6,N6)); band=np.array([0,1,2,0,1,2])
        for k,p in enumerate(positions):
            d=np.asarray(p['rel'])**2+s_ind**2+(np.asarray(p['extra'])**2 if p.get('extra') is not None else 0.)
            C[6*k:6*k+6,6*k:6*k+6]+=np.diag(d)
        for b,sb in enumerate(sig_b):
            v=np.tile((band==b).astype(float),n); C+=sb**2*np.outer(v,v)
        if s_line>0:
            for j in range(6):
                v=np.zeros(N6); v[j::6]=1; C+=s_line**2*np.outer(v,v)
        shift=np.zeros(N6)
        if use_rec:
            for c in sorted(set(p.get('cls') for p in positions if p.get('cls'))):
                mu,sd=rec_stats(c)
                for cols in ((0,3),(2,5)):
                    vv=np.zeros(N6)
                    for k,p in enumerate(positions):
                        if p.get('cls')==c:
                            for j in cols: vv[6*k+j]=sd[j]; shift[6*k+j]+=mu[j]
                    C+=np.outer(vv,vv)
        s.y=np.concatenate([p['y'] for p in positions])+shift; s.Ci=np.linalg.inv(C); s.cap=cap
        Ci=s.Ci; s.A=np.array([[np.ones(6)@Ci[6*k:6*k+6,6*l:6*l+6]@np.ones(6) for l in range(n)] for k in range(n)]); s.Ai=np.linalg.inv(s.A)
        M=F.M; s.MC={(k,l):M@Ci[6*k:6*k+6,6*l:6*l+6].T for k in range(n) for l in range(n)}
        s.yc={k:[Ci[6*k:6*k+6,6*l:6*l+6]@s.y[6*l:6*l+6] for l in range(n)] for k in range(n)}
    def lnlike(s):
        n=s.n; g=[sum(s.yc[k][l][None,:]-s.MC[(k,l)] for l in range(n)) for k in range(n)]
        q0=0; bb=np.empty((s.F.nm,n))
        for k in range(n):
            r=s.y[6*k:6*k+6][None,:]-s.F.M; q0=q0+np.einsum('ij,ij->i',r,g[k]); bb[:,k]=g[k].sum(1)
        eh=bb@s.Ai.T; chi=q0-np.einsum('ij,ij->i',bb@s.Ai,bb)
        sd=np.sqrt(np.diag(s.Ai)); lt=log_ndtr((s.cap-eh)/sd[None,:]).sum(1)
        return -0.5*chi+lt,eh,sd,chi
def gfit(F,positions,keep_samples=True,K=40000,seed=1,**kw):
    E=GEngine(F,positions,**kw); l,eh,sd,chi=E.lnlike(); out=dict(lnpost=l,chi=chi,eh=eh,sd=sd)
    if keep_samples:
        from scipy.stats import norm
        rng=np.random.default_rng(seed); p=np.exp(l-l.max()); p/=p.sum(); idx=rng.choice(len(p),K,p=p); e=np.empty((K,E.n))
        for k in range(E.n):
            m=eh[idx,k]; ub=(E.cap-m)/sd[k]; u=rng.random(K); e[:,k]=m+sd[k]*norm.ppf(np.clip(u*norm.cdf(ub),1e-300,1))
        out.update(sidx=idx,seta=e,sw=np.full(K,1.0/K))
    return out
