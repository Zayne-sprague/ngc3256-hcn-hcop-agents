# Q3: non-LTE (RADEX) fit of the HCN and HCO+ J=2-1,3-2,4-3 ladders at a common 10" beam; posteriors, region comparison, systematics.
# Data: out/s30_ratios.json (Xb55 J=2 / Xb4f J=3 / Xb47 J=4 smoothed to 10" circular, |v-vsys|<350 km/s, Jy/beam km/s) -> K km/s (RJ, rest frame).
#       zone (deblended) fluxes from out/s37_deblend2.json x Xb4d 3-2 template fluxes.
# Model: single RADEX zone (grid out/q3/grid_<geom>.npz from s52): Tkin, n(H2), N(HCO+), N(HCN) [per 150 km/s], common beam/velocity filling eta (nuisance, flat in ln eta, eta<=1).
#   I_model(J,mol) = eta * 1.0645*150 km/s * T_R(J,mol).   Flat priors in log n (3-8), log Tk (10-300), log N(HCO+) (11.5-19), log A=log N(HCN)/N(HCO+) (-1..1.5).
# Likelihood: Gaussian in ln I, covariance = diag(stat^2 + s_ind^2) + sum_band s_b^2 (1_b 1_b^T)  (band-scale nuisance shared by HCN and HCO+ and, in joint fits, by all regions).
# usage: python3 -I s53_ladder_fit.py   -> out/q3/*.txt, *.json, *.npz, *.png
import json,sys,itertools,numpy as np
from scipy.special import logsumexp
from scipy.stats import chi2 as chi2d
OUT='/workspace/out/q3/'
C=299792.458; VSYS=2750.; Z=VSYS/C
REST={'HCN':{2:177.26111,3:265.88618,4:354.50547},'HCO+':{2:178.37506,3:267.55763,4:356.73424}}
JS=(2,3,4); NU=lambda m,J:REST[m][J]/(1+Z)
def tk(S,m,J,omega):   # Jy km/s over solid angle omega [arcsec^2] -> K km/s, RJ, rest-frame (T_rest=(1+Z) T_obs)
    return 1.222e6*S/(NU(m,J)**2*omega)*(1+Z)
OM10=1.1331*100.
# ---------------------------------------------------------------- data
R30=json.load(open('/workspace/out/s30_ratios.json')); R30={r['tag']:r for r in R30}
TAGS={2:'Xb55_J2_10as',3:'Xb4f_J3_10as',4:'Xb47_J4_10as'}
def region_data(row,ap=None):
    I=[];e=[]
    for m,k,ek in (('HCN','HCN','eHCN'),('HCO+','HCOp','eHCOp')):
        for J in JS:
            r=R30[TAGS[J]]['rows'][row]; om=OM10 if ap is None else ap
            I.append(tk(r[k],m,J,om)); e.append(tk(r[ek],m,J,om)/tk(r[k],m,J,om))
    return np.array(I),np.array(e)     # order HCN2,3,4, HCO2,3,4 ; e = relative stat error
ROWS={'N':'N_pos','S':'S_pos','Bridge':'Bridge_pos','E1':'E1_pos','E2':'E2_pos','W':'W_pos','Sblue':'Sblue_pos'}
DATA={k:region_data(v) for k,v in ROWS.items()}
DATA['ap<8"']=region_data('ap_r<8',ap=np.pi*64.)
DATA['ap<4"']=region_data('ap_r<4',ap=np.pi*16.)
# disk = mean of the three disk positions (E1,E2,W) (independent noise, calibration shared)
I_=np.mean([DATA[k][0] for k in ('E1','E2','W')],0); DATA['disk(E1,E2,W)']=(I_,np.full(6,0.01))
# zones (deblended N,S) from s37: flux_J = f_J * template(zone,mol) [Jy km/s, total], same solid angle for all J
D37=json.load(open('/workspace/out/s37_deblend2.json'))['nominal']
TEMPL={'N':{'HCN':22.0,'HCO+':39.5},'S':{'HCN':15.3,'HCO+':40.8}}
ZT={2:'Xb55_J2_10',3:'Xb4f_J3_10',4:'Xb47_J4_10'}
for zi,z in ((0,'N'),(1,'S')):
    I=[];e=[]
    for m,key in (('HCN','HCN'),('HCO+','HCO+')):
        for J in JS:
            f=D37[ZT[J]][key]['f'][zi]; ef=D37[ZT[J]][key]['ef'][zi]
            I.append(tk(f*TEMPL[z][m],m,J,OM10)); e.append(ef/f)
    DATA['zoneN' if z=='N' else 'zoneS']=(np.array(I),np.array(e))
# N core subtraction (two-zone): core flux r<0.5" Xb4b 3-2 (s43), added to the 10" N beam at J=3, ladder L=I_K(J)/I_K(3)
import pandas as pd
b=pd.read_csv('/workspace/out/s43_profiles.csv'); n=b[b.region=='N'].reset_index(drop=True); bm=1.1331*0.2456*0.2373/0.047**2
CORE_JY={'HCN':float((n.I_HCN*n.npix/bm)[:2].sum()),'HCO+':float((n.I_HCO*n.npix/bm)[:2].sum())}
def two_zone(Lc43,Lc23,relerr=0.15):
    I,e=DATA['N']; I=I.copy(); e=e.copy()
    for i,(m,J) in enumerate([(m,J) for m in ('HCN','HCO+') for J in JS]):
        c3=tk(CORE_JY[m],m,3,OM10); c=c3*{2:Lc23,3:1.0,4:Lc43}[J]
        env=I[i]-c; err=np.hypot(I[i]*e[i],c*relerr); I[i]=env; e[i]=err/env
    return I,e
# ---------------------------------------------------------------- grid
def load_grid(geom):
    d=np.load(OUT+f'grid_{geom}.npz'); TK,LN,LNC=d['TK'],d['LN'],d['LNC']
    def lnI(a):   # (T,n,N,J,3) -> ln(1.0645*150*T_R) with T_R>0 and tau>0
        TR=a[...,0][...,:]; ok=(TR>0)&(a[...,2]>0)
        with np.errstate(all='ignore'): return np.where(ok,np.log(1.0645*150*np.where(ok,TR,1)),np.nan),a[...,1],a[...,2]
    H=lnI(d['hcn']); O=lnI(d['hcop'])
    return TK,LN,LNC,H,O
class Model:
    def __init__(s,geom='lvg',coarse=False):
        TK,LN,LNC,H,O=load_grid(geom); sl=(slice(None),slice(None,None,2),slice(None,None,2)) if coarse else (slice(None),)*3
        s.TK,s.LN,s.LNC=TK,LN[sl[1]],LNC[sl[2]]; hs=(slice(None),sl[1],sl[2])
        lH=H[0][hs];lO=O[0][hs]                                  # (T,n,N,J)
        s.lH,s.lO=lH,lO; s.TexH,s.tauH=H[1][hs],H[2][hs]; s.TexO,s.tauO=O[1][hs],O[2][hs]
        T,nn,N=len(TK),len(s.LN),len(s.LNC)
        # model ln I[T,n,NO,NH,6]
        s.shape=(T,nn,N,N)
        mod=np.empty((T,nn,N,N,6));
        for k in range(3): mod[...,k]=lH[:,:,None,:,k]; mod[...,3+k]=lO[:,:,:,None,k]
        s.mod=mod.reshape(-1,6); s.good=np.isfinite(s.mod).all(1)
        NO,NH=np.meshgrid(s.LNC,s.LNC,indexing='ij'); lA=(NH-NO)
        s.lA=np.broadcast_to(lA,(T,nn,N,N)).reshape(-1); s.inprior=(s.lA>=-1.0)&(s.lA<=1.5)
        s.w=(s.good&s.inprior)
        ii=np.indices(s.shape); s.iT,s.iN,s.iO,s.iH=[a.reshape(-1) for a in ii]
def cov(rel,sig_b,s_ind=0.07):
    """6x6 (or 12x12 stacked for 2 regions with rel=(rel1,rel2)) covariance in ln I"""
    rel=np.atleast_2d(rel); nr=rel.shape[0]; Cm=np.diag((rel**2+s_ind**2).ravel())
    band=np.tile([0,1,2,0,1,2],nr)
    for b,sb in enumerate(sig_b):
        v=(band==b).astype(float); Cm+=sb**2*np.outer(v,v)
    return Cm
SCEN={'B7free':(0.20,0.10,10.),'B5free':(10.,0.10,0.20),'cal10':(0.10,0.10,0.10),'cal20':(0.20,0.10,0.20),'cal35':(0.35,0.15,0.35),'ratios-only':(10.,10.,10.)}
LNETA_MAX=0.0
def marg(M,y,Cm,nreg=1):
    """profile over per-region ln eta (flat) analytically. y:(6*nreg), model ln I M:(Ng,6) (same theta for all regions). returns chi2_marg, ln eta_hat (Ng,nreg)"""
    Ci=np.linalg.inv(Cm); D=np.zeros((6*nreg,nreg))
    for k in range(nreg): D[6*k:6*k+6,k]=1
    Mm=D.T@Ci@D; Mi=np.linalg.inv(Mm)
    r=y[None,:]-np.tile(M,(1,nreg)); Cr=r@Ci
    q=np.einsum('ij,ij->i',Cr,r); bb=Cr@D; mu=bb@Mi.T
    return q-np.einsum('ij,ij->i',bb@Mi,bb),mu
def post_single(mod,I,rel,scn,tau_cut=None,X=None):
    y=np.log(I); Cm=cov(rel,SCEN[scn]); w=mod.w.copy()
    c2,mu=marg(mod.mod[w],y,Cm); mu=mu[:,0]
    ok=mu<=LNETA_MAX
    idx=np.where(w)[0][ok]; c2=c2[ok]; mu=mu[ok]
    if tau_cut is not None:   # H13CN limit on tau(HCN 4-3) of the model zone
        tH4=mod.tauH[mod.iT[idx],mod.iN[idx],mod.iH[idx],2]; keep=tH4<tau_cut; idx=idx[keep]; c2=c2[keep]; mu=mu[keep]
    lw=-0.5*c2; Zev=logsumexp(lw)-np.log(w.sum())
    p=np.exp(lw-lw.max()); p/=p.sum()
    return dict(idx=idx,p=p,c2=c2,mu=mu,lnZ=Zev,ndof_eff=6-5)
def quant(mod,P,q=(16,50,84)):
    idx=P['idx'];p=P['p']
    def wq(x,qs=q):
        o=np.argsort(x); cs=np.cumsum(p[o]); return [float(x[o][np.searchsorted(cs,qq/100)]) for qq in qs]
    out={}
    out['log n']=wq(mod.LN[mod.iN[idx]]); out['Tk']=wq(mod.TK[mod.iT[idx]])
    out['log N(HCO+)']=wq(mod.LNC[mod.iO[idx]]); out['log N(HCN)']=wq(mod.LNC[mod.iH[idx]]); out['log A']=wq(mod.lA[idx])
    out['log10 eta']=[v/np.log(10) for v in wq(P['mu'])]
    out['tau(HCN3-2)']=wq(np.log10(mod.tauH[mod.iT[idx],mod.iN[idx],mod.iH[idx],1]))
    out['tau(HCO+3-2)']=wq(np.log10(mod.tauO[mod.iT[idx],mod.iN[idx],mod.iO[idx],1]))
    out['log P/k']=wq(mod.LN[mod.iN[idx]]+np.log10(mod.TK[mod.iT[idx]]))
    out['log N(HCO+)/dv']=[v-np.log10(150) for v in out['log N(HCO+)']]
    return out
def probs(mod,P):
    idx=P['idx'];p=P['p']; n=mod.LN[mod.iN[idx]]; T=mod.TK[mod.iT[idx]]; A=10**mod.lA[idx]
    tH=mod.tauH[mod.iT[idx],mod.iN[idx],mod.iH[idx],1]
    return {'P(n>1e5)':float(p[n>5].sum()),'P(n>1e6)':float(p[n>6].sum()),'P(n<1e4)':float(p[n<4].sum()),'P(Tk>50)':float(p[T>50].sum()),'P(Tk<30)':float(p[T<30].sum()),
            'P(A>1)':float(p[A>1].sum()),'P(A>3)':float(p[A>3].sum()),'P(tauHCN3-2>1)':float(p[tH>1].sum()),'P(tauHCN3-2>10)':float(p[tH>10].sum())}
def nT_map(mod,P):
    m=np.zeros((len(mod.TK),len(mod.LN))); np.add.at(m,(mod.iT[P['idx']],mod.iN[P['idx']]),P['p']); return m
# ---------------------------------------------------------------- pair comparison (N vs S etc.) on the coarse grid
def pair_compare(modc,A,B,scn):
    (IA,eA),(IB,eB)=A,B; y=np.log(np.r_[IA,IB]); Cm=cov(np.array([eA,eB]),SCEN[scn]); Ci=np.linalg.inv(Cm)
    D=np.zeros((12,2)); D[:6,0]=1; D[6:,1]=1; Mi=np.linalg.inv(D.T@Ci@D); Pm=Ci-Ci@D@Mi@D.T@Ci   # 12x12 projected precision
    K=Mi@D.T@Ci                                                                                    # eta_hat = K r (2x12)
    w=modc.w; idx=np.where(w)[0]; M=modc.mod[idx]; r=[y[:6][None,:]-M, y[6:][None,:]-M]
    PAA,PAB,PBB=Pm[:6,:6],Pm[:6,6:],Pm[6:,6:]
    qA=np.einsum('ij,ij->i',r[0]@PAA,r[0]); qB=np.einsum('ij,ij->i',r[1]@PBB,r[1])
    muA_a=r[0]@K[0,:6]; muA_b=r[1]@K[0,6:]; muB_a=r[0]@K[1,:6]; muB_b=r[1]@K[1,6:]
    Ng=len(idx); iT=modc.iT[idx]; iN=modc.iN[idx]
    # H0: same theta
    q0=qA+qB+2*np.einsum('ij,ij->i',r[0]@PAB,r[1]); ok0=(muA_a+muA_b<=0)&(muB_a+muB_b<=0)
    lnZ0=logsumexp(-0.5*q0[ok0])-np.log(Ng); chi_min0=q0[ok0].min(); th0=idx[ok0][np.argmin(q0[ok0])]
    # H1: same (T,n), separate (N,A); H2: all separate. full pair matrices in chunks
    lw1=[];lw2=[];chi2min=1e9; chi1min=1e9
    cell=iT*1000+iN
    chunk=400; tot2=[]; blocks={}
    cells=np.unique(cell); cs={c:np.where(cell==c)[0] for c in cells}
    for s in range(0,Ng,chunk):
        sl=slice(s,s+chunk)
        cross=r[0][sl]@PAB@r[1].T                                  # (c,Ng)
        q=qA[sl,None]+qB[None,:]+2*cross
        ok=((muA_a[sl,None]+muA_b[None,:])<=0)&((muB_a[sl,None]+muB_b[None,:])<=0)
        q=np.where(ok,q,np.inf); chi2min=min(chi2min,q.min())
        tot2.append(logsumexp(-0.5*q))
        same=(cell[sl,None]==cell[None,:]); q1=np.where(same,q,np.inf)
        chi1min=min(chi1min,q1.min())
        # per cell sums for H1
        for c in np.unique(cell[sl]):
            rows=np.where(cell[sl]==c)[0]; lw1.append((c,logsumexp(-0.5*q1[rows])))
    lnZ2=logsumexp(tot2)-2*np.log(Ng)
    d={}
    for c,v in lw1: d.setdefault(c,[]).append(v)
    # H1: for each (T,n) cell, sum over pairs within the cell / Ncell^2, averaged over cells
    ncell={c:len(cs[c]) for c in cs}
    lnZ1=logsumexp([logsumexp(v)-2*np.log(ncell[c]) for c,v in d.items()])-np.log(len(cs))
    return dict(lnZ0_shared=float(lnZ0),lnZ1_sharedTn=float(lnZ1),lnZ2_separate=float(lnZ2),lnB_20=float(lnZ2-lnZ0),lnB_10=float(lnZ1-lnZ0),lnB_21=float(lnZ2-lnZ1),
                chi2min_shared=float(chi_min0),chi2min_sharedTn=float(chi1min),chi2min_separate=float(chi2min),
                dchi2_shared_vs_sep=float(chi_min0-chi2min),dchi2_sharedTn_vs_sep=float(chi1min-chi2min))
# ---------------------------------------------------------------- main
if __name__=='__main__':
    log=[];res={}
    def P(s=''): print(s); log.append(s)
    mod=Model('lvg'); modc=Model('lvg',coarse=True)
    P(f'grid {mod.shape} models; usable (finite RADEX, A in prior): {mod.w.sum()}/{len(mod.w)} ({100*mod.w.mean():.0f}%)')
    P(f'N core flux (Xb4b r<0.5", |v|<350): HCN {CORE_JY["HCN"]:.2f}, HCO+ {CORE_JY["HCO+"]:.2f} Jy km/s')
    P('\n# Observed ladders at the 10" beam: K km/s (RJ, rest frame) [HCN 2-1,3-2,4-3 | HCO+ 2-1,3-2,4-3] and R_K=HCN/HCO+')
    for k,(I,e) in DATA.items():
        P(f'{k:16s} '+' '.join(f'{v:7.2f}' for v in I)+f' | R_K {I[0]/I[3]:.3f} {I[1]/I[4]:.3f} {I[2]/I[5]:.3f}  | T(J)/T(3-2): HCN {I[0]/I[1]:.2f},{I[2]/I[1]:.2f}  HCO+ {I[3]/I[4]:.2f},{I[5]/I[4]:.2f}')
    regions=['N','S','Bridge','disk(E1,E2,W)','Sblue','ap<8"','zoneN','zoneS']
    # (1) single-region fits, scenarios
    P('\n# Posteriors (median [16,84]) single zone LVG; flat priors; scenarios of band-to-band calibration (sigma per band 2-1,3-2,4-3)')
    store={}
    for scn in ('cal20','cal10','cal35','ratios-only'):
        for rg in regions:
            I,e=DATA[rg]; Pp=post_single(mod,I,e,scn); q=quant(mod,Pp); pr=probs(mod,Pp)
            # best fit: min chi2 of marginal
            ib=np.argmin(Pp['c2']); idb=Pp['idx'][ib]
            store[(scn,rg)]=(Pp,q,pr)
            res[f'{scn}|{rg}']=dict(q=q,probs=pr,chi2min=float(Pp['c2'].min()),lnZ=Pp['lnZ'],best=dict(Tk=float(mod.TK[mod.iT[idb]]),logn=float(mod.LN[mod.iN[idb]]),logNO=float(mod.LNC[mod.iO[idb]]),logA=float(mod.lA[idb]),log10eta=float(Pp['mu'][ib]/np.log(10))))
            P(f'[{scn:11s}] {rg:15s} chi2min={Pp["c2"].min():5.2f} (dof<=1) log n {q["log n"][1]:.2f} [{q["log n"][0]:.2f},{q["log n"][2]:.2f}]  Tk {q["Tk"][1]:.0f} [{q["Tk"][0]:.0f},{q["Tk"][2]:.0f}]  logN(HCO+) {q["log N(HCO+)"][1]:.2f} [{q["log N(HCO+)"][0]:.2f},{q["log N(HCO+)"][2]:.2f}]  logA {q["log A"][1]:.2f} [{q["log A"][0]:.2f},{q["log A"][2]:.2f}]  log10 eta {q["log10 eta"][1]:.2f}  logP/k {q["log P/k"][1]:.1f} [{q["log P/k"][0]:.1f},{q["log P/k"][2]:.1f}] | P(n>1e5)={pr["P(n>1e5)"]:.2f} P(Tk>50)={pr["P(Tk>50)"]:.2f} P(A>1)={pr["P(A>1)"]:.2f} P(tau3>1)={pr["P(tauHCN3-2>1)"]:.2f}')
        P('')
    # (2) prior-vs-posterior information: marginal posterior of n and T for N (cal20) as 2D maps
    maps={f'{rg}':nT_map(mod,store[('cal20',rg)][0]) for rg in regions}
    np.savez(OUT+'s53_nT_posteriors.npz',TK=mod.TK,LN=mod.LN,**{k.replace('"','').replace('<','lt').replace('(','_').replace(')','').replace(',','_'):v for k,v in maps.items()})

    # (2b) extra calibration scenarios: Band 7 or Band 5 absolute scale free (only the other bands + the same-J ratios constrain)
    P('# Band-7 (4-3) or Band-5 (2-1) absolute scale left free (sigma=1000%)')
    for scn in ('B7free','B5free'):
        for rg in ('N','S','Bridge','disk(E1,E2,W)'):
            I,e=DATA[rg]; Pp=post_single(mod,I,e,scn); q=quant(mod,Pp); pr=probs(mod,Pp); res[f'{scn}|{rg}']=dict(q=q,probs=pr,chi2min=float(Pp['c2'].min()))
            P(f'[{scn}] {rg:15s} chi2min={Pp["c2"].min():.2f} log n {q["log n"][1]:.2f} [{q["log n"][0]:.2f},{q["log n"][2]:.2f}] Tk {q["Tk"][1]:.0f} [{q["Tk"][0]:.0f},{q["Tk"][2]:.0f}] logA {q["log A"][1]:.2f} [{q["log A"][0]:.2f},{q["log A"][2]:.2f}] | P(n>1e5)={pr["P(n>1e5)"]:.2f} P(Tk>50)={pr["P(Tk>50)"]:.2f} P(A>1)={pr["P(A>1)"]:.2f}')
    # (2c) profile chi2 (min over the other parameters) vs n and vs Tk, and the allowed (Tk,n) fraction (chi2<chi2min+4)
    P('\n# Profile chi2 vs n(H2) and vs Tk (min over the others), cal20; allowed cell fraction = (Tk,n) cells with profile chi2 < min+4')
    for scn in ('cal20','ratios-only'):
        for rg in ('N','S','disk(E1,E2,W)'):
            Pp=store[(scn,rg)][0]; idx=Pp['idx']; c2=Pp['c2']; iT=mod.iT[idx]; iN=mod.iN[idx]
            prof=np.full((len(mod.TK),len(mod.LN)),np.inf); np.minimum.at(prof,(iT,iN),c2)
            m0=prof.min(); allowed=(prof<m0+4)
            pn=[(float(mod.LN[j]),float(prof[:,j].min()-m0),float(mod.TK[prof[:,j].argmin()])) for j in range(0,len(mod.LN),4)]
            pT=[(float(mod.TK[i]),float(prof[i].min()-m0),float(mod.LN[prof[i].argmin()])) for i in range(len(mod.TK))]
            res[f'profile|{scn}|{rg}']=dict(allowed_frac=float(allowed.mean()),prof_n=pn,prof_T=pT,chi2min=float(m0))
            P(f'[{scn}] {rg:15s} chi2min {m0:.2f}; allowed (Tk,n) cells {100*allowed.mean():.0f}% of {allowed.size}')
            P('    d chi2 vs log n (best Tk): '+'  '.join(f'{a:.1f}:{b:.1f}(T{c:.0f})' for a,b,c in pn))
            P('    d chi2 vs Tk (best log n): '+'  '.join(f'{a:.0f}:{b:.1f}(n{c:.1f})' for a,b,c in pT))
            # A at fixed n: best A
            Ab=[]
            for j in range(0,len(mod.LN),4):
                sel=iN==j
                if sel.any(): ib=np.argmin(np.where(sel,c2,np.inf)); Ab.append((float(mod.LN[j]),float(mod.lA[idx[ib]]),float(mod.LNC[mod.iO[idx[ib]]])))
            P('    best log A, log N(HCO+) at fixed log n: '+'  '.join(f'{a:.1f}:A{b:+.2f},N{c:.2f}' for a,b,c in Ab))
    # derived: N(H2), path length, for assumed X(HCO+)
    P('\n# Derived with assumed X(HCO+)=1e-9 (factor 3-10 uncertain): N(H2)=N(HCO+)/X per 150 km/s and path L=N(H2)/n(H2) [pc], cal20')
    for rg in ('N','S','Bridge','disk(E1,E2,W)'):
        Pp=store[('cal20',rg)][0]; idx=Pp['idx']; p=Pp['p']; NH2=mod.LNC[mod.iO[idx]]+9.; L=NH2-mod.LN[mod.iN[idx]]-np.log10(3.0857e18)
        def wq(x):
            o=np.argsort(x); cs=np.cumsum(p[o]); return [round(float(x[o][np.searchsorted(cs,qq/100)]),2) for qq in (16,50,84)]
        res[f'derived|{rg}']=dict(logNH2=wq(NH2),logL_pc=wq(L)); P(f'  {rg:15s} log N(H2)/150km/s {wq(NH2)}  log L[pc] {wq(L)}')
    # (3) H13CN opacity prior: tau(HCN 4-3) < 3.8 / 6.4 / 12.8 (X=30,50,100)  (Q4, beam-averaged limit at N; applied to the model zone for every region)
    P('# With the Q4 H13CN limit on tau(HCN 4-3) of the zone (X=12C/13C=30,50,100 -> <3.8,6.4,12.8), cal20')
    for tc,X in ((3.8,30),(6.4,50),(12.8,100)):
        for rg in ('N','S','Bridge','disk(E1,E2,W)'):
            I,e=DATA[rg]; Pp=post_single(mod,I,e,'cal20',tau_cut=tc); q=quant(mod,Pp); pr=probs(mod,Pp)
            res[f'tau<{tc}|{rg}']=dict(q=q,probs=pr,lnZ=Pp['lnZ'])
            P(f'  X={X:3d} {rg:15s} log n {q["log n"][1]:.2f} [{q["log n"][0]:.2f},{q["log n"][2]:.2f}] Tk {q["Tk"][1]:.0f} [{q["Tk"][0]:.0f},{q["Tk"][2]:.0f}] logA {q["log A"][1]:.2f} [{q["log A"][0]:.2f},{q["log A"][2]:.2f}] logNHCO {q["log N(HCO+)"][1]:.2f} | P(n>1e5)={pr["P(n>1e5)"]:.2f} P(A>1)={pr["P(A>1)"]:.2f}')
    # (4) geometry: slab instead of LVG
    P('\n# Geometry: plane-parallel slab (radex_slab) instead of expanding sphere, cal20')
    mods=Model('slab')
    for rg in ('N','S','Bridge','disk(E1,E2,W)'):
        I,e=DATA[rg]; Pp=post_single(mods,I,e,'cal20'); q=quant(mods,Pp); pr=probs(mods,Pp)
        res[f'slab|{rg}']=dict(q=q,probs=pr)
        P(f'  {rg:15s} chi2min={Pp["c2"].min():.2f} log n {q["log n"][1]:.2f} [{q["log n"][0]:.2f},{q["log n"][2]:.2f}] Tk {q["Tk"][1]:.0f} [{q["Tk"][0]:.0f},{q["Tk"][2]:.0f}] logA {q["log A"][1]:.2f} [{q["log A"][0]:.2f},{q["log A"][2]:.2f}] | P(n>1e5)={pr["P(n>1e5)"]:.2f} P(A>1)={pr["P(A>1)"]:.2f}')
    # (5) two-zone: subtract the N core (thick/mid/thin ladders) before fitting the envelope
    P('\n# N envelope after subtracting the compact core (Xb4b r<0.5") at the 10" beam; core K-ladder L(4-3/3-2),L(2-1/3-2)')
    for nm,(L43,L23) in {'thick':(1.0,1.0),'mid':(0.8,1.1),'thin':(0.57,1.3)}.items():
        I,e=two_zone(L43,L23); Pp=post_single(mod,I,e,'cal20'); q=quant(mod,Pp); pr=probs(mod,Pp)
        res[f'N_env_{nm}']=dict(q=q,probs=pr,I=list(I),R=[I[0]/I[3],I[1]/I[4],I[2]/I[5]])
        P(f'  core {nm:5s}: env R_K {I[0]/I[3]:.2f} {I[1]/I[4]:.2f} {I[2]/I[5]:.2f}  log n {q["log n"][1]:.2f} [{q["log n"][0]:.2f},{q["log n"][2]:.2f}] Tk {q["Tk"][1]:.0f} [{q["Tk"][0]:.0f},{q["Tk"][2]:.0f}] logA {q["log A"][1]:.2f} [{q["log A"][0]:.2f},{q["log A"][2]:.2f}] | P(n>1e5)={pr["P(n>1e5)"]:.2f} P(A>1)={pr["P(A>1)"]:.2f}')
    # (6) region comparison
    P('\n# Region-to-region comparison: shared (n,Tk,N,A) [H0] vs shared (n,Tk) with separate (N,A) [H1] vs all separate [H2]; band-scale nuisance common to both regions; eta separate; coarse grid')
    pairs=[('N','S'),('zoneN','zoneS'),('N','Bridge'),('S','Bridge'),('Bridge','disk(E1,E2,W)'),('N','disk(E1,E2,W)'),('S','disk(E1,E2,W)')]
    for scn in ('cal20','cal10','ratios-only'):
        for a,bb in pairs:
            if scn!='cal20' and (a,bb) not in (('N','S'),('zoneN','zoneS'),('Bridge','disk(E1,E2,W)')): continue
            d=pair_compare(modc,DATA[a],DATA[bb],scn); res[f'pair|{scn}|{a}|{bb}']=d
            P(f'[{scn:11s}] {a:8s} vs {bb:15s} lnB(sep/shared)={d["lnB_20"]:6.2f}  lnB(sep/sharedTn)={d["lnB_21"]:6.2f}  dchi2(shared-sep)={d["dchi2_shared_vs_sep"]:6.2f} (4 extra pars; p={chi2d.sf(max(d["dchi2_shared_vs_sep"],0),4):.3f})  dchi2(sharedTn-sep)={d["dchi2_sharedTn_vs_sep"]:6.2f} (2 extra; p={chi2d.sf(max(d["dchi2_sharedTn_vs_sep"],0),2):.3f})')
    # (7) model-free ratio comparison (calibration-free same-J HCN/HCO+ ratios) N vs S etc.
    P('\n# HCN/HCO+ (K km/s) ratio per J and ratio N/S, errors = stat + 7% independent per line (same-band calibration cancels)')
    for a,bb in pairs:
        IA,eA=DATA[a]; IB,eB=DATA[bb]; txt=[]
        for k,J in enumerate(JS):
            rA=IA[k]/IA[3+k]; rB=IB[k]/IB[3+k]; sA=np.hypot(np.hypot(eA[k],eA[3+k]),0.07*np.sqrt(2)); sB=np.hypot(np.hypot(eB[k],eB[3+k]),0.07*np.sqrt(2))
            txt.append(f'J={J}: {rA:.2f}/{rB:.2f}={rA/rB:.2f}+-{rA/rB*np.hypot(sA,sB):.2f} ({np.log(rA/rB)/np.hypot(sA,sB):.1f}sig)')
        P(f'  {a:8s}/{bb:15s} '+'  '.join(txt))

    # (6b) consistency with the 0.25" peak T_b (Q4 s43/s48): can the 10" single zone reach T_b(HCO+ 3-2) of the compact peaks (eta<=1 at 0.25")?
    P('\n# Posterior of T_R(3-2) of the 10" zone vs observed 0.25" peak T_b (N: HCN 1.46, HCO+ 1.56 K; S: HCN 0.79, HCO+ 1.87 K); cal20')
    for rg,(tH,tO) in (('N',(1.461,1.556)),('S',(0.792,1.869))):
        Pp=store[('cal20',rg)][0]; idx=Pp['idx']; p=Pp['p']
        TRO=np.exp(mod.lO[mod.iT[idx],mod.iN[idx],mod.iO[idx],1])/(1.0645*150); TRH=np.exp(mod.lH[mod.iT[idx],mod.iN[idx],mod.iH[idx],1])/(1.0645*150)
        def wq(x):
            o=np.argsort(x); cs=np.cumsum(p[o]); return [round(float(x[o][np.searchsorted(cs,qq/100)]),3) for qq in (16,50,84)]
        d=dict(TR_HCO=wq(TRO),TR_HCN=wq(TRH),P_TRO_ge_peak=float(p[TRO>=tO].sum()),P_TRH_ge_peak=float(p[TRH>=tH].sum()),P_both=float(p[(TRO>=tO)&(TRH>=tH)].sum()))
        res[f'peakTb|{rg}']=d; P(f'  {rg}: T_R(HCO+ 3-2) {d["TR_HCO"]} K, T_R(HCN 3-2) {d["TR_HCN"]} K; P(T_R(HCO+)>={tO}) = {d["P_TRO_ge_peak"]:.2f}, P(T_R(HCN)>={tH}) = {d["P_TRH_ge_peak"]:.2f}, both {d["P_both"]:.2f}')
    # figure
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig,ax=plt.subplots(2,3,figsize=(13,7.5),sharex=True,sharey=True)
    for a,(scn,rg) in zip(ax.ravel(),[('cal20','N'),('cal20','S'),('cal20','Bridge'),('cal20','disk(E1,E2,W)'),('ratios-only','N'),('cal20','zoneN')]):
        m=nT_map(mod,store[(scn,rg)][0]); a.pcolormesh(mod.LN,mod.TK,m/m.max(),cmap='viridis',shading='nearest'); a.set_yscale('log'); a.set_title(f'{rg}  [{scn}]',fontsize=10)
    for a in ax[1]: a.set_xlabel('log10 n(H2) [cm-3]')
    for a in ax[:,0]: a.set_ylabel('Tkin [K]')
    fig.suptitle('Marginal posterior of (n, Tkin): single-zone RADEX LVG fit of the HCN+HCO+ J=2-1,3-2,4-3 ladder at the 10" beam (peak-normalised)'); fig.tight_layout(); fig.savefig(OUT+'s53_nT_posteriors.png',dpi=110)
    json.dump(res,open(OUT+'s53_results.json','w'),indent=1,default=float)
    open(OUT+'s53_ladder_fit.txt','w').write('\n'.join(log))
