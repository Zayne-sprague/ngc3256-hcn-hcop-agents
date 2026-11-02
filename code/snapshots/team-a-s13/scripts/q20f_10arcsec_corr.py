# Q20f: team-b p0028 10" ladders refit with a justified, CORRELATED error model (reviewer fix 3).
# data: 6 ln intensities per region (HCN, HCO+ x 2-1, 3-2, 4-3; team-b p0028 10" values, aperture not independent of ours).
# error model: independent per-line term = stat (3%) (+) aperture/beam-matching [2-1 10% (B5 beam 9.4x6.2" ~ the 10" aperture, loses extended flux), 3-2 3%, 4-3 5%];
#   common per-band flux-scale factor g_b (B5 for 2-1, B6 for 3-2, B7 for 4-3), Gaussian sigma 10% in ln, SHARED by HCN and HCO+ (nuisance, marginalised analytically: covariance += sigma_g^2 a_b a_b^T).
# model: ln I = ln TR(T,n,N_m; RADEX LVG grids out/q4c) + ln eta_m + ln(1.0645 w); w uniform 100-300 km/s; eta<=1; shared eta or eta_HCN/eta_HCO+ within x10; flat priors as Q20c.
# usage: q20f_10arcsec_corr.py [sg=0.10] [ap_scale=1] [indep_only=0]
import numpy as np, sys
SG=float(sys.argv[1]) if len(sys.argv)>1 else 0.10; APS=float(sys.argv[2]) if len(sys.argv)>2 else 1.0; INDEP=int(sys.argv[3]) if len(sys.argv)>3 else 0
O='/workspace/out/q4c/'
I={'N':{'hcn':[12.38,4.81,1.62],'hco+':[20.27,10.04,4.51]},'S':{'hcn':[11.66,4.26,1.35],'hco+':[20.89,9.95,4.28]}}
ap=np.array([0.10,0.03,0.05])*APS; st=0.03
s_ind=np.hypot(st,ap); s_ind=np.tile(s_ind,2)   # order: hcn(2-1,3-2,4-3), hco(2-1,3-2,4-3)
if INDEP:   # old-style uncorrelated 15%
    C=np.diag(np.full(6,0.15**2))
else:
    C=np.diag(s_ind**2)
    for b in range(3):
        a=np.zeros(6); a[b]=1; a[3+b]=1; C+=SG**2*np.outer(a,a)
Ci=np.linalg.inv(C)
LE=np.linspace(np.log(1e-3),0,41); W=np.linspace(100,300,5)
for suf in ('','_tbg5'):
    G={m:np.load(O+f'grid4_{m}{suf}.npz') for m in ('hcn','hco+')}; TK=G['hcn']['TK']
    TRh=G['hcn']['TR'][...,1:].astype(float); TRo=G['hco+']['TR'][...,1:].astype(float); TAU=G['hcn']['tau'][...,2]
    shp=(len(TK),21,23,23)
    Th=np.broadcast_to(TRh[:,:,:,None,:],shp+(3,)); To=np.broadcast_to(TRo[:,:,None,:,:],shp+(3,))
    ok=(Th>0).all(-1)&(To>0).all(-1)&np.isfinite(Th).all(-1)&np.isfinite(To).all(-1)
    lm=np.log(np.where(ok[...,None],np.concatenate([Th,To],-1),1.))
    wT=np.log(TK[1:]/TK[:-1]); wT=np.concatenate([[wT[0]],wT]); 
    for reg in ('N','S'):
        d=np.log(np.array(I[reg]['hcn']+I[reg]['hco+']))
        r=(d-lm).reshape(-1,6); A=r@Ci; rCr=(A*r).sum(1); ok1=ok.reshape(-1)
        Tidx=np.broadcast_to(np.arange(len(TK))[:,None,None,None],shp).reshape(-1)
        for case in ('shared','sep'):
            etas=[(e,0.) for e in LE] if case=='shared' else [(e,q) for e in LE for q in np.linspace(-np.log(10),np.log(10),9)]
            etas=[(e,q) for e,q in etas if e+q<=0]
            res=[]; info=[]
            for e,q in etas:
                eh,eo=e+q,e
                terms=[]
                for w in W:
                    o=np.array([eh]*3+[eo]*3)+np.log(1.0645*w)
                    c2=rCr-2*(A@o)+o@Ci@o
                    terms.append(-0.5*c2)
                L=np.logaddexp.reduce(np.stack(terms),0)-np.log(len(W))
                res.append(np.where(ok1,L,-np.inf)+np.log(wT[Tidx]))
            R=np.stack(res); mx=R.max(); P=np.exp(R-mx); Z=P.sum()
            pT=np.array([P[:,Tidx==i].sum() for i in range(len(TK))])/Z
            cdf=np.cumsum(pT); q_=lambda p: np.interp(p,cdf,TK)
            eo=np.exp(np.array([e for e,q in etas])); pe=P.sum(1)/Z; eta_med=np.exp((pe*np.log(eo)).sum())
            # profile chi2 vs T for eta_o>=0.05 and all
            prof={}
            for lab,sel in (('all',eo>0),('eta>=0.05',eo>=0.05)):
                c=-2*(R[sel]-np.log(wT[Tidx])[None,:]).max(0)
                prof[lab]={int(t):c[Tidx==i].min() for i,t in enumerate(TK) if t in (20,25,30,50,80,150)}
            tau=(P.sum(0)*TAU[:,:,:,None].repeat(23,3).reshape(-1)).sum()/Z
            print(f"SG={SG} APS={APS} INDEP={INDEP} {reg} tbg{suf or '2.73'} {case:6s} lnZ={mx+np.log(Z):.2f} T {q_(.5):.0f} [{q_(.16):.0f},{q_(.84):.0f}] P(T>50)={pT[TK>50].sum():.2f} P(T<=30)={pT[TK<=30].sum():.2f} eta_med {eta_med:.3f} <tau32H>={tau:.0f}")
            print('     chi2min all',round(min(prof['all'].values()),1),' chi2 by T (all eta):',{k:round(v,1) for k,v in prof['all'].items()},' (eta>=0.05):',{k:round(v,1) for k,v in prof['eta>=0.05'].items()})
