# Q3 round 2, fixes 4,5: (a) envelope + dense compact core fit at N and S; (b) tests of the A (=N(HCN)/N(HCO+)) result: separate filling for HCN and HCO+, free line width, fine-grid intervals.
# Core: Xb4b r<0.5" flux (3-2) scaled by g in [0,3] (g=1 measured; g>1 allows more compact gas than the 0.25" annuli recover, g<1 less) and a RADEX core ladder L_J=T(J)/T(3-2) drawn from LVG models with
#   n 10^5.5-10^7, Tk 50-150 K, N(HCO+) 14-18, A in -1..1.5 whose T_R(3-2) >= observed 0.25" peak (1.46 HCN / 1.56 HCO+ K, Q4 s43) and whose HCN/HCO+ T_R(3-2) ratio is 0.94+-0.15 (observed core ratio, equal filling).
#   Core flux error 15% (added to the envelope error). Envelope = observed 10" ladder minus core, fitted with the single-zone model (nominal calibration, no recovery, plus the one-sided variant).
# usage: python3 -I s55d_core_A_tests.py -> out/q3b/s55d_core_A_tests.txt/json
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s55_lib as S, s55b_recovery_fit as B
from scipy.special import logsumexp
import pandas as pd
L=S.L; F=B.F; log=[]; res={}
def P(s=''): print(s,flush=True); log.append(s)
# ---- core fluxes (Jy km/s) from s43
b=pd.read_csv('/workspace/out/s43_profiles.csv'); bm=1.1331*0.2456*0.2373/0.047**2
CORE={z:{m:float((b[b.region==z].reset_index(drop=True)[c]*b[b.region==z].reset_index(drop=True).npix/bm)[:2].sum()) for m,c in (('HCN','I_HCN'),('HCO+','I_HCO'))} for z in ('N','S')}
PEAK={'N':(1.461,1.556),'S':(0.792,1.869)}     # 0.25" peak Tb (HCN,HCO+) K, s43
# ---- core ladder population
TR=lambda col: np.exp(F.M[:,col])/(1.0645*150)
tr=[TR(c) for c in range(6)]
def core_pop(z):
    pk=PEAK[z]; ratio=tr[1]/tr[4]; target=(CORE[z]['HCN']/CORE[z]['HCO+'])*(L.NU('HCO+',3)/L.NU('HCN',3))**2
    m=(F.ln>=5.5)&(F.ln<=7)&(10**F.lT>=50)&(10**F.lT<=150)&(F.lNO>=14)&(F.lNO<=18)&(tr[1]>=pk[0])&(tr[4]>=pk[1])&(np.abs(np.log(ratio/target))<0.15)
    return m,target
LAD={}
for z in ('N','S'):
    m,target=core_pop(z); idx=np.where(m)[0]
    lad=np.stack([tr[0][idx]/tr[1][idx],tr[2][idx]/tr[1][idx],tr[3][idx]/tr[4][idx],tr[5][idx]/tr[4][idx]],1)   # L2H,L4H,L2O,L4O
    P(f'core ladder population {z}: {len(idx)} models (n 5.5-7, Tk 50-150, T_R(3-2)>=0.25" peak, HCN/HCO+ ratio {target:.2f}+-15%); L=T(J)/T(3-2) 16/50/84%: '+
      ', '.join(f'{nm} {np.percentile(lad[:,i],16):.2f}/{np.percentile(lad[:,i],50):.2f}/{np.percentile(lad[:,i],84):.2f}' for i,nm in enumerate(('L2 HCN','L4 HCN','L2 HCO+','L4 HCO+')))+
      f'; fraction with tau(HCN3-2)>1: {np.mean(10**F.tauH[idx,1]>1):.2f}')
    o=np.argsort(lad[:,1]); sel=o[np.linspace(0,len(o)-1,9).astype(int)]; LAD[z]=lad[sel]
# ---- envelope + core fit
def fit_env(z,sig_b=(0.15,0.10,0.15),gs=(0,0.5,1,1.5,2,2.5,3),rec=None):
    I,e=L.DATA[z]; c3={'HCN':L.tk(CORE[z]['HCN'],'HCN',3,L.OM10),'HCO+':L.tk(CORE[z]['HCO+'],'HCO+',3,L.OM10)}
    tot=None; lnZg={g:[] for g in gs}; fr={g:[] for g in gs}; ok_all=[]
    for g in gs:
        for lad in LAD[z]:
            c=np.array([c3['HCN']*lad[0],c3['HCN'],c3['HCN']*lad[1],c3['HCO+']*lad[2],c3['HCO+'],c3['HCO+']*lad[3]])*g
            env=I-c
            if (env<=0.05*I).any(): lnZg[g].append(-np.inf); continue
            rel=np.hypot(I*e,c*0.15)/env; En=S.Engine(F,S.covmat(rel,sig_b))
            lnL,mu,chi=En.lnlike(np.log(env) if rec is None else np.log(env)-np.log(rec)); lnZg[g].append(logsumexp(lnL)-np.log(len(lnL)))
            tot=lnL if tot is None else np.logaddexp(tot,lnL); fr[g].append(c/I)
    return tot,lnZg,fr
for z in ('N','S'):
    P(f'\n## {z}: envelope + compact core (core 3-2 flux Xb4b r<0.5": HCN {CORE[z]["HCN"]:.2f}, HCO+ {CORE[z]["HCO+"]:.2f} Jy km/s)')
    tot,lnZg,fr=fit_env(z); Pst=S.Post(F,tot-np.log(1)); d=Pst.summary(); res[f'core+env|{z}']=d
    base=S.Post(F,S.Engine(F,S.covmat(L.DATA[z][1])).lnlike(np.log(L.DATA[z][0]))[0]).summary()
    P('  single zone, no core:        '+B.line(base)); P('  envelope with core removed:  '+B.line(d)); P('  '+B.pline(d))
    zs=np.array([logsumexp(lnZg[g]) for g in lnZg]); gsl=list(lnZg); pz=np.exp(zs-zs.max()); pz/=pz.sum()
    P('  posterior on core scale g (g=1: measured Xb4b core): '+', '.join(f'g={g}: {p:.2f}' for g,p in zip(gsl,pz)))
    cum=np.cumsum(pz); g95=gsl[int(np.searchsorted(cum,0.95))]; res[f'g95|{z}']=g95
    # core flux fractions of the 10" ladder for g=1 and g95 (median over ladders)
    for g in (1,g95):
        f=np.array(fr[g]); 
        if len(f): P(f'  core fraction of the 10" ladder at g={g}: HCN 2-1/3-2/4-3 = '+'/'.join(f'{np.median(f[:,i]):.2f}' for i in range(3))+' ; HCO+ = '+'/'.join(f'{np.median(f[:,3+i]):.2f}' for i in range(3))+f' (range over core ladders {f.min():.2f}-{f.max():.2f})')
    res[f'frac|{z}']={str(g):np.median(np.array(fr[g]),0).tolist() for g in gsl if len(fr[g])}
    # positivity-only limit (the largest g for which any core ladder leaves envelope>5% in all 6 lines)
    P('  g accepted by positivity (envelope >5% of the 10" flux in all 6 lines): '+', '.join(f'g={g}: {sum(np.isfinite(v) for v in lnZg[g])}/9 ladders' for g in gsl))
# ---------------------------------------------------------------- (b) A, eta, width
P('\n## A tests at N (nominal calibration 15/10/15%, no recovery; fine grid: step 0.1 dex in A and Tk 0.05 dex, so intervals are not grid-pinned)')
En={}
for z in ('N','S','Bridge'):
    I,e=L.DATA[z]; En[z]=S.Engine(F,S.covmat(e)); y=np.log(I)
    d=S.Post(F,En[z].lnlike(y)[0]).summary(); res[f'nominal|{z}']=d; P(f'[{z}] nominal: '+B.line(d,('log n','Tk','log A','log N(HCO+)/dv'))+' ; '+B.pline(d))
P('# HCN more compact than HCO+: eta_HCN = q * eta_HCO+ (q<=1; HCN data shifted by -ln q)')
for z in ('N','S'):
    I,e=L.DATA[z]
    for q in (1,0.5,0.25,0.1):
        y=np.log(I).copy(); y[:3]-=np.log(q); d=S.Post(F,En[z].lnlike(y)[0]).summary(); res[f'q{q}|{z}']=d
        P(f'[{z}] q={q}: log A {d["log A"][2]:.2f} [68% {d["log A"][1]:.2f},{d["log A"][3]:.2f}; 95% {d["log A"][0]:.2f},{d["log A"][4]:.2f}], log n {d["log n"][2]:.2f} [{d["log n"][1]:.2f},{d["log n"][3]:.2f}], Tk {d["Tk"][2]:.0f} [{d["Tk"][1]:.0f},{d["Tk"][3]:.0f}]; P(A>1)={d["P"]["A>1"]:.2f}')
    # q marginalised log-uniform 0.1-1 and 0.1-3
    for lab,qs in (('q log-uniform 0.1-1',np.logspace(-1,0,10)),('q log-uniform 0.1-3 (HCN may fill more)',np.logspace(-1,np.log10(3),12))):
        tot=None
        for q in qs:
            y=np.log(I).copy(); y[:3]-=np.log(q); l=En[z].lnlike(y)[0]; tot=l if tot is None else np.logaddexp(tot,l)
        d=S.Post(F,tot).summary(); res[f'{lab}|{z}']=d; P(f'[{z}] {lab}: log A {d["log A"][2]:.2f} [68% {d["log A"][1]:.2f},{d["log A"][3]:.2f}; 95% {d["log A"][0]:.2f},{d["log A"][4]:.2f}], log n {d["log n"][2]:.2f} [{d["log n"][1]:.2f},{d["log n"][3]:.2f}], Tk {d["Tk"][2]:.0f} [{d["Tk"][1]:.0f},{d["Tk"][3]:.0f}]')
P('# Free line width w (the model flux is eta*w*T_R(N/dv); with N/dv free, w only changes the eta<=1 bound: ln eta_max = ln(w/150)); N is quoted per km/s so is unaffected')
for z in ('N','S'):
    I,e=L.DATA[z]; y=np.log(I)
    for w in (60,100,150,200,300):
        Eg=S.Engine(F,S.covmat(e),cap=np.log(w/150.)); d=S.Post(F,Eg.lnlike(y)[0]).summary(); res[f'w{w}|{z}']=d
        P(f'[{z}] w={w} km/s: log n {d["log n"][2]:.2f} [{d["log n"][1]:.2f},{d["log n"][3]:.2f}], Tk {d["Tk"][2]:.0f} [{d["Tk"][1]:.0f},{d["Tk"][3]:.0f}], log A {d["log A"][2]:.2f} [{d["log A"][1]:.2f},{d["log A"][3]:.2f}], log N(HCO+)/dv {d["log N(HCO+)/dv"][2]:.2f}')
json.dump(res,open(S.OUT+'s55d_core_A_tests.json','w'),indent=1,default=float); open(S.OUT+'s55d_core_A_tests.txt','w').write('\n'.join(log))
