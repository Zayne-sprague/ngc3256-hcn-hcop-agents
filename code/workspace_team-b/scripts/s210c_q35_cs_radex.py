# Q35 fix 1+2: CS excitation at N, S, C(Bridge) and column-density / abundance ratios N(HCN)/N(HCO+), N(CS)/N(HCN), N(CS)/N(HCO+).
# Model: the team RADEX LVG grid (s62_lib.Grid: Tk 10-300 K, log n 3-8, log N(HCO+), log A=N(HCN)/N(HCO+), dv=150 km/s, HCN/HCO+ J=1..4) + a CS grid (s210a, same conventions; log N(CS) 11.5-17.5 flat)
# with ONE common Tk, n and beam filling eta<=1 (flat in ln eta, marginalised analytically) for all species (same gas).
# Data at ONE beam (Harada+2018 Table 13, 1.7", peak pixel, K km/s; lit/harada_table13_parsed.csv): HCN 1-0, 3-2; HCO+ 1-0, 3-2; CS 2-1, 5-4. Errors: table rms + 5% per line; shared band scales (1-0 band 3: 10%, 3-2/5-4 band 6: 15%) marginalised with 7 Sobol draws as s186.
# Variants: 'hh' = HCN+HCO+ lines only (CS column free; then CS 5-4/2-1, CS 7-6/HCN 4-3 are PREDICTIONS compared with observation = excitation tension test);
#           'ladder' = all six 1.7" lines; 'ladder+4.3' = plus the 4.3" K ratios CS7-6/HCN4-3, HCN4-3/HCO+4-3 (s200 deblend; errors from the s210b injection; eta cancels); 'noCS5' = CS 2-1 only.
# usage: python3 -I s210c_q35_cs_radex.py -> out/q35/s210c.txt/json
import sys, json, numpy as np, pandas as pd
sys.path.insert(0,'/workspace/scripts')
import s62_lib as S
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.special import log_ndtr
from scipy.stats import qmc, norm
O='/workspace/out/q35/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
F=S.Grid(0.2); Tk=10**F.lT; valid=(Tk<=300.01)
d=np.load(O+'grid_cs.npz'); TKg=np.log10(d['TK']); LN0=d['LN']; LC0=d['LNC']; cs=d['cs']
NC=np.round(np.arange(11.5,17.51,0.2),2); nC=len(NC)
def ci(k):
    a=cs[:,:,:,k,:]; ok=(a[...,0]>0)&np.isfinite(a[...,0])
    ls={'I':np.where(ok,np.log(1.0645*150*np.where(ok,a[...,0],1)),np.nan),'T':np.where(ok,a[...,1],np.nan),'tau':np.where(ok,np.log10(np.where(ok,a[...,2],1)),np.nan)}
    out={q:np.empty((nC,F.nm),np.float32) for q in ls}
    for ic,nc in enumerate(NC):
        pts=np.stack([F.lT,F.ln,np.full(F.nm,nc)],1)
        for q in ls: out[q][ic]=RGI((TKg,LN0,LC0),ls[q],bounds_error=False)(pts)
    return out
C2,C5,C7=ci(1),ci(4),ci(6)
H=pd.read_csv('/workspace/lit/harada_table13_parsed.csv'); H['sp']=H.species_transition.str.replace(' ','')
def g(p,s): r=H[(H.pos==p)&(H.sp.str.startswith(s))].iloc[0]; return r.I_Kkms,r.eI_Kkms
LINES=['HCN1','HCN3','HCO+1','HCO+3','CS2','CS5']
R43={'N':((0.157,0.21),(0.552,0.14)),'S':((0.179,0.38),(0.331,0.24)),'C':((0.222,0.60),(0.371,0.33))}   # (CS7-6/HCN4-3), (HCN4-3/HCO+4-3), K ratios, fractional errors (s210b box)
base=np.stack([F.M1[:,0],F.M[:,1],F.M1[:,1],F.M[:,4]],1).astype(np.float64)   # HCN1,HCN3,HCO1,HCO3 ln I
lHCN4=F.M[:,2]; lHCO4=F.M[:,5]
def wq(x,w,p=(16,50,84)):
    m=np.isfinite(x)&(w>0); x=x[m]; w=w[m]; o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return [float(x[o][min(np.searchsorted(c,a/100),len(x)-1)]) for a in p]
def fit(pos,variant,nd=7,sb=(0.10,0.15),sl=0.05,seed=3):
    yv=np.array([np.log(g(pos,s)[0]) for s in LINES]); rel=np.array([g(pos,s)[1]/g(pos,s)[0] for s in LINES])
    use=np.ones(6,bool)
    if variant=='hh': use[4:]=False
    if variant=='noCS5': use[5]=False
    u=qmc.Sobol(d=2,scramble=True,seed=seed).random(8)[:nd]; z=norm.ppf(np.clip(u,1e-9,1-1e-9))*np.array(sb)
    W=np.zeros((nC,F.nm)); var=(rel**2+sl**2)[use]; Ci=1/var; a=Ci.sum()
    for zz in z:
        yy=(yv-np.array([zz[0],zz[1],zz[0],zz[1],zz[0],zz[1]]))[use]
        for ic in range(nC):
            Mm=np.concatenate([base,C2['I'][ic][:,None],C5['I'][ic][:,None]],1)[:,use]
            ok=valid&np.isfinite(Mm).all(1)&np.isfinite(C7['I'][ic])
            r=yy[None,:]-np.where(ok[:,None],Mm,0); b=(r*Ci).sum(1); chi=(r*r*Ci).sum(1)-b*b/a; l=-0.5*chi+log_ndtr(np.sqrt(a)*(0-b/a))
            if variant=='ladder+4.3':
                (R1,f1),(R2,f2)=R43[pos]; m1=np.exp(C7['I'][ic]-lHCN4); m2=np.exp(lHCN4-lHCO4)
                l=l-0.5*(np.log(R1/m1)/f1)**2-0.5*(np.log(R2/m2)/f2)**2
            W[ic]+=np.where(ok,np.exp(np.where(ok,l,-1e30)),0)
    # hh variant: N(CS) is unconstrained by the data -> flat prior on log N(CS)/N(HCN) in [-1.5, 0] (CS/HCN abundance ratio 0.03-1; Harada CS/HCN line ratios 0.2-0.3)
    if variant=='hh': W*=((NC[:,None]-F.lNH[None,:]>=-1.5-1e-6)&(NC[:,None]-F.lNH[None,:]<=1e-6))
    # best-fit chi2 (no band offsets) for the diagnostics: n_data - n_param(T,n,N_O,A,N_CS,eta)
    ch=1e30
    for ic in range(nC):
        Mm=np.concatenate([base,C2['I'][ic][:,None],C5['I'][ic][:,None]],1)[:,use]; ok=valid&np.isfinite(Mm).all(1)&np.isfinite(C7['I'][ic])
        r=(yv[use])[None,:]-np.where(ok[:,None],Mm,0); b=(r*Ci).sum(1); chi=(r*r*Ci).sum(1)-b*b/a
        eta_ok=(b/a<=0)   # eta<=1 : mean offset ln eta <=0
        chi=np.where(eta_ok,chi,(r*r*Ci).sum(1)); chi=np.where(ok,chi,1e30)
        if variant=='ladder+4.3':
            (R1,f1),(R2,f2)=R43[pos]; chi=chi+(np.log(R1/np.exp(C7['I'][ic]-lHCN4))/f1)**2+(np.log(R2/np.exp(lHCN4-lHCO4))/f2)**2
        ch=min(ch,float(np.nanmin(chi)))
    global CHI; CHI=ch
    return W/W.sum()
res={}
for pos,lab in (('N','N'),('S','S'),('C','Bridge(C)')):
    P(f'\n######## {lab}: Harada 1.7" K km/s: '+', '.join(f'{s} {g(pos,s)[0]:.1f}' for s in LINES)+f' | CS5-4/2-1 = {g(pos,"CS5")[0]/g(pos,"CS2")[0]:.3f}')
    for var in ('hh','ladder','ladder+4.3','noCS5'):
        W=fit(pos,var); w=W.sum(0); wf=W.reshape(-1); chi_=CHI
        NCg=np.repeat(NC[:,None],F.nm,1)
        lNH=np.broadcast_to(F.lNH,(nC,F.nm)); lNO=np.broadcast_to(F.lNO,(nC,F.nm)); lA=np.broadcast_to(F.lA,(nC,F.nm))
        f=lambda x: wq(np.asarray(x,float).reshape(-1),wf)
        out=dict(Tk=wq(Tk,w),logn=wq(F.ln,w),logA=wq(F.lA,w),logNO=wq(F.lNO,w),PTk50=float(w[Tk>50].sum()),
                 logNCS_NHCN=f(NCg-lNH),logNCS_NHCO=f(NCg-lNO),logNCS=f(NCg),
                 TexCS21=f(C2['T']),TexCS54=f(C5['T']),TexCS76=f(C7['T']),tauCS54=f(C5['tau']),
                 r54_21=f(np.exp(C5['I']-C2['I'])),CS76_HCN43=f(np.exp(C7['I']-lHCN4[None,:])),HCN43_HCO43=wq(np.exp(lHCN4-lHCO4),w),
                 TexCS21_HCN=None)
        out.pop('TexCS21_HCN'); res[f'{pos}|{var}']=out; out['chi2min']=chi_
        P(f'  [{var:10s}] Tk {out["Tk"][1]:.0f} [{out["Tk"][0]:.0f},{out["Tk"][2]:.0f}] P(>50)={out["PTk50"]:.2f}; log n {out["logn"][1]:.2f} [{out["logn"][0]:.2f},{out["logn"][2]:.2f}]; logA {out["logA"][1]:.2f} [{out["logA"][0]:.2f},{out["logA"][2]:.2f}]'
          f' | log N(CS)/N(HCN) {out["logNCS_NHCN"][1]:.2f} [{out["logNCS_NHCN"][0]:.2f},{out["logNCS_NHCN"][2]:.2f}]; CS/HCO+ {out["logNCS_NHCO"][1]:.2f} [{out["logNCS_NHCO"][0]:.2f},{out["logNCS_NHCO"][2]:.2f}]'
          f' | Tex CS 2-1/5-4/7-6 {out["TexCS21"][1]:.0f}/{out["TexCS54"][1]:.0f}/{out["TexCS76"][1]:.0f} K; pred CS5-4/2-1 {out["r54_21"][1]:.2f} [{out["r54_21"][0]:.2f},{out["r54_21"][2]:.2f}]; pred CS7-6/HCN4-3 {out["CS76_HCN43"][1]:.3f} [{out["CS76_HCN43"][0]:.3f},{out["CS76_HCN43"][2]:.3f}]; HCN/HCO+ 4-3 pred {out["HCN43_HCO43"][1]:.2f}; best-fit chi2 {chi_:.2f}')
        res[f'{pos}|{var}|W']=None
# N/S contrasts of column ratios (independent posteriors; samples drawn from the weights)
rng=np.random.default_rng(1)
def samp(pos,var,n=4000):
    W=fit(pos,var); i=rng.choice(W.size,n,p=W.reshape(-1)); ic,im=np.unravel_index(i,W.shape)
    J=lambda: rng.uniform(-0.1,0.1,n)   # +-0.1 dex: the A and N grids have 0.2 dex cells (discretisation of the posterior)
    return dict(NCS_NHCN=NC[ic]-F.lNH[im]+J()+J(),NCS_NHCO=NC[ic]-F.lNO[im]+J()+J(),A=F.lA[im]+J()+J()*0)
P('\n######## N/S contrasts of column ratios (posterior medians and 16-84%, ratio of independent posteriors)')
for var in ('hh','ladder','ladder+4.3'):
    sN=samp('N',var); sS=samp('S',var); sC=samp('C',var)
    for nm,key in (('N(HCN)/N(HCO+)','A'),('N(CS)/N(HCN)','NCS_NHCN'),('N(CS)/N(HCO+)','NCS_NHCO')):
        if var=='hh' and key!='A': continue
        for a,b,tag in ((sN,sS,'N/S'),(sC,sS,'C/S')):
            x=10**(a[key]-b[key]); q=np.percentile(x,[16,50,84]); res[f'ratio|{var}|{nm}|{tag}']=q.tolist()
            P(f'  [{var:10s}] {nm} {tag} = {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}]; P(>1) = {(x>1).mean():.2f}')
    if var!='hh':
        for nm,key in (('N(HCN)/N(CS)','NCS_NHCN'),('N(HCO+)/N(CS)','NCS_NHCO')):
            for a,b,tag in ((sN,sS,'N/S'),):
                x=10**(-(a[key]-b[key])); q=np.percentile(x,[16,50,84]); res[f'ratio|{var}|{nm}|{tag}']=q.tolist(); P(f'  [{var:10s}] {nm} {tag} = {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}]; P(>1) = {(x>1).mean():.2f}')
res={k:v for k,v in res.items() if v is not None}
open(O+'s210c.txt','w').write('\n'.join(log)); json.dump(res,open(O+'s210c.json','w'),indent=1,default=float)
