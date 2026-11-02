"""Q59 part 2: (a) cross-EB flux-scale check (X95d vs X961, 12C and 13C come from different EBs) from line-free continuum at N, S, 2.2";
(b) inference from the s350 budget: double ratio D = [H13X/HX]_N / [H13X/HX]_S, 12C-line opacity from H13X/HX for a 12C/13C prior,
(c) likelihood of two hypotheses from team-a p0059 (column+abundance: D~2.5; excitation/common column: D~0.6-1.0).
usage: python3 -I s351_q59_inference.py  (needs out/q59/s350.json) -> out/q59/s351.json, s351.txt"""
import sys, json, os, numpy as np
sys.path.insert(0,'/workspace/scripts')
from s350_q59_isotopologue_budget import *
from scipy.optimize import brentq
log=[]
def P(s=''): print(s,flush=True); log.append(s)
# ---------- (a) cross-EB flux scale ----------
# The archive cubes are continuum-subtracted (line-free Tb ~ 0 +- 0.5 mK), so a continuum bridge between X95d (13C) and X961 (12C) is NOT possible.
# X95d and X961 are consecutive EBs of one night (2016-03-07, 04:26 and 05:30 UT) with the same project calibration; X961 matches X145 (2015.1.00412) and
# Harada+18 Table 13 to <=5% (Q29 s183). A 5% scale term per ratio is carried (s350 ratio_err); it is negligible next to the 15-45% statistical errors.
P('== (a) continuum bridge between EBs impossible (cubes are continuum-subtracted); 5% cross-EB scale carried ==')
# ---------- (b) inference ----------
R=json.load(open(O+'s350.json'))
def stats(line,rn):
    s=[r for r in R if r['line']==line and r['reg']==rn]; x=np.array([r['ratio'] for r in s]); e=np.median([r['ratio_err'] for r in s])
    sd=np.sqrt(e**2+x.std()**2); return float(x.mean()),float(sd),float(e),float(x.std())
tab={}
P('\n== (b) ratios (mean over 96 variants per region/line; total error = sqrt(median stat^2 + variant spread^2)) ==')
for line in ('H13CN','H13COp'):
    for rn in ('N','S'):
        m,sd,e,s=stats(line,rn); tab[f'{line}_{rn}']=dict(mean=m,err=sd,stat=e,spread=s)
        P(f' {line:7s} {rn}: {m:.4f} +- {sd:.4f} (stat {e:.4f}, spread {s:.4f});  12C/13C-equivalent 1/ratio = {1/m if m>0 else np.inf:.0f}; 3sig upper = {m+3*sd:.4f}')
def D_post(line):
    rng=np.random.default_rng(1); n=400000
    mN,eN=tab[f'{line}_N']['mean'],tab[f'{line}_N']['err']; mS,eS=tab[f'{line}_S']['mean'],tab[f'{line}_S']['err']
    rN=rng.normal(mN,eN,n); rN=rN[rN>0]
    # posterior on D with flat prior in log D over 0.2..20 ; predicted rS = rN/D ; Gaussian likelihood of the S measurement
    lD=np.linspace(np.log(0.2),np.log(20),400); D=np.exp(lD)
    L=np.array([np.mean(np.exp(-0.5*((mS-rN[:20000]/d)/eS)**2)) for d in D]); post=L/np.trapezoid(L,lD)
    cdf=np.cumsum(post)*(lD[1]-lD[0]); q=lambda p: float(np.exp(np.interp(p,cdf/cdf[-1],lD)))
    def Lh(lo,hi): m=(D>=lo)&(D<=hi); return float(np.trapezoid(L[m],lD[m])/(lD[m][-1]-lD[m][0]))
    return dict(D_median=q(.5),D_16=q(.16),D_84=q(.84),D_05=q(.05),P_D_gt_2=float(np.trapezoid(post[D>=2],lD[D>=2])),P_D_lt_1p2=float(np.trapezoid(post[D<=1.2],lD[D<=1.2])),
                marg_L_excitation_0p6_1p0=Lh(0.6,1.0),marg_L_column_2_3=Lh(2.0,3.0),BF_col_over_exc=Lh(2.0,3.0)/Lh(0.6,1.0),
                chi2_D1=float(((mS-mN/1.0)/np.hypot(eS,eN))**2),chi2_D2p5=float(((mS-mN/2.5)/np.hypot(eS,eN/2.5))**2),ratio_N_over_S_direct=mN/mS if mS>0 else None)
post={}
P('\n== double ratio D = (H13X/HX)_N/(H13X/HX)_S: posterior from the S measurement given the N value (flat in log D 0.2-20) ==')
P(' hypotheses: exc/common column (D 0.6-1.0) vs column+abundance (D 2-3) as in p0059')
for line in ('H13CN','H13COp'):
    post[line]=D_post(line); P(f' {line}: '+', '.join(f'{k}={v:.3g}' for k,v in post[line].items()))
# combined (independent datasets of the same cubes: treat as two independent probes of the same D; HCO+ and HCN could differ in abundance, so only joint as an exploratory product)
# opacity from 13/12 ratio: r = (1/X) * tau/(1-exp(-tau)) * (Tex13 = Tex12 assumed; tau13=tau/X, thin 13C)
def tau_of(r,X):
    t=r*X
    if t<=1: return 0.0
    return brentq(lambda tt: tt/(1-np.exp(-tt))-t,1e-6,60)
P('\n== (b2) 12C opacity from the 13/12 line ratio under equal Tex and thin 13C: r*X = tau/(1-exp(-tau)) (X = 12C/13C) ==')
op={}
for line,nm in (('H13CN','HCN'),('H13COp','HCO+')):
    for rn in ('N','S'):
        m,sd=tab[f'{line}_{rn}']['mean'],tab[f'{line}_{rn}']['err']; row={}
        for X in (40,60,80,100):
            t=[tau_of(max(r,0),X) for r in (m-sd,m,m+sd)]; row[X]=t
        op[f'{nm}_{rn}']=row
        P(f' {nm:4s} {rn}: '+' | '.join(f'X={X}: tau {v[1]:.2f} [{v[0]:.2f},{v[2]:.2f}]' for X,v in row.items()))
# Monte Carlo: tau posterior under log-flat X prior 40-100 (Tang+19 range for starbursts ~40-100; Martin+19 NGC 253 ~ 40-?)
rng=np.random.default_rng(2); mc={}
P('\n== tau(12C 1-0) posterior, X ~ logflat(40,100), ratio Gaussian truncated to the physical region r>=0 ==')
for line,nm in (('H13CN','HCN'),('H13COp','HCO+')):
    for rn in ('N','S'):
        m,sd=tab[f'{line}_{rn}']['mean'],tab[f'{line}_{rn}']['err']; r=rng.normal(m,sd,20000); X=np.exp(rng.uniform(np.log(40),np.log(100),20000))
        t=np.array([tau_of(max(a,0),b) for a,b in zip(r,X)]); mc[f'{nm}_{rn}']=[float(np.percentile(t,p)) for p in (5,16,50,84,95)]
        P(f' {nm:4s} {rn}: tau percentiles 5/16/50/84/95 = '+' '.join('%.2f'%v for v in mc[f'{nm}_{rn}'])+f'  P(tau>1)={np.mean(t>1):.2f}  P(tau>3)={np.mean(t>3):.2f}')
# difference in tau N vs S
for nm in ('HCN','HCO+'):
    a=np.array(mc[f'{nm}_N']); b=np.array(mc[f'{nm}_S']); P(f' {nm}: median tau N {a[2]:.2f}, S {b[2]:.2f}')
# Harada Table 13 comparison (1.7", peak pixel, K km/s): H13CO+ 9.4+-0.95, HCO+ 205.25 -> ratio
P('\n== comparison with Harada+2018 Table 13 (N, 1.7", peak pixel): H13CO+/HCO+ = 9.4/205.25 = %.3f +- %.3f; H13CN not listed =='%(9.4/205.25,9.4/205.25*0.1))

# ---------- (c) 13C-based N(HCN)/N(HCO+) contrast (optically thin 13C, same 12C/13C for both species) vs the 12C contrast ----------
rng=np.random.default_rng(3); n=300000
g=lambda k: rng.normal(tab[k]['mean'],tab[k]['err'],n)
rNc,rNo,rSo,rSc=g('H13CN_N'),g('H13COp_N'),g('H13COp_S'),g('H13CN_S')
ok=(rNc>0)&(rNo>0)&(rSo>0)
Cd=np.exp(np.linspace(np.log(0.3),np.log(30),300)); lC=np.log(Cd)
pred=(rNc[ok]/rNo[ok])[:,None]*rSo[ok][:,None]/Cd[None,:]     # predicted S H13CN/HCN
sg=tab['H13CN_S']['err']; L=np.mean(np.exp(-0.5*((tab['H13CN_S']['mean']-pred)/sg)**2),axis=0); pst=L/np.trapezoid(L,lC)
cdf=np.cumsum(pst)*(lC[1]-lC[0]); cdf/=cdf[-1]; q=lambda p: float(np.exp(np.interp(p,cdf,lC)))
C12=1.45
res13=dict(C13_median=q(.5),C13_16=q(.16),C13_84=q(.84),C13_05=q(.05),C13_lower95=q(.05),P_C13_gt_C12=float(np.trapezoid(pst[Cd>C12],lC[Cd>C12])),P_C13_gt_2p5=float(np.trapezoid(pst[Cd>2.5],lC[Cd>2.5])),
           P_C13_lt_1p0=float(np.trapezoid(pst[Cd<1.0],lC[Cd<1.0])),N_ratio_13=float(np.median(rNc[ok]/rNo[ok])))
P('\n== (c) 13C-based (thin) N(HCN)/N(HCO+) N/S contrast C13 from [H13CN/HCN]_S measured vs predicted [rNc/rNo]/C * rSo (flat in log C 0.3-30) ==')
P(' 12C intensity contrast (HCN/HCO+ N/S at 2.2-2.5", Q56/p0059) = 1.45; excitation-only/thin: C13 ~ C12; opacity compression at N: C13 > C12')
P(' '+', '.join(f'{k}={v:.3g}' for k,v in res13.items()))
P(' N: H13CN/H13CO+ = %.2f; HCN/HCO+ (12C) = %.2f => (12C ratio / 13C ratio)_N = %.2f  [=1 if both thin & same X]'%(rNc[ok].mean()/rNo[ok].mean(),0.83,0.83/(rNc[ok].mean()/rNo[ok].mean())))
json_extra=dict(contrast13=res13)
json.dump(dict(table=tab,double_ratio=post,tau_fixedX=op,tau_mc=mc,**json_extra),open(O+'s351.json','w'),indent=1,default=float)
open(O+'s351.txt','w').write('\n'.join(log))
