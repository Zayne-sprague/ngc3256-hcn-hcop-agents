"""Q29: does adding J=1-0 change the N/S excitation-vs-abundance conclusion?  RADEX (LVG, LAMDA) grid of Q3/Q27 (s62_lib.Grid: Tk 10-500 K, log n 3-8, log N(HCO+), log A=N(HCN)/N(HCO+)),
data = HCN, HCO+ 1-0 (2015.1.00412 / 2015.1.00993, +-350 km/s) and 3-2 (Xb4d) at the SAME circular beam (1.5" and 2.2"), aperture means (K km/s, x(1+z)).
Per position: Gaussian in ln I, 5% per-line + statistical error, band scales 1-0 (sig_b10) and 3-2 (sig_b32) shared by all positions (MC over 41 Sobol draws), common beam filling eta<=1 (flat in ln eta, marginalised analytically), Tk<=300 K, flat priors on grid.
Hypotheses for the pair (P1,P2): free (independent), H_Tn (same Tk and n; columns and A free = abundance/column only), H_A (same A=N(HCN)/N(HCO+); Tk, n, N free = excitation/column only), H_ANO (same A and N(HCO+): excitation only),
H_T, H_n, H_all (same Tk,n,N,A: only filling differs).  Output ln Bayes factors vs free, with 3-2 only and with 1-0 + 3-2.
usage: python3 -I s186_q29_radex_evidence.py -> out/q29/s186.txt/json"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from scipy.special import log_ndtr, logsumexp
from scipy.stats import qmc, norm
import s62_lib as S
O='/workspace/out/q29/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
F=S.Grid(); step=F.step
iT=F.iT; iN=F.iN; iO=np.round((F.lNO-11.5)/step).astype(int); iA=np.round((F.lA+1)/step).astype(int)
nT,nN,nO,nA=iT.max()+1,iN.max()+1,iO.max()+1,iA.max()+1
valid=(10**F.lT<=300.01)&np.isfinite(F.M1).all(1)
Z=0.00927; K=1+Z
cols_all=np.stack([F.M1[:,0],F.M[:,1],F.M1[:,1],F.M[:,4]],1)     # HCN1, HCN3, HCO1, HCO3
r3=json.load(open(O+'s183.json')); r182=json.load(open(O+'s182.json'))
def data(th,reg):
    if th=='1.50': a=r3[f'2015.1.00412 X145|1.5|{reg}']; x=r3['xb4d_32']['1.50'][reg]['win']
    elif th=='2.20': a=r3[f'2015.1.00993 X961|2.2|{reg}']; x=r3['xb4d_32']['2.20'][reg]['win']
    elif th=='2.20x': a=r3[f'2015.1.00412 X145|2.2|{reg}']; x=r3['xb4d_32']['2.20'][reg]['win']
    I=np.array([a['HCN']['w350_ap'][0],x['HCN32'][0],a['HCOp']['w350_ap'][0],x['HCOp32'][0]])*K
    e=np.array([a['HCN']['w350_ap'][1],x['HCN32'][1],a['HCOp']['w350_ap'][1],x['HCOp32'][1]])*K
    return I,e/I
def dense(l):
    D=np.zeros((nT,nN,nO,nA)); D[iT,iN,iO,iA]=np.where(valid,np.exp(l),0.0); return D
def lnl(y,rel,use,sl=0.05):
    var=rel[use]**2+sl**2; r=(y[use]-0)[None,:]-cols_all[:,use]; Ci=1/var; a=Ci.sum()
    b=(r*Ci).sum(1); eh=b/a; chi=(r*r*Ci).sum(1)-b*b/a
    l=-0.5*chi+log_ndtr(np.sqrt(a)*(0-eh)); return np.where(valid&np.isfinite(l),l,-1e30)
def evid(th,pair,use,nd=41,sb=(0.10,0.15),sl=0.05,seed=2):
    u=qmc.Sobol(d=2,scramble=True,seed=seed).random(64)[:nd]; z=norm.ppf(np.clip(u,1e-9,1-1e-9))*np.array(sb)
    d1,e1=data(th,pair[0]); d2,e2=data(th,pair[1])
    acc={k:[] for k in ('free','Tn','A','ANO','T','n','all')}; post={}
    for zz in z:
        off=np.array([zz[0],zz[1],zz[0],zz[1]])   # HCN1,HCN3,HCO1,HCO3 : band scales shared by both lines and both positions
        Ds=[]
        for y,e in ((d1,e1),(d2,e2)):
            l=lnl(np.log(y)-off,e,use,sl); Ds.append(dense(l-l.max()))
        A1,A2=Ds
        m=lambda X,ax: X.mean(axis=ax)
        acc['free'].append(A1.mean()*A2.mean())
        acc['Tn'].append((m(A1,(2,3))*m(A2,(2,3))).mean())
        acc['A'].append((m(A1,(0,1,2))*m(A2,(0,1,2))).mean())
        acc['ANO'].append((m(A1,(0,1))*m(A2,(0,1))).mean())
        acc['T'].append((m(A1,(1,2,3))*m(A2,(1,2,3))).mean())
        acc['n'].append((m(A1,(0,2,3))*m(A2,(0,2,3))).mean())
        acc['all'].append((A1*A2).mean())
    return {k:float(np.log(np.mean(v))-np.log(np.mean(acc['free']))) for k,v in acc.items()}
def post1(th,reg,use,sb=(0.1,0.15),sl=0.05,nd=41,seed=2):
    """single-position posterior (marginal over band draws) of Tk, n, log A, log N(HCO+)"""
    u=qmc.Sobol(d=2,scramble=True,seed=seed).random(64)[:nd]; z=norm.ppf(np.clip(u,1e-9,1-1e-9))*np.array(sb)
    y,e=data(th,reg); w=np.zeros(F.nm)
    for zz in z:
        off=np.array([zz[0],zz[1],zz[0],zz[1]]); l=lnl(np.log(y)-off,e,use,sl); w+=np.where(valid,np.exp(l-l.max()),0)
    w/=w.sum()
    def q(x,p=(16,50,84)):
        o=np.argsort(x); c=np.cumsum(w[o]); return [float(x[o][min(np.searchsorted(c,a/100),len(x)-1)]) for a in p]
    T=10**F.lT
    return dict(Tk=q(T),logn=q(F.ln),logA=q(F.lA),logNO=q(F.lNO),tauHCN32=q(F.tauH[:,2]),tauHCO32=q(F.tauO[:,2]),tauHCO10=q(F.tauO[:,0]),tauHCN10=q(F.tauH[:,0]),PTk50=float(w[T>50].sum()))
res={}
USE={'3-2 only':np.array([False,True,False,True]),'1-0 + 3-2':np.array([True,True,True,True]),'1-0 only':np.array([True,False,True,False])}
for th in ('1.50','2.20'):
    P(f'\n######## beam {th}" (1-0: {"2015.1.00412 X145" if th=="1.50" else "2015.1.00993 X961"}; 3-2: Xb4d). Data K km/s (x1.0093):')
    for g in ('N','S','Bridge'):
        y,e=data(th,g); P(f'  {g}: HCN1 {y[0]:.1f} HCN3 {y[1]:.1f} HCO1 {y[2]:.1f} HCO3 {y[3]:.1f}  rel.stat {np.round(e,3).tolist()}')
    for un,use in USE.items():
        for g in ('N','S','Bridge'):
            p=post1(th,g,use); res[f'post|{th}|{un}|{g}']=p
            P(f'  [{un}] {g}: Tk {p["Tk"][1]:.0f} [{p["Tk"][0]:.0f},{p["Tk"][2]:.0f}] P(Tk>50)={p["PTk50"]:.2f}; log n {p["logn"][1]:.2f} [{p["logn"][0]:.2f},{p["logn"][2]:.2f}]; log A {p["logA"][1]:.2f} [{p["logA"][0]:.2f},{p["logA"][2]:.2f}]; log N(HCO+)/150km/s {p["logNO"][1]:.1f}; tau(HCN 3-2) {10**p["tauHCN32"][1]:.1f} [{10**p["tauHCN32"][0]:.1f},{10**p["tauHCN32"][2]:.1f}]; tau(HCO+ 1-0) {10**p["tauHCO10"][1]:.1f}; tau(HCN 1-0) {10**p["tauHCN10"][1]:.1f}')
    for pair in (('N','S'),('N','Bridge'),('S','Bridge')):
        for un in ('3-2 only','1-0 + 3-2'):
            E=evid(th,pair,USE[un]); res[f'ev|{th}|{pair[0]}-{pair[1]}|{un}']=E
            P(f'  ln BF vs free [{pair[0]} vs {pair[1]}, {un}]: '+'  '.join(f'{k} {v:+.2f}' for k,v in E.items() if k!='free'))
# sensitivity to systematics (1.5", N vs S): band scale sigmas and per-line
P('\n######## sensitivity (1.5", N vs S, 1-0 + 3-2): ln BF vs free for (sig_b10, sig_b32, per-line)')
for sb,sl in (((0.05,0.05),0.03),((0.10,0.15),0.05),((0.2,0.3),0.05),((0.1,0.15),0.10),((0.3,0.3),0.10)):
    E=evid('1.50',('N','S'),USE['1-0 + 3-2'],sb=sb,sl=sl); res[f'sens|{sb}|{sl}']=E
    P(f'  sb={sb} sl={sl}: '+'  '.join(f'{k} {v:+.2f}' for k,v in E.items() if k!='free'))
open(O+'s186.txt','w').write('\n'.join(log)); json.dump(res,open(O+'s186.json','w'),indent=1,default=float)
