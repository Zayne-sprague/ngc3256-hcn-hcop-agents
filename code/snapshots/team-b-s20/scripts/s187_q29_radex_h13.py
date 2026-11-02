"""Q29 s187: RADEX ladder evidence (s186 machinery, 2.2" beam, X961 1-0 + Xb4d 3-2) with the H13CO+ 1-0 / HCO+ 1-0 integral ratio as an additional likelihood term:
R_obs = I(HCO+)/I(H13CO+) (s185 template amplitudes, X95d spw0+spw1 combined), model R = (1-exp(-tau))/(1-exp(-tau/X)) with tau = RADEX tau(HCO+ 1-0) of the grid cell, X=12C/13C log-flat 40-100 (marginalised),
sigma(lnR) = hypot(stat/null error, 0.08 cross-project flux scale).  Also reports the posterior with and without the term (Tk, n, A, tau, N(HCO+), N(HCN)) and ln BF N vs S / N vs Bridge / S vs Bridge.
Variants: X range 25-100 (low 12C/13C like Harada R=22), sigma inflated to 0.4.  usage: python3 -I s187_q29_radex_h13.py -> out/q29/s187.txt/json"""
import sys, json, numpy as np
src=open('/workspace/scripts/s186_q29_radex_evidence.py').read(); src=src[:src.index("res={}")]
src=src.replace("O='/workspace/out/q29/'; log=[]","O='/workspace/out/q29/'; log=[]")
exec(src)
J=json.load(open(O+'s185.json')); T=J['table']
def Robs(reg):
    a=T[f'{reg}|H13COp|comb']; ref=J['HCOp_int'][reg][0]
    return ref/a['amp'], np.hypot(a['err']/a['amp'],0.08)
def lnLR(reg,Xlo=40,Xhi=100,inflate=1.0):
    R,s=Robs(reg); tau=10**F.tauO[:,0]; X=np.exp(np.linspace(np.log(Xlo),np.log(Xhi),9)); out=[]
    for x in X:
        Rm=np.expm1(-tau)/np.expm1(-tau/x); out.append(-0.5*((np.log(Rm)-np.log(R))/(s*inflate))**2)
    r=np.logaddexp.reduce(np.array(out),axis=0)-np.log(len(X)); return np.where(np.isfinite(r),r,-1e30)
def lnl2(y,rel,use,sl,extra):
    return lnl(y,rel,use,sl)+extra
def post2(th,reg,use,extra,sb=(0.1,0.15),sl=0.05,nd=41,seed=2):
    u=qmc.Sobol(d=2,scramble=True,seed=seed).random(64)[:nd]; z=norm.ppf(np.clip(u,1e-9,1-1e-9))*np.array(sb)
    y,e=data(th,reg); w=np.zeros(F.nm)
    for zz in z:
        off=np.array([zz[0],zz[1],zz[0],zz[1]]); l=lnl2(np.log(y)-off,e,use,sl,extra); w+=np.where(valid,np.exp(l-l.max()),0)
    w/=w.sum()
    def q(x,p=(16,50,84)):
        o=np.argsort(x); c=np.cumsum(w[o]); return [float(x[o][min(np.searchsorted(c,a/100),len(x)-1)]) for a in p]
    Tk=10**F.lT
    return dict(Tk=q(Tk),logn=q(F.ln),logA=q(F.lA),logNO=q(F.lNO),logNH=q(F.lNH),tauHCO10=q(F.tauO[:,0]),tauHCN10=q(F.tauH[:,0]),PTk50=float(w[Tk>50].sum()))
def evid2(th,pair,use,extras,nd=41,sb=(0.10,0.15),sl=0.05,seed=2):
    u=qmc.Sobol(d=2,scramble=True,seed=seed).random(64)[:nd]; z=norm.ppf(np.clip(u,1e-9,1-1e-9))*np.array(sb)
    d1,e1=data(th,pair[0]); d2,e2=data(th,pair[1]); acc={k:[] for k in ('free','Tn','A','ANO','T','n','all')}
    for zz in z:
        off=np.array([zz[0],zz[1],zz[0],zz[1]]); Ds=[]
        for y,e,ex in ((d1,e1,extras[0]),(d2,e2,extras[1])):
            l=lnl2(np.log(y)-off,e,use,sl,ex); Ds.append(dense(l-l.max()))
        A1,A2=Ds; m=lambda X,ax: X.mean(axis=ax)
        acc['free'].append(A1.mean()*A2.mean()); acc['Tn'].append((m(A1,(2,3))*m(A2,(2,3))).mean()); acc['A'].append((m(A1,(0,1,2))*m(A2,(0,1,2))).mean())
        acc['ANO'].append((m(A1,(0,1))*m(A2,(0,1))).mean()); acc['T'].append((m(A1,(1,2,3))*m(A2,(1,2,3))).mean()); acc['n'].append((m(A1,(0,2,3))*m(A2,(0,2,3))).mean()); acc['all'].append((A1*A2).mean())
    return {k:float(np.log(np.mean(v))-np.log(np.mean(acc['free']))) for k,v in acc.items()}
res={}; ALL=np.array([True]*4); th='2.20'
P('H13CO+ constraint: R=I(HCO+)/I(H13CO+) (template integrals, 2.2" aperture):')
for g in ('N','S','Bridge'):
    R,s=Robs(g); P(f'  {g}: R = {R:.1f}  (ln-sigma {s:.2f})  -> {R*np.exp(-s):.0f}-{R*np.exp(s):.0f}')
VAR={'X40-100':dict(Xlo=40,Xhi=100,inflate=1.0),'X25-100':dict(Xlo=25,Xhi=100,inflate=1.0),'X40-100, sigma x2.5':dict(Xlo=40,Xhi=100,inflate=2.5)}
for vn,kw in VAR.items():
    P(f'\n#### variant {vn}')
    ex={g:lnLR(g,**kw) for g in ('N','S','Bridge')}
    for g in ('N','S','Bridge'):
        p0=post2(th,g,ALL,0.0); p1=post2(th,g,ALL,ex[g]); res[f'post|{vn}|{g}']=dict(no_H13=p0,with_H13=p1)
        for lab,p in (('without',p0),('with   ',p1)):
            P(f'  {g} {lab} H13CO+: Tk {p["Tk"][1]:.0f} [{p["Tk"][0]:.0f},{p["Tk"][2]:.0f}] P(Tk>50)={p["PTk50"]:.2f}; log n {p["logn"][1]:.2f} [{p["logn"][0]:.2f},{p["logn"][2]:.2f}]; log A {p["logA"][1]:.2f} [{p["logA"][0]:.2f},{p["logA"][2]:.2f}]; tau(HCO+ 1-0) {10**p["tauHCO10"][1]:.1f} [{10**p["tauHCO10"][0]:.1f},{10**p["tauHCO10"][2]:.1f}]; tau(HCN 1-0) {10**p["tauHCN10"][1]:.1f} [{10**p["tauHCN10"][0]:.1f},{10**p["tauHCN10"][2]:.1f}]; log N(HCO+) {p["logNO"][1]:.2f} [{p["logNO"][0]:.2f},{p["logNO"][2]:.2f}]; log N(HCN) {p["logNH"][1]:.2f} [{p["logNH"][0]:.2f},{p["logNH"][2]:.2f}]')
    for pair in (('N','S'),('N','Bridge'),('S','Bridge')):
        E0=evid2(th,pair,ALL,(0.0,0.0)); E1=evid2(th,pair,ALL,(ex[pair[0]],ex[pair[1]])); res[f'ev|{vn}|{pair[0]}-{pair[1]}']=dict(no_H13=E0,with_H13=E1)
        P(f'  ln BF vs free [{pair[0]} vs {pair[1]}] without: '+'  '.join(f'{k} {v:+.2f}' for k,v in E0.items() if k!='free'))
        P(f'  ln BF vs free [{pair[0]} vs {pair[1]}] WITH   : '+'  '.join(f'{k} {v:+.2f}' for k,v in E1.items() if k!='free'))
P('\n#### sensitivity of N vs S ln BF (with H13CO+, X40-100) to band-scale and per-line noise')
ex={g:lnLR(g) for g in ('N','S')}
for sb,sl in (((0.05,0.05),0.03),((0.10,0.15),0.05),((0.2,0.3),0.05),((0.1,0.15),0.10),((0.3,0.3),0.10)):
    E=evid2(th,('N','S'),ALL,(ex['N'],ex['S']),sb=sb,sl=sl); res[f'sens|{sb}|{sl}']=E
    P(f'  sb={sb} sl={sl}: '+'  '.join(f'{k} {v:+.2f}' for k,v in E.items() if k!='free'))
open(O+'s187.txt','w').write('\n'.join(log)); json.dump(res,open(O+'s187.json','w'),indent=1,default=float)
