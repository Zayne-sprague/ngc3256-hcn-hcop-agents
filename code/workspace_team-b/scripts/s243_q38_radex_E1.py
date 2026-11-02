"""Q38 fix 1: RADEX (s186/s193 evidence machinery, LVG LAMDA grid) contrast E1 vs N with the UV-COMPLETE MATCHED 2.2" data only: 1-0 (X961 12m) + 3-2 (Q34 12m+7m comb), both lines, K km/s from s240.json block C
(r=2.2" apertures re-centred on the HCO+ peak).  Replaces the 10" (blended) ladder as the basis of the E1-vs-centre Tk contrast.  Only 4 data per position (HCN/HCO+ x 1-0/3-2): density and Tk are weakly constrained,
the test is the same-Tk hypothesis BF and P(Tk_E1<Tk_N) with shared band-scale draws (10% 1-0, 15% 3-2), 5% per-line, Tk<=300 K.  usage: python3 -I s243_q38_radex_E1.py -> out/q38/s243.txt/json"""
import sys, json, numpy as np
src=open('/workspace/scripts/s186_q29_radex_evidence.py').read(); src=src[:src.index("res={}")]
exec(src)
C=json.load(open('/workspace/out/q38/s240.json'))['C']
def data(th,reg):
    a=C[f'1-0 X961 2.2"|{reg}']; b=C[f'3-2 comb 2.2"|{reg}']
    I=np.array([a['HCN'][0],b['HCN'][0],a['HCOp'][0],b['HCOp'][0]])*K; e=np.array([a['HCN'][1],b['HCN'][1],a['HCOp'][1],b['HCOp'][1]])*K
    return I,e/I
out={}
USE={'1-0 + 3-2':np.array([True,True,True,True])}
T=10**F.lT
for g in ('N','E1','W','E2'):
    y,e=data('2.20',g); P(f' {g}: HCN1 {y[0]:.1f} HCN3 {y[1]:.1f} HCO1 {y[2]:.1f} HCO3 {y[3]:.1f} rel.stat {np.round(e,3).tolist()}')
    p=post1('2.20',g,USE['1-0 + 3-2']); out[f'post|{g}']=p
    P(f'  {g}: Tk {p["Tk"][1]:.0f} [{p["Tk"][0]:.0f},{p["Tk"][2]:.0f}] P(Tk>50)={p["PTk50"]:.2f}; log n {p["logn"][1]:.2f} [{p["logn"][0]:.2f},{p["logn"][2]:.2f}]; log A {p["logA"][1]:.2f} [{p["logA"][0]:.2f},{p["logA"][2]:.2f}]; log N(HCO+) {p["logNO"][1]:.1f}')
def pcontrast(reg,nd=41,sb=(0.10,0.15),sl=0.05,seed=2):
    u=qmc.Sobol(d=2,scramble=True,seed=seed).random(64)[:nd]; z=norm.ppf(np.clip(u,1e-9,1-1e-9))*np.array(sb)
    d1,e1=data('2.20','N'); d2,e2=data('2.20',reg); use=USE['1-0 + 3-2']; pt=[]; pn=[]; pa=[]
    for zz in z:
        off=np.array([zz[0],zz[1],zz[0],zz[1]]); ws=[]
        for y,e in ((d1,e1),(d2,e2)):
            l=lnl(np.log(y)-off,e,use,sl); w=np.where(valid,np.exp(l-l.max()),0); ws.append(w/w.sum())
        # P(Tk_reg < Tk_N) etc. over independent joint within the draw
        def pless(x):
            o=np.argsort(x); xs=x[o]; w1=ws[0][o]; w2=ws[1][o]; c1=np.cumsum(w1)
            # P(x2 < x1) = sum_j w2_j * P(x1 > x2_j)
            return float(np.sum(w2*(1-c1)))
        pt.append(pless(T)); pn.append(pless(F.ln)); pa.append(pless(F.lA))
    return np.mean(pt),np.mean(pn),np.mean(pa)
for reg in ('E1','W'):
    a,b,c=pcontrast(reg); out[f'P|{reg}']=[a,b,c]
    P(f' P(Tk_{reg}<Tk_N)={a:.2f}  P(n_{reg}<n_N)={b:.2f}  P(logA_{reg}<logA_N)={c:.2f}   (shared band-scale draws; ties count as half-open)')
for un in ('1-0 + 3-2',):
    for reg in ('E1','W'):
        E=evid('2.20',('N',reg),USE[un]); out[f'ev|N-{reg}']=E
        P(f' ln BF vs free [N vs {reg}, {un}]: '+'  '.join(f'{k} {v:+.2f}' for k,v in E.items() if k!='free'))
open('/workspace/out/q38/s243.txt','w').write('\n'.join(log)); json.dump(out,open('/workspace/out/q38/s243.json','w'),indent=1,default=float)
