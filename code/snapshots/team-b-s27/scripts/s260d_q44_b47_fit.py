# Q44(c): one-zone RADEX fit (s220b machinery, via s260b) of HCN, HCO+ 1-0, 3-2, 4-3 at the Xb47 beam (4.14x2.70", s260c.json), 12C only (no 13C at this beam), no CS.
# 4-3 lines (Band 7) use the band-6 offset o6 plus an extra independent 15% (Band-7 scale) in the variance.  Variants (VAR): 
#   B6  : r=0.6" beam-peak, lines 1-0,3-2,4-3 (6 lines)      B6n : same without 4-3 (1-0,3-2 only; like-for-like control)
#   B22 : r=2.2" aperture, 6 lines                            B22n: same w/o 4-3
# usage: python3 -I s260d_q44_b47_fit.py POS[,POS] VAR [k...]   -> out/q44/fit_<pos>_<var>_s<k>_lX160.npz
import sys,os,json,time,numpy as np
sys.path.insert(0,'/workspace/scripts'); _a=sys.argv; sys.argv=['x']
import s260b_q44_fit as F
sys.argv=_a
from scipy.interpolate import CubicSpline
B=F.B; m=F.m; s=F.s; X=m.X; S235=F.S235
R=json.load(open('/workspace/out/q44/s260c.json'))
def lnI3(a):
    TR=np.nan_to_num(a[...,:3,0],nan=1e-3); l=np.log(1.0645*150*np.maximum(TR,1e-3)); return CubicSpline(m.LNC,l,axis=2)(m.NF)
MOD3={'H':lnI3(m.d['hcn']),'O':lnI3(m.d['hcop'])}
VAR={'B6':('0.6',True),'B6n':('0.6',False),'B22':('2.2',True),'B22n':('2.2',False)}
CF={'var':None,'pos':None}
def obs(pos,rad):
    rn={'N':'N','S':'S','C':'Bridge'}[pos]; g=lambda k: R[f'{rn}|{rad}|{k}']
    return {'H':[(0,'3',g('HCN1')),(1,'6',g('HCN3')),(2,'7',g('HCN4'))],'O':[(0,'3',g('HCO1')),(1,'6',g('HCO3')),(2,'7',g('HCO4'))]}
orig0=m.chi_s0
def chi_new(it,jn,sp,o3=None,o6=None):
    if CF['var'] is None: return orig0(it,jn,sp,o3,o6)
    if sp not in MOD3: return None
    rad,use4=VAR[CF['var']]; M=MOD3[sp][it,jn]
    if np.isnan(M).any(): return None
    o3=X.o3 if o3 is None else o3; o6=X.o6 if o6 is None else o6
    f=0
    for li,b,(I,e) in obs(CF['pos'],rad)[sp]:
        if li==2 and not use4: continue
        y=np.log(I); v=(e/I)**2+0.05**2+(0.01**2) + (0.15**2 if b=='7' else 0)
        off=o3 if b=='3' else o6
        f=f+(y-M[None,:,li]-off[:,None])**2/v
    return f
m.chi_s0=chi_new
def run(pos,var,k):
    CF['var']=var; CF['pos']=pos; B.CUR['terms']=[]; B.RAT.clear(); F.CFG['var']='P'; F.CFG['pos']=pos
    return F.run_one(pos,1.6,'N',S235.SIGS[k])
if __name__=='__main__':
    poss=sys.argv[1].split(','); var=sys.argv[2]; ks=[int(x) for x in sys.argv[3:]] or [0,1,2]
    for pos in poss:
        for k in ks:
            fn=F.O+f'fit_{pos}_{var}_s{k}_lX160.npz'
            if os.path.exists(fn): continue
            t=time.time(); o=run(pos,var,k); np.savez(fn,**o); cm=o['chimin']; Tk=m.TK
            print(pos,var,S235.SIGS[k],f"chi2 {o['chibest']:.2f} cold-hot {cm[Tk<=30].min()-cm[Tk>=100].min():.2f} {time.time()-t:.0f}s",flush=True)
