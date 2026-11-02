# Q40: add 13C-isotopologue constraints (H13CO+ 1-0, 3-2; H13CN 1-0) to the Q39 fine-grid one-gas RADEX fit (s220b machinery, imported and monkeypatched).
# 13C species = same Tk, n, eta, band scale as the 12C species, column N/X (X = 12C/13C = 40, 60, 90 assumed; the column shift is rounded to 0.05 dex: X_eff 39.8/63.1/89.1).
# Isotopologue data enter in LINEAR intensity (K km/s, so that non-detections are fine): chi = (I_mod*exp(o) - I_obs)^2/(e^2 + (0.10 I_obs)^2 [+ (0.10 I_mod)^2 for model/LAMDA-scaling]).
# Data (1.7" peak pixel matched to Harada Table 13 12C lines):
#  primary 'A': N H13CO+ 1-0 9.4+-0.95, 3-2 2.6+-0.49 (Harada 1.7"); N H13CN 1-0 = 0.030+-0.006 x HCN (Q29 s185, peak pixel of the 2.2" map)
#               S H13CO+ 1-0 0.014+-0.009 x HCO+ (Q29 peak pixel; not tabulated by Harada = not detected there), S H13CN 0.0057+-0.0059 x HCN (Q29 2.2" aperture fit)
#               C H13CO+ 1-0 3.44+-0.95 (Harada); C H13CN 0.0237+-0.0044 x HCN (Q29 Bridge 2.2" aperture)
#  variant 'B' (aperture ratios of Q29 for H13CO+ and H13CN at all positions: N 0.0231, 0.0207; S 0.0184,0.0057; C 0.0185,0.0237) with Harada 12C intensities
#  variant 'H' : Harada H13CO+ only (N 1-0+3-2; C 1-0; S none).  Variant 'N' = no isotopologue (reproduces Q39).
# usage: python3 -I s230_q40_isotopologue_fit.py  (about 10 min with 16 procs) -> out/q40/*.npz, out/q40/fit.log
import sys,os,json,time,numpy as np
sys.path.insert(0,'/workspace/scripts'); os.makedirs('/workspace/out/q40',exist_ok=True)
import s220b_q39_fine_fit as m
from scipy.interpolate import CubicSpline
from multiprocessing import Pool
O='/workspace/out/q40/'
NF,LNC=m.NF,m.LNC; d=m.d
def lnI13(a,logX):
    TR=np.nan_to_num(a[...,:2,0],nan=1e-3); l=np.log(1.0645*150*np.maximum(TR,1e-3))
    return CubicSpline(LNC,l,axis=2,extrapolate=True)(np.clip(NF-logX,11.5,19))
H=m.H
def hv(p,s): r=H[(H.pos==p)&(H.sp.str.startswith(s))].iloc[0]; return r.I_Kkms
def iso(variant):
    """returns {pos:{species:[(lineindex 0=1-0 band3,1=3-2 band6, Iobs, err)]}}"""
    q={'N':dict(O=[(0,9.4,0.95),(1,2.6,0.49)],H=[(0,0.030*hv('N','HCN1'),0.006*hv('N','HCN1'))]),
       'S':dict(O=[(0,0.014*hv('S','HCO+1'),0.009*hv('S','HCO+1'))],H=[(0,0.0057*hv('S','HCN1'),0.0059*hv('S','HCN1'))]),
       'C':dict(O=[(0,3.44,0.95)],H=[(0,0.0237*hv('C','HCN1'),0.0044*hv('C','HCN1'))])}
    if variant=='B':
        for p,(ro,rh,eo,eh) in {'N':(0.0231,0.0207,0.0031,0.0034),'S':(0.0184,0.0057,0.0040,0.0059),'C':(0.0185,0.0237,0.0033,0.0044)}.items():
            q[p]['O']=[(0,ro*hv(p,'HCO+1'),eo*hv(p,'HCO+1'))]+([(1,ro*hv(p,'HCO+3'),eo*hv(p,'HCO+3'))] if p!='S' and False else [])
            q[p]['H']=[(0,rh*hv(p,'HCN1'),eh*hv(p,'HCN1'))]
    if variant=='H':
        q['S']=dict(O=[],H=[]); q['N']['H']=[]; q['C']['H']=[]
    if variant=='N': q={p:dict(O=[],H=[]) for p in q}
    return q
MOD13={}
def chi_s(it,jn,s,o3=None,o6=None):
    f=m.chi_s0(it,jn,s,o3,o6)
    if f is None: return None
    for (li,I,e) in X_iso[s]:
        M=MOD13[s][it,jn][:,li]; o=m.X.o3 if li==0 else m.X.o6
        mod=np.exp(M[None,:]+o[:,None]); sig2=e**2+(0.10*I)**2+(0.10*mod)**2
        f=f+(mod-I)**2/sig2*(1.0)
    return f
m.chi_s0=m.chi_s; m.chi_s=chi_s
def run_one(pos,logX,variant,sig=(0.10,0.15)):
    global X_iso
    q=iso(variant)[pos]; X_iso={'H':[(li,I,e) for li,I,e in q['H']],'O':[(li,I,e) for li,I,e in q['O']],'C':[]}
    MOD13['H']=lnI13(d['hcn'],logX); MOD13['O']=lnI13(d['hcop'],logX)
    # note 1.0645*150 conversion identical to s220b (K km/s from T_R with dv=150 km/s)
    cells=[(i,j) for i in range(len(m.TK)) for j in range(len(m.LN))]
    m.setup(pos,sig[0],sig[1],['H','O','C'])
    with Pool(min(16,os.cpu_count())) as P:
        mm=np.array(P.map(m.cell_min,cells,chunksize=20)); cb=mm.min()
        res=P.map(m.cell_pass2,[(c,cb,45.) for c in cells],chunksize=20)
    Wc=np.array([r['Wc'] if r else 0 for r in res]).reshape(len(m.TK),len(m.LN))
    out=dict(chibest=cb,Wc=Wc,chimin=mm.reshape(len(m.TK),len(m.LN)))
    for key in m.LAGS:
        hs=[r['hist'][key] for r in res if r and key in r['hist']]
        out['h_'+key]=np.sum(hs,0); out['p_'+key]=np.min([r['prof'][key] for r in res if r and key in r['prof']],0)
    return out
if __name__=='__main__':
    variants=sys.argv[1].split(',') if len(sys.argv)>1 else ['N','A','B','H']
    Xs=[40,60,90]
    for v in variants:
      for X in (Xs if v!='N' else [60]):
        for pos in 'NSC':
            t=time.time(); o=run_one(pos,np.log10(X),v); tag=f'{pos}_{v}_X{X}'
            np.savez(O+f'fit_{tag}.npz',**o)
            cm=o['chimin']; Tk=m.TK; hot=cm[Tk>=100].min(); cold=cm[Tk<=30].min()
            print(tag,f"chi2min {o['chibest']:.2f}  hot(Tk>=100) {hot:.2f}  cold(Tk<=30) {cold:.2f}  dchi2 cold-hot {cold-hot:.2f}  {time.time()-t:.0f}s",flush=True)
