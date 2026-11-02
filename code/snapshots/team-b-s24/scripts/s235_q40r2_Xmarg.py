# Q40 round 2: H13CO+ (Harada 1.7" 1-0, 3-2; beam-matched to the 12C lines) added to the Q39 one-gas RADEX fit, with
#  (i) 12C/13C X stepped on a fine log grid (0.05 dex, X=22..100) so that it can be marginalised with a prior (done in s236),
#  (ii) three band-scale prior widths (s3,s6) = (0.10,0.15) [Q39 default], (0.20,0.25), (0.40,0.40) [~free],
#  (iii) optional extra fractional systematic FR on the 13C lines (s230 already has 10% of I + 10% of model; FR adds in quadrature).
# Variants: 'H' = H13CO+ only (matched 1.7"); 'A2' = H13CO+ + H13CN with the H13CN (2.2" Q29 maps scaled to Harada 12C) given a 30% beam-mismatch systematic.
# usage: python3 -I s235_q40r2_Xmarg.py POS[,POS] VARIANT   -> out/q40r2/fit_<pos>_<var>_s<k>_lX<logX*100>.npz (+ stdout log)
import sys,os,numpy as np,time
sys.path.insert(0,'/workspace/scripts'); _a=sys.argv; sys.argv=['x']
import s230_q40_isotopologue_fit as s
sys.argv=_a
m=s.m; O='/workspace/out/q40r2/'
FR={'H':0.0,'A2':0.0}
BEAMSYS={}            # species -> extra fractional systematic
def chi_s(it,jn,sp,o3=None,o6=None):
    f=m.chi_s0(it,jn,sp,o3,o6)
    if f is None: return None
    for (li,I,e) in s.X_iso[sp]:
        M=s.MOD13[sp][it,jn][:,li]; o=m.X.o3 if li==0 else m.X.o6
        mod=np.exp(M[None,:]+o[:,None]); fr=BEAMSYS.get(sp,0.0)
        s2=e**2+(0.10*I)**2+(0.10*mod)**2; s2o=s2+(fr*I)**2+(fr*mod)**2
        f=f-(mod-I)**2/s2+(mod-I)**2/s2o     # s.chi_s already added the base term; swap in the augmented variance
    return f
# NB: m.chi_s is s.chi_s (adds base 13C term). Wrap it.
base=m.chi_s
def wrapped(it,jn,sp,o3=None,o6=None):
    f=base(it,jn,sp,o3,o6)
    if f is None or not BEAMSYS.get(sp): return f
    for (li,I,e) in s.X_iso[sp]:
        M=s.MOD13[sp][it,jn][:,li]; o=m.X.o3 if li==0 else m.X.o6
        mod=np.exp(M[None,:]+o[:,None]); fr=BEAMSYS[sp]
        s2=e**2+(0.10*I)**2+(0.10*mod)**2; s2o=s2+(fr*I)**2+(fr*mod)**2
        f=f-(mod-I)**2/s2+(mod-I)**2/s2o
    return f
m.chi_s=wrapped
SIGS=[(0.10,0.15),(0.20,0.25),(0.40,0.40)]
LOGX=np.round(np.arange(1.35,2.0001,0.05),2)
def iso2(variant,pos):
    q=s.iso('A' if variant=='A2' else 'H')[pos]; return q
if __name__=='__main__':
    poss=sys.argv[1].split(','); var=sys.argv[2]
    if var=='A2': BEAMSYS['H']=0.30
    sigs=SIGS if var=='H' else SIGS[:2]
    for pos in poss:
      for k,sg in enumerate(sigs):
        for lx in LOGX:
            t=time.time(); fn=O+f'fit_{pos}_{var}_s{k}_lX{int(round(lx*100))}.npz'
            if os.path.exists(fn): continue
            # run_one reads s.iso(variant) -> use 'A' for A2
            q=iso2(var,pos); s.X_iso={'H':q['H'],'O':q['O'],'C':[]}
            vv='A' if var=='A2' else 'H'
            o=s.run_one(pos,lx,vv,sg); np.savez(fn,**o)
            cm=o['chimin']; Tk=m.TK
            print(pos,var,sg,lx,f"chi2 {o['chibest']:.2f} cold-hot {cm[Tk<=30].min()-cm[Tk>=100].min():.2f} {time.time()-t:.0f}s",flush=True)
