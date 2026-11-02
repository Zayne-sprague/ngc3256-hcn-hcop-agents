"""Q53 r2: single-zone RADEX (q36 'cmb' grids, 150 km/s) HCN/HCO+ ladder over (Tk, n, N_HCO+, X=N_HCN/N_HCO+):
(1) how the ratio-of-ratios R43/R32 depends on R32 over the model space (is a constant R43/R32 an allowed extrapolation?),
(2) cells reproducing the observed N and S R10 (+R32) at one beam -> predicted R21, R43.
R_J = TR_HCN(J)/TR_HCO+(J) in K. Output out/q53/radex_ratio.json"""
import numpy as np, json
def L(m,r): 
    G=np.load(f'/workspace/out/q36/grid_{m}_{r}_cmb.npz'); return G
Gh,Go=L('hcn','N'),L('hco+','N')
TK,LN,LNC=Gh['TK'],Gh['LN'],Gh['LNC']; dN=LNC[1]-LNC[0]
TRh,TRo=Gh['TR'],Go['TR']            # (T,n,N,J)
tauh=Gh['tau']
res={}
# build all combos: for each (T,n): N_HCO index k, HCN index k+shift (logX = shift*dN)
shifts=np.arange(-4,9)   # logX = -1 .. +2
rows=[]
for i,T in enumerate(TK):
  for j,n in enumerate(LN):
    for k,lNo in enumerate(LNC):
      for s in shifts:
        kk=k+s
        if kk<0 or kk>=len(LNC): continue
        a=TRh[i,j,kk,:]; b=TRo[i,j,k,:]
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and np.all(b>1e-3) and a[0]>1e-3): continue
        rows.append((T,n,lNo,s*dN,*(a/b),tauh[i,j,kk,0],tauh[i,j,kk,2],a[2]))
A=np.array(rows); T,n,lNo,lX,R10,R21,R32,R43,t10,t32,Th32=A.T
print(len(A),'cells')
# restrict to physically plausible region: n 1e4-1e7, T 15-150, HCN TR(3-2)>1 K (detectable)
ok=(n>=4)&(n<=7)&(T>=15)&(Th32>0.5)
res['ncell']=int(ok.sum())
# (1) R43/R32 vs R32 bins
bins=[(0.3,0.5),(0.5,0.8),(0.8,1.2),(1.2,2.0)]
res['R43_over_R32_by_R32bin']={}
for lo,hi in bins:
    m=ok&(R32>=lo)&(R32<hi)
    q=np.percentile(R43[m]/R32[m],[5,16,50,84,95]) if m.sum()>20 else None
    res['R43_over_R32_by_R32bin'][f'{lo}-{hi}']=dict(n=int(m.sum()),pct5_16_50_84_95=None if q is None else [round(float(x),3) for x in q])
# (2) fits
def fit(name,obs,errs,extra=None):
    chi=sum(((arr-o)/e)**2 for arr,o,e in zip(*[ [R10 if k=='R10' else R21 if k=='R21' else R32 for k in obs],[obs[k] for k in obs],[errs[k] for k in obs]]))
    w=np.exp(-0.5*(chi-chi[ok].min()))*ok
    w/=w.sum()
    def wq(x,q): 
        o=np.argsort(x); c=np.cumsum(w[o]); return [float(x[o][np.searchsorted(c,qq)]) for qq in q]
    out=dict(chi2min=float(chi[ok].min()),ess=float(1/np.sum(w**2)),
        R21=wq(R21,[.16,.5,.84]),R43=wq(R43,[.16,.5,.84]),R43_over_R32=wq(R43/R32,[.16,.5,.84]),
        logT=wq(np.log10(T),[.16,.5,.84]),logn=wq(n,[.16,.5,.84]),logX=wq(lX,[.16,.5,.84]),tau10_HCN=wq(t10,[.16,.5,.84]))
    res[name]=out; print(name,json.dumps(out))
# N and S at the 1.7" (320 pc) beam: R10 Harada, R32 from out/q53/matched32.csv (1.7" interpolated), 10% model error
import pandas as pd
M=pd.read_csv('/workspace/out/q53/matched32.csv'); M=M[M.ous=='Xb4b']
def r32(reg,bm): g=M[M.reg==reg]; return float(np.interp(bm,g.beam_as,g.ratio))
for reg,r10 in (('N',0.961),('S',0.645)):
    o=dict(R10=r10,R32=r32(reg,1.7)); e=dict(R10=0.1*r10,R32=0.1*o['R32'])
    fit(f'{reg}_1p7',o,e)
# core 0.3": R32 only
for reg in 'NS':
    o=dict(R32=r32(reg,0.3)); e=dict(R32=0.1*o['R32']); fit(f'{reg}_core_R32only',o,e)
json.dump(res,open('/workspace/out/q53/radex_ratio.json','w'),indent=1)
print(json.dumps(res['R43_over_R32_by_R32bin'],indent=1))
