# Q47: summarise N-only NS runs with/without the 13C ratios in the likelihood (out/q47/runs): per-seed T quantiles, P(T<=30), P(T>=40), logZ, ESS; mean/sd across seeds.
import glob,json,numpy as np,pandas as pd
def wq(x,w,q): o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return np.interp(q,c,x[o])
rows=[]
for f in sorted(glob.glob('/workspace/out/q47/runs/q47_N_*.npz')):
    d=np.load(f,allow_pickle=True); j=json.load(open(f[:-4]+'.json')); nm=list(d['names']); T=10**d['S'][:,nm.index('lTc_N')]; w=d['w']; q=wq(T,w,[.16,.5,.84])
    rows.append(dict(mode=j['mode'],data=j['data'],seed=j['seed'],logZ=j['logZ'],ESS=j['ess'],T16=q[0],T50=q[1],T84=q[2],P_le30=w[T<=30].sum(),P_ge40=w[T>=40].sum(),chi2min=j['chi2min']))
df=pd.DataFrame(rows); df.round(3).to_csv('/workspace/out/q47/ns13_runs.csv',index=False); pd.set_option('display.width',200)
print(df.round(3).to_string()); print(df.groupby(['mode','data'])[['logZ','T50','P_le30','P_ge40']].agg(['mean','std']).round(3).to_string())
