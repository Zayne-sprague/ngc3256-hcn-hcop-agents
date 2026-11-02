# Q12c: posterior mass of the 'compact core thick' / 'both thick, abundance-matched' branches and the shape change from the H13CN term (two-comp shared-T model, r7, 3 seeds pooled).
# thick: tau32(HCN, compact) > 10 (computed on the RADEX grid as in q12c_twocomp.py); matched: |log10 N(HCN)/N(HCO+)|_c < 0.3; n_c>=1e6: log n_c>=6. Output out/q12c/branches.csv
import numpy as np, pandas as pd
from scipy.ndimage import map_coordinates
O='/workspace/out/q4c/'; D='/workspace/out/q12c/runs/'
G=np.load(O+'grid4_hcn.npz'); TK,LN,LNC=G['TK'],G['LN'],G['LNC']; TAU=G['tau'][...,2].astype(float)
def tau(lT,ln,lN): return map_coordinates(TAU,np.array([np.interp(lT,np.log10(TK),np.arange(len(TK))),np.interp(ln,LN,np.arange(len(LN))),np.interp(lN,LNC,np.arange(len(LNC)))]),order=1,mode='nearest')
rows=[]
for reg in 'NS':
  for model in ('1z','2cT'):
    for tag in ('base','geo','geo30','floor25'):
      r={}
      for pre in ('','h_'):
        X=[]
        for s in (1,2,3):
            try: z=np.load(D+f'{pre}{tag}_{model}_{reg}_r7_tbg2.73_s{s}.npz'); nm=list(z['names']); X.append(z['X'])
            except Exception: pass
        if not X: continue
        X=np.concatenate(X); t=tau(X[:,nm.index('lTc')],X[:,nm.index('nc')],X[:,nm.index('NHc')]); lx=X[:,nm.index('NHc')]-X[:,nm.index('NOc')]
        thick=t>10; match=np.abs(lx)<0.3; hi=X[:,nm.index('nc')]>=6
        r.update({f'{pre}P_thick10':thick.mean(),f'{pre}P_matched':match.mean(),f'{pre}P_thick&matched':(thick&match).mean(),f'{pre}P_nc>=1e6':hi.mean(),f'{pre}tau32_16/50/84':'%.1f/%.1f/%.1f'%tuple(np.percentile(t,[16,50,84]))})
      if r: rows.append(dict(region=reg,model=model,tag=tag,**r))
B=pd.DataFrame(rows); B.round(3).to_csv('/workspace/out/q12c/branches.csv',index=False); pd.set_option('display.width',250,'display.max_columns',30); print(B.round(3).to_string())
