# validate q36_rt against the team RADEX grids (Tbg 2.73, 5, 20 K) over the full Q22 prior box; report agreement fraction
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts'); import numpy as np, q36_rt as R
for mol,fn in (('hcn','/opt/lamda/hcn.dat'),('hco+','/opt/lamda/hco+.dat')):
    M=R.parse(fn)
    for suf,Tb in (('',2.73),('_tbg5',5.),('_tbg20',20.)):
        G=np.load(f'/workspace/out/q4c/grid4_{mol}{suf}.npz'); err=[]; bad=[]
        for iT,T in enumerate(G['TK'][::2]):
          for iN,ln in enumerate(G['LN'][::2]):
            for lN in G['LNC'][::3]:
                ref=G['TR'][2*iT,2*iN,list(G['LNC']).index(lN)]
                ok=np.isfinite(ref)&(ref>0.05)&(ref<1e3)
                if not ok.all(): continue
                o=R.solve(M,T,10**ln,10**lN,150.,Tb)['TR'][:4]
                e=np.abs(o/ref-1).max(); err.append(e)
                if e>0.05: bad.append((T,ln,lN,round(e,2)))
        err=np.array(err); print(mol,Tb,'n=%d median %.4f  frac<2%% %.3f frac<5%% %.3f'%(len(err),np.median(err),(err<.02).mean(),(err<.05).mean()),bad[:6])
