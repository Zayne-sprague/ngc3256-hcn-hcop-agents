# Q36 step 4: HCN/HCO+ RADEX-equivalent grids (same axes/format as out/q4c/grid4_*.npz) with the frequency-dependent background of q36_bg.
# usage: q36_grids.py  -> out/q36/grid_{mol}_{reg}_{scen}.npz  (and CMB-only control for self-consistency vs RADEX)
import sys,itertools,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, q36_rt as R, q36_bg as B
from concurrent.futures import ProcessPoolExecutor
TK=np.array([10,15,20,25,30,40,50,65,80,100,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25); FN={'hcn':'/opt/lamda/hcn.dat','hco+':'/opt/lamda/hco+.dat'}
def job(a):
    mol,reg,scen=a; M=R.parse(FN[mol]); f=(lambda nu: B.Jbg(reg,scen,nu/1e9)); out=np.full((len(TK),len(LN),len(LNC),4,3),np.nan)
    for i,T in enumerate(TK):
        for j,ln in enumerate(LN):
            for k,lN in enumerate(LNC):
                o=R.solve(M,T,10**ln,10**lN,150.,f); out[i,j,k,:,0]=o['TR'][:4]; out[i,j,k,:,1]=o['tau'][:4]; out[i,j,k,:,2]=o['Tex'][:4]
    np.savez(f'/workspace/out/q36/grid_{mol}_{reg}_{scen}.npz',TR=out[...,0],tau=out[...,1],Tex=out[...,2],TK=TK,LN=LN,LNC=LNC); return a
if __name__=='__main__':
    jobs=[(m,r,s) for m in FN for r in 'NS' for s in ('nominal','dusty','high','cmb')]+[(m,'E',s) for m in FN for s in B.ENVTB]
    with ProcessPoolExecutor(12) as ex:
        for a in ex.map(job,jobs): print('done',a,flush=True)
