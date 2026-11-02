"""Q54 r2 fix2: RADEX-equivalent (q36_rt, CMB background, 150 km/s, H2 partner (CN: He x1.37 proxy, LAMDA)) grids for HNC and CN on the Q36/Q47 axes.
HNC: transitions 1-0..4-3 (first 4). CN: transitions (0-based) 0,1 = N=1-0 J=1/2,3/2 groups; 6,7 = N=3-2 J=5/2-3/2, 7/2-5/2 groups (hyperfine unresolved: sum of components in one 150 km/s line).
usage: q54h_grids.py -> out/q54/grid_{hnc,cn}_cmb.npz"""
import sys,warnings,time; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, q36_rt as R
from concurrent.futures import ProcessPoolExecutor
TK=np.array([10,15,20,25,30,40,50,65,80,100,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25)
FN={'hnc':('/workspace/out/q54/hnc.norm.dat',[0,1,2,3]),'cn':('/workspace/out/q54/cn.norm.dat',[0,1,5,6,7])}
def job(a):
    mol,i=a; fn,idx=FN[mol]; M=R.parse(fn); T=TK[i]; out=np.full((len(LN),len(LNC),len(idx),2),np.nan)
    for j,ln in enumerate(LN):
        for k,lN in enumerate(LNC):
            o=R.solve(M,T,10**ln,10**lN,150.,2.725); out[j,k,:,0]=o['TR'][idx]; out[j,k,:,1]=o['tau'][idx]
    return mol,i,out
if __name__=='__main__':
    t0=time.time(); jobs=[(m,i) for m in FN for i in range(len(TK))]; res={m:np.full((len(TK),len(LN),len(LNC),len(FN[m][1]),2),np.nan) for m in FN}
    with ProcessPoolExecutor(16) as ex:
        for m,i,o in ex.map(job,jobs): res[m][i]=o; print('done',m,i,round(time.time()-t0),flush=True)
    for m in FN: np.savez(f'/workspace/out/q54/grid_{m}_cmb.npz',TR=res[m][...,0],tau=res[m][...,1],TK=TK,LN=LN,LNC=LNC)
