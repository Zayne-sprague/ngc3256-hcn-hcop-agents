# Q55d: low-column RADEX grids (HCN, HCO+ collisional data, same T/n grid as q40h) used as 13C isotopologues:
# N(H13CN)=N(HCN)/X13 etc. Output out/q55/radex_low_{hcn,hco+}.npy: shape (nT,nn,nN,4 J,2 [TR,tau]) for logN = 10.0..12.75 (hcn), 10.0..12.75 (hco+)
import sys,os,itertools; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import q40h_radex as R
LNL=np.arange(10.0,12.76,0.25)
os.makedirs('/workspace/out/q55',exist_ok=True)
if __name__=='__main__':
    for mol in ('hcn','hco+'):
        with ProcessPoolExecutor(12) as ex: res=list(ex.map(R.run,list(itertools.product([mol],R.TK,R.LN,LNL)),chunksize=40))
        a=np.array(res).reshape(len(R.TK),len(R.LN),len(LNL),4,2); np.save(f'/workspace/out/q55/radex_low_{mol}.npy',a); print(mol,a.shape,np.isnan(a).sum())
