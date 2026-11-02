# Q2 sensitivity: HCN/HCO+ RADEX grids (same T,n,N,dv grid as q4c_grid.py) with a different background blackbody Tbg (arg, K). Output out/q4c/grid4_{mol}_tbg{T}.npz
import sys,os,itertools,numpy as np
sys.path.insert(0,'/workspace/scripts')
import q4c_grid as g
from concurrent.futures import ProcessPoolExecutor
tb=float(sys.argv[1]); suf='_tbg%d'%tb
def run(a):
    mol,T,ln,lN=a; fn,rng,nl=g.CFG[mol]; wd=f'/tmp/rxb_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    import subprocess
    subprocess.run(['radex'],input=f'{fn}\n{wd}/o.txt\n{rng}\n{T}\n1\nH2\n{10**ln:.4e}\n{tb}\n{10**lN:.4e}\n{g.DV}\n0\n',text=True,capture_output=True)
    out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out.append((float(p[3]),g.fl(p[5]),g.fl(p[6]),g.fl(p[7])))
    out=sorted(out)[:nl]
    while len(out)<nl: out.append((np.nan,)*4)
    return [(o[3],o[2],o[1]) for o in out]
if __name__=='__main__':
    for mol in ('hcn','hco+'):
        jobs=list(itertools.product([mol],g.TK,g.LN,g.LNC))
        with ProcessPoolExecutor(14) as ex: r=list(ex.map(run,jobs,chunksize=100))
        a=np.array(r,float).reshape(len(g.TK),len(g.LN),len(g.LNC),4,3)
        np.savez(f'/workspace/out/q4c/grid4_{mol}{suf}.npz',TR=a[...,0],tau=a[...,1],Tex=a[...,2],TK=g.TK,LN=g.LN,LNC=g.LNC); print(mol,suf,np.isnan(a).sum(),flush=True)
