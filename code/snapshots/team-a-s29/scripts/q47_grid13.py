# Q47(c): RADEX H13CN J=1-0 grid (own LAMDA h13cn@xpol; Tbg 2.73 K, dv=150 km/s, uniform sphere), axes as in out/q36 grids (T,n) with N(H13CN) log 11..18.5 step 0.25.
# Output out/q47/grid_h13cn_10.npz (TR in K for the 1-0 line, tau).
import os,re,subprocess,itertools,numpy as np
from concurrent.futures import ProcessPoolExecutor
G=np.load('/workspace/out/q36/grid_hcn_N_cmb.npz'); TK=G['TK']; LN=G['LN']; LNC=np.arange(11,18.51,0.25); DV=150.
def run(a):
    T,ln,lN=a; wd=f'/tmp/rx47_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    subprocess.run(['radex'],input=f'/opt/lamda/h13cn@xpol.dat\n{wd}/o.txt\n80 90\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n',text=True,capture_output=True)
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit() and p[0]=='1':
            f=lambda x:float(re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)); return (f(p[7]),f(p[5]))
    return (np.nan,np.nan)
if __name__=='__main__':
    jobs=list(itertools.product(TK,LN,LNC))
    with ProcessPoolExecutor(8) as ex: r=list(ex.map(run,jobs,chunksize=100))
    a=np.array(r,float).reshape(len(TK),len(LN),len(LNC),2); os.makedirs('/workspace/out/q47',exist_ok=True)
    np.savez('/workspace/out/q47/grid_h13cn_10.npz',TR=a[...,0],tau=a[...,1],TK=TK,LN=LN,LNC=LNC); print(a.shape,np.isnan(a).sum())
