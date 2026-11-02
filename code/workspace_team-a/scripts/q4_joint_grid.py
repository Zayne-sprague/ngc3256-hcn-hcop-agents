# Q4b-2: RADEX grids (LAMDA hcn.dat, hco+.dat; H2 partner, Tbg 2.73 K, uniform sphere, dv=150 km/s) with Tex, tau and T_R for J=2-1,3-2,4-3. Axes (Tkin, log n, log N). Output out/q4b/grid_{hcn,hco+}.npz
import os,re,subprocess,itertools,numpy as np
from concurrent.futures import ProcessPoolExecutor
TK=np.array([10,15,20,25,30,40,50,65,80,100,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25); DV=150.
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(a):
    mol,T,ln,lN=a; wd=f'/tmp/rxj_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    subprocess.run(['radex'],input=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n150 400\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n',text=True,capture_output=True)
    out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out.append((int(p[1]),fl(p[5]),fl(p[6]),fl(p[7])))
    out=sorted(out)[:3]; return [(o[3],o[2],o[1]) for o in out]
if __name__=='__main__':
    for mol in ('hcn','hco+'):
        path=f'/workspace/out/q4b/grid_{mol}.npz'
        if os.path.exists(path): continue
        jobs=list(itertools.product([mol],TK,LN,LNC))
        with ProcessPoolExecutor(8) as ex: r=list(ex.map(run,jobs,chunksize=100))
        a=np.array(r,float).reshape(len(TK),len(LN),len(LNC),3,3); np.savez(path,TR=a[...,0],tau=a[...,1],Tex=a[...,2],TK=TK,LN=LN,LNC=LNC); print(mol,a.shape,np.isnan(a).sum())
