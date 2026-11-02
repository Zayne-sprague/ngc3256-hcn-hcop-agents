# Q4c-2: RADEX grids (LAMDA, H2 collider, Tbg 2.73 K, uniform sphere, dv=150 km/s; N/dv scaling used for other dv) incl. J=1-0: HCN, HCO+ (1-0..4-3), CS (7-6, 8-7), H13CN (4-3), H13CO+ (not run).
# Output out/q4c/grid4_{mol}.npz  with TR,tau,Tex arrays (T,n,N,line) ; lines ordered by frequency.
import os,re,subprocess,itertools,sys,numpy as np
from concurrent.futures import ProcessPoolExecutor
TK=np.array([10,15,20,25,30,40,50,65,80,100,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25); DV=150.
CFG={'hcn':('/opt/lamda/hcn.dat','50 400',4),'hco+':('/opt/lamda/hco+.dat','50 400',4),'cs':('/opt/lamda/cs.dat','300 400',2),'h13cn':('/opt/lamda/h13cn@xpol.dat','340 350',1)}
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(a):
    mol,T,ln,lN=a; fn,rng,nl=CFG[mol]; wd=f'/tmp/rxq_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    subprocess.run(['radex'],input=f'{fn}\n{wd}/o.txt\n{rng}\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n',text=True,capture_output=True)
    out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out.append((float(p[3]),fl(p[5]),fl(p[6]),fl(p[7])))
    out=sorted(out)[:nl]
    while len(out)<nl: out.append((np.nan,)*4)
    return [(o[3],o[2],o[1]) for o in out]
if __name__=='__main__':
    for mol in CFG:
        path=f'/workspace/out/q4c/grid4_{mol}.npz'
        if os.path.exists(path): continue
        jobs=list(itertools.product([mol],TK,LN,LNC)); nl=CFG[mol][2]
        with ProcessPoolExecutor(4) as ex: r=list(ex.map(run,jobs,chunksize=200))
        a=np.array(r,float).reshape(len(TK),len(LN),len(LNC),nl,3); np.savez(path,TR=a[...,0],tau=a[...,1],Tex=a[...,2],TK=TK,LN=LN,LNC=LNC); print(mol,a.shape,np.isnan(a).sum(),flush=True)
