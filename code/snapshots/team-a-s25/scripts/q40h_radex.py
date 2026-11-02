# Q40h: RADEX grids (LAMDA hcn.dat, hco+.dat; H2 collider; Tbg 2.73 K; dv=150 km/s) for J=1-0..4-3, T_R and tau; then chi2 of (T,n,N(HCN),X=N(HCN)/N(HCO+)) vs R10,R21,R32 (T_B scale) of N,S,E1,W
# Output: out/q40b/radex_{hcn,hco+}_J4.npy, q40h_radex.log/json
import sys,os,subprocess,itertools,re,json; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
TK=np.array([20,30,50,80,150.]); LN=np.arange(3,7.51,0.25); LNH=np.arange(13,17.51,0.25); LNO=np.arange(12,18.51,0.25); DV=150.
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(a):
    mol,T,ln,lN=a; wd=f'/tmp/rx40_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n60 400\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True); out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out.append((int(p[1]),fl(p[7]),fl(p[6])))
    out=sorted(out)[:4]; return [[o[1],o[2]] for o in out]
def table(mol,lns):
    path=f'/workspace/out/q40b/radex_{mol}_J4.npy'
    if os.path.exists(path): return np.load(path)
    with ProcessPoolExecutor(8) as ex: res=list(ex.map(run,list(itertools.product([mol],TK,LN,lns)),chunksize=40))
    a=np.array(res).reshape(len(TK),len(LN),len(lns),4,2); np.save(path,a); return a
if __name__=='__main__':
    H=table('hcn',LNH); O=table('hco+',LNO); print(H.shape,O.shape)
