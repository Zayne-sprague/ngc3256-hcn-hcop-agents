# Q25: HCN/HCO+ 1-0 and 2-1 per annulus (10" beam, Tb units) and RADEX translation (LAMDA hcn/hco+, H2, Tbg 2.73 K, dv=150): which (Tk,n,X=N(HCN)/N(HCO+)) reproduce R10,R21, and HCN/HCO+ 2-1/1-0 Tb ratios, inner (4-8") vs outer (8-12").
import sys,os,subprocess,itertools,re,json; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
TK=np.array([10,15,20,30,50,80,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25); DV=150.
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(a):
    mol,T,ln,lN=a; wd=f'/tmp/rx25_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n60 200\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True); out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out.append((int(p[1]),fl(p[7]),fl(p[6])))
    out=sorted(out); return [o[1] for o in out[:2]]
def table(mol):
    path=f'/workspace/out/q25/radex_{mol}.npy'
    if os.path.exists(path): return np.load(path)
    with ProcessPoolExecutor(4) as ex: res=list(ex.map(run,list(itertools.product([mol],TK,LN,LNC)),chunksize=50))
    a=np.array(res).reshape(len(TK),len(LN),len(LNC),2); np.save(path,a); return a
if __name__=='__main__':
    H=table('hcn'); O=table('hco+')   # [T,n,N,J(1-0,2-1)] T_R
    d=json.load(open('/workspace/out/q25/q25_ratio10.json'))
    f=1.0126  # Tb ratio = flux ratio * (nu_HCO+/nu_HCN)^2 for same line J
    def obs(src):
        h10=np.array(d[f'{src}|HCN 1-0']['F']); o10=np.array(d[f'{src}|HCO+ 1-0']['F']); h21=np.array(d['Xb55|HCN 2-1']['F']); o21=np.array(d['Xb55|HCO+ 2-1']['F'])
        return dict(R10=h10/o10*f,R21=h21/o21*f,LH=h21/h10*0.25,LO=o21/o10*0.25)
    lab=['0-4','4-8','6-9','8-12']
    for src in ('Michiyama','Harada'):
        o=obs(src); print(f'== {src} (Tb units): annuli',lab)
        for k,v in o.items(): print(f'  {k}: {np.round(v,3)}')
        # gradient: inner (4-8") vs outer (8-12"): which (T,n) satisfy ladders 2-1/1-0 for both HCN and HCO+ within 15%, and corresponding X from R10
        for ia,name in ((1,'4-8'),(3,'8-12')):
            sols=[]
            for it,T in enumerate(TK):
                for iN,n in enumerate(LN):
                    # ladder ratio per species vs N: need both HCN and HCO+ ladder match
                    lh=H[it,iN,:,1]/H[it,iN,:,0]; lo=O[it,iN,:,1]/O[it,iN,:,0]
                    mH=np.abs(lh/o['LH'][ia]-1)<0.15; mO=np.abs(lo/o['LO'][ia]-1)<0.15
                    if mH.any() and mO.any():
                        for a in np.where(mH)[0]:
                            for b in np.where(mO)[0]:
                                R=H[it,iN,a,0]/O[it,iN,b,0]
                                if abs(R/o['R10'][ia]-1)<0.15: sols.append((T,n,LNC[a]-LNC[b]))
            if sols:
                s=np.array(sols); print(f'  annulus {name}": {len(s)} grid solutions; Tk {sorted(set(s[:,0]))}; log n {s[:,1].min():.2f}..{s[:,1].max():.2f}; logX {s[:,2].min():+.2f}..{s[:,2].max():+.2f}; median logX {np.median(s[:,2]):+.2f}')
            else: print(f'  annulus {name}": no solutions within 15%')
    # sensitivity: fixed Tk=30 K, N(HCO+)=1e14.5: R10,R21,LH,LO vs n
    it=3; ia=list(LNC).index(14.5)
    print('\nSensitivity Tk=30K, logN(HCO+)=14.5, X=1 (logN HCN=14.5): logn -> R10,R21 (Tb), HCN 2-1/1-0, HCO+ 2-1/1-0')
    for iN,n in enumerate(LN):
        if 3.5<=n<=7: print(f' {n:.2f}',np.round([H[it,iN,ia,0]/O[it,iN,ia,0],H[it,iN,ia,1]/O[it,iN,ia,1],H[it,iN,ia,1]/H[it,iN,ia,0],O[it,iN,ia,1]/O[it,iN,ia,0]],3))
