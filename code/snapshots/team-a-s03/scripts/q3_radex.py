# Q3 step 6: RADEX (LAMDA hcn.dat/hco+.dat, H2, Tbg 2.73 K, dv 150 km/s) for the Soutflow blue lobe: constraints from 3-2/1-0 Tb ratios (Xb4d 3-2 at 2.1" vs Michiyama+2018 1-0)
# and HCN/HCO+ 3-2. Output: allowed n(H2) per Tkin for each species, and joint (n,T,X=N(HCN)/N(HCO+)) region.
import sys,subprocess,os,itertools,re,warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
TK=np.array([15,30,50,100,150.]); LN=np.arange(3,7.51,0.25); LNC=np.arange(12.5,17.51,0.25); DV=150.
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(a):
    mol,T,ln,lN=a; wd=f'/tmp/rq_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n80 280\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True); out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out.append((int(p[1]),fl(p[7]),fl(p[6])))
    out=sorted(out); return [(o[1],o[2]) for o in out[:3]]    # 1-0, 2-1, 3-2: T_R, tau
def table(mol):
    path=f'/workspace/out/q3/radex_{mol}.npz'
    if os.path.exists(path): d=np.load(path); return d['TR'],d['tau']
    jobs=list(itertools.product([mol],TK,LN,LNC))
    with ProcessPoolExecutor(4) as ex: r=list(ex.map(run,jobs,chunksize=40))
    A=np.array(r).reshape(len(TK),len(LN),len(LNC),3,2); np.savez(path,TR=A[...,0],tau=A[...,1]); return A[...,0],A[...,1]
if __name__=='__main__':
    TRH,TH=table('hcn'); TRO,TO=table('hco+')
    r31H=TRH[...,2]/TRH[...,0]; r31O=TRO[...,2]/TRO[...,0]
    # observables (Vsys=2775 case of q3_r31_vs_M18.csv; error widened to span Vsys 2740-2810): 
    O=pd.read_csv('/workspace/out/q3/q3_r31_vs_M18.csv'); b=O[O.pos=='M18 South-OF-Blue']
    def ob(sp): x=b[b.line==sp]; m=x[x.vsys==2775].iloc[0]; lo,hi=x.r31.min(),x.r31.max(); return m.r31,np.hypot(m.e_r31,(hi-lo)/2)
    (oH,eH),(oO,eO)=ob('HCN'),ob('HCO+'); R32,eR=1.7,0.4   # HCN/HCO+ 3-2 Tb at Soutflow blue wing (q3_velocity_ratio / q1_soutflow_wing_variants: 1.3-1.8)
    print(f'observables: r31(HCN)={oH:.2f}+-{eH:.2f}, r31(HCO+)={oO:.2f}+-{eO:.2f}, HCN/HCO+ 3-2={R32}+-{eR}')
    rows=[]
    for it,T in enumerate(TK):
        for sp,r31,o,e in (('HCN',r31H,oH,eH),('HCO+',r31O,oO,eO)):
            # for each n: range of column over which r31 within 1 sigma (min over N), report n range for which any N works
            ok=np.abs(r31[it]-o)<e   # (n,N)
            ns=LN[ok.any(1)]
            rows.append(dict(T=T,line=sp,logn_lo=ns.min() if len(ns) else np.nan,logn_hi=ns.max() if len(ns) else np.nan))
    D=pd.DataFrame(rows); print(D.to_string()); D.to_csv('/workspace/out/q3/q3_radex_r31_ranges.csv',index=False)
    # joint: chi2 over (T,n,NH,NO): r31H, r31O, R32 (3-2 Tb ratio)
    R32m=TRH[:,:,:,None,2]/TRO[:,:,None,:,2]
    chi=((r31H[:,:,:,None]-oH)/eH)**2+((r31O[:,:,None,:]-oO)/eO)**2+((R32m-R32)/eR)**2
    chi=np.where(np.isfinite(chi),chi,1e9); i=np.unravel_index(np.argmin(chi),chi.shape)
    print('best joint: T=%g logn=%.2f logNH=%.2f logNO=%.2f chi2=%.2f (3 obs, 4 params -> unconstrained)'%(TK[i[0]],LN[i[1]],LNC[i[2]],LNC[i[3]],chi[i]))
    rows=[]
    for it,T in enumerate(TK):
        m=chi[it]<=chi.min()+1.0
        if m.any():
            nn=LN[m.any((1,2))]; xx=(LNC[:,None]-LNC[None,:])[None] ; X=np.where(m,np.broadcast_to(LNC[None,:,None]-LNC[None,None,:],m.shape),np.nan)
            rows.append(dict(T=T,chi2min=chi[it].min(),logn_lo=nn.min(),logn_hi=nn.max(),logX_lo=np.nanmin(X),logX_hi=np.nanmax(X)))
        else: rows.append(dict(T=T,chi2min=chi[it].min()))
    J=pd.DataFrame(rows);print(J.round(2).to_string());J.round(3).to_csv('/workspace/out/q3/q3_radex_joint.csv',index=False)
