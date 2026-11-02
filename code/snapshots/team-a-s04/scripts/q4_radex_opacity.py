# Q4c: direct RADEX opacity test. H13CN (LAMDA h13cn@xpol.dat) is run on the same (Tkin, n, N) grid as the Q1 HCN/HCO+ ladder fits with N(H13CN)=N(HCN)/A (A=12C/13C = 30, 50, 90 assumed);
# the model Tb ratio HCN(4-3)/H13CN(4-3) (same file, same dv=150, Tbg 2.73 K) must exceed the 3-sigma lower limit measured in q4_h13cn.py (Tb ratio = flux ratio*0.949), replacing the tau<=cut shortcut of q1_radex.py.
# Reports, for N and S (the two regions with Q1 ladder observables), the allowed Tkin / n / N(HCN) / tau(3-2), tau(4-3) of HCN and HCO+ with and without the H13CN condition, delta-chi2 contours of the 4-parameter Q1 fit.
import sys,os,re,itertools,subprocess,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
TK=np.array([10,15,20,30,50,80,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25); DV=150.
OUT='/workspace/out/q4/'
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(args):
    T,ln,lN=args; wd=f'/tmp/rx4_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/h13cn@xpol.dat\n{wd}/o.txt\n150 400\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True)
    out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out.append((int(p[1]),fl(p[6]),fl(p[7])))
    out=sorted(out,key=lambda x:x[0]); return [(o[2],o[1]) for o in out[:3]]
def grid(shift):   # log10 N(H13CN) = LNC - shift
    path=OUT+f'radex_h13cn_shift{shift:.2f}.npz'
    if os.path.exists(path): d=np.load(path); return d['TR'],d['tau']
    jobs=[(T,ln,lN-shift) for T in TK for ln in LN for lN in LNC]
    with ProcessPoolExecutor(4) as ex: res=list(ex.map(run,jobs,chunksize=60))
    a=np.array(res).reshape(len(TK),len(LN),len(LNC),3,2); np.savez(path,TR=a[...,0],tau=a[...,1]); return a[...,0],a[...,1]
def main():
    # --- observables from the Q1 deblended ladders (same as q1_radex.py)
    L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv'); E=pd.read_csv('/workspace/out/q1/q1_deblend_summary.csv')
    def lad(c,sp,q): r=L[(L.comp==c)&(L.species==sp)&(L.quantity==q)].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
    def r32(c): r=E[(E.model=='nomissing')&(E.comp==c)&(E.quantity.str.startswith('HCN/HCO+ 3-2 (Tb) from 1'))].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
    OBS={}
    for c in ('N','S'):
        o=[]
        for sp in ('HCN','HCO+'):
            for q in ('3-2 [Xb4f]/2-1 (Tb)','4-3/3-2 [Xb4f 3-2] (Tb)'):
                v,s=lad(c,sp,q); o.append((v,np.hypot(0.14*v,s)))
        v,s=r32(c); o.append((v,np.hypot(0.05*v,s))); OBS[c]=np.array(o)
    d1=np.load('/workspace/out/q1/radex_hcn.npz'); TRH,TAUH=d1['TR'],d1['tau']; d2=np.load('/workspace/out/q1/radex_hco+.npz'); TRO,TAUO=d2['TR'],d2['tau']
    h=TRH; o=TRO
    P=np.empty((len(TK),len(LN),len(LNC),len(LNC),5))
    P[...,0]=(h[...,1]/h[...,0])[:,:,:,None]; P[...,1]=(h[...,2]/h[...,1])[:,:,:,None]; P[...,2]=(o[...,1]/o[...,0])[:,:,None,:]; P[...,3]=(o[...,2]/o[...,1])[:,:,None,:]; P[...,4]=h[:,:,:,None,1]/o[:,:,None,:,1]
    P=np.where(np.isfinite(P),P,1e9)
    # H13CN-HCN model Tb ratio at 4-3 (same molecular file): R(T,n,N_HCN;A)
    TR0,_=grid(0.0)
    lim=pd.read_csv(OUT+'h13cn_tau_limits.csv') if os.path.exists(OUT+'h13cn_tau_limits.csv') else None
    rows=[]
    LIMS={'N':{'conservative (min over adequate wing models)':None},'S':{}}
    RL=pd.read_csv(OUT+'h13cn_Rlimits.csv').set_index('region')
    for A in (30,50,90):
        TRA,_=grid(np.log10(A)); R43=TR0[...,2]/TRA[...,2]      # (T,n,N): HCN/H13CN Tb ratio of the model
        for c in ('N','S'):
            chi0=(((P-OBS[c][:,0])/OBS[c][:,1])**2).sum(-1)
            for lname,Rl in (('no H13CN condition',0.),('H13CN conservative limit',RL.loc[c,'R_Tb_lower_conservative']),('H13CN best-model limit',RL.loc[c,'R_Tb_lower_bestmodel'])):
                ok=(R43>=Rl)[:,:,:,None]*np.ones((1,1,1,len(LNC)),bool)
                chi=np.where(ok,chi0,1e9); m=chi.min(); i=np.unravel_index(np.argmin(chi),chi.shape)
                for lab,dc in (('1sig(4par)',4.72),('2sig',9.70)):
                    sel=chi<=m+dc; iT,iN_,iH,iO=np.where(sel)
                    t32=TAUH[iT,iN_,iH,1]; t43=TAUH[iT,iN_,iH,2]; o32=TAUO[iT,iN_,iO,1]; o43=TAUO[iT,iN_,iO,2]
                    rows.append(dict(A_12C13C=A,region=c,condition=lname,R_Tb_lower=Rl,level=lab,chi2min=m,best_Tkin=TK[i[0]],best_logn=LN[i[1]],best_logN_HCN_cm2=LNC[i[2]],best_logX=LNC[i[2]]-LNC[i[3]],
                       best_tau32_HCN=TAUH[i[0],i[1],i[2],1],best_tau43_HCN=TAUH[i[0],i[1],i[2],2],best_R43_model=R43[i[0],i[1],i[2]],
                       Tkin_min=TK[iT].min(),Tkin_max=TK[iT].max(),logn_min=LN[iN_].min(),logn_max=LN[iN_].max(),logNHCN_min=LNC[iH].min(),logNHCN_max=LNC[iH].max(),
                       tau32_HCN_min=t32.min(),tau32_HCN_max=t32.max(),tau43_HCN_min=t43.min(),tau43_HCN_max=t43.max(),R43_model_min=R43[iT,iN_,iH].min(),R43_model_max=R43[iT,iN_,iH].max(),tau32_HCO_max=o32.max(),tau43_HCO_max=o43.max(),npts=int(sel.sum())))
    R=pd.DataFrame(rows); R.to_csv(OUT+'radex_opacity.csv',index=False)
    pd.set_option('display.width',300,'display.max_columns',60); print(R[R.level=='1sig(4par)'].round(2).drop(columns=['level']).to_string())
if __name__=='__main__': main()
