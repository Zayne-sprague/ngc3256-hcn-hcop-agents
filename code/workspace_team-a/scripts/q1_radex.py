# Q1 fix 4: quick RADEX grid (LAMDA hcn.dat, hco+.dat, H2 partner, Tbg 2.73 K, uniform sphere, dv=150 km/s) for the deblended N and S ladders.
# Observables per region (Tb, component-integrated, same beam filling assumed for all lines): HCN & HCO+ 3-2[Xb4f]/2-1, 4-3/3-2[Xb4f], and HCN/HCO+ 3-2 (1" templates).
# Parameters: Tkin, n(H2), N(HCN)/dv, N(HCO+)/dv (common Tkin, n).  chi2 with 14% cross-band calibration + the ensemble spread of the deblending in quadrature.
# Tests: best fit / allowed ranges per region; X=N(HCN)/N(HCO+) required; equal-X test (N and S forced to the same X, free n,T,N); optical-depth consistency with HCN/H13CN limits.
import sys,subprocess,os,warnings,itertools; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
TK=np.array([10,15,20,30,50,80,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25); DV=150.
import re
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(args):
    mol,T,ln,lN=args; wd=f'/tmp/rx_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n150 400\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True)
    out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit():
            out.append((int(p[1]),fl(p[5]),fl(p[6]),fl(p[7])))      # lower J, Tex, tau, T_R
    out=sorted(out,key=lambda x:x[0]); return [(o[3],o[2]) for o in out[:3]]   # T_R,tau for 2-1,3-2,4-3
def table(mol):
    path=f'/workspace/out/q1/radex_{mol}.npz'
    if os.path.exists(path): d=np.load(path); return d['TR'],d['tau']
    jobs=list(itertools.product([mol],TK,LN,LNC))
    with ProcessPoolExecutor(4) as ex: res=list(ex.map(run,jobs,chunksize=50))
    arr=np.array(res).reshape(len(TK),len(LN),len(LNC),3,2); np.savez(path,TR=arr[...,0],tau=arr[...,1]); return arr[...,0],arr[...,1]
TRH,TAUH=table('hcn'); TRO,TAUO=table('hco+')
# ---- observables from the deblending ensemble ----
L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv'); E=pd.read_csv('/workspace/out/q1/q1_deblend_summary.csv')
def lad(comp,sp,q):
    r=L[(L.comp==comp)&(L.species==sp)&(L.quantity==q)].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
def r32(comp):
    r=E[(E.model=='nomissing')&(E.comp==comp)&(E.quantity.str.startswith('HCN/HCO+ 3-2 (Tb) from 1'))].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
XB=0.14
OBS={}
for c in ('N','S'):
    o=[]
    for sp in ('HCN','HCO+'):
        for q in ('3-2 [Xb4f]/2-1 (Tb)','4-3/3-2 [Xb4f 3-2] (Tb)'):
            v,s=lad(c,sp,q); o.append((v,np.hypot(XB*v,s)))
    v,s=r32(c); o.append((v,np.hypot(0.05*v,s))); OBS[c]=np.array(o)
    print(c,'observables [HCN 3-2/2-1, HCN 4-3/3-2, HCO+ 3-2/2-1, HCO+ 4-3/3-2, HCN/HCO+ 3-2]:',np.round(o,3).tolist())
# ---- model predictions on the grid: axes (T,n,NH,NO) ----
def pred():
    h=TRH; o=TRO
    hA=h[:,:,:,1]/h[:,:,:,0]; hB=h[:,:,:,2]/h[:,:,:,1]; oA=o[:,:,:,1]/o[:,:,:,0]; oB=o[:,:,:,2]/o[:,:,:,1]
    P=np.empty((len(TK),len(LN),len(LNC),len(LNC),5))
    P[...,0]=hA[:,:,:,None]; P[...,1]=hB[:,:,:,None]; P[...,2]=oA[:,:,None,:]; P[...,3]=oB[:,:,None,:]
    P[...,4]=h[:,:,:,None,1]/o[:,:,None,:,1]
    return P
P=pred(); P=np.where(np.isfinite(P),P,1e9)
logX=(LNC[:,None]-LNC[None,:])      # log10 N_HCN/N_HCO+
res={}
rows=[]; EX=[]
CUTS=[('no tau cut',1e9),('tau43(HCN)<=5.5 [H13CN>9.5, A=50]',5.5),('tau43(HCN)<=28.5 [H13CN>2.3, A=50]',28.5)]
tau43=TAUH[:,:,:,2][:,:,:,None]*np.ones((1,1,1,len(LNC)))
for cutname,cut in CUTS:
    res={}
    for c,ob in OBS.items():
        chi=(((P-ob[:,0])/ob[:,1])**2).sum(-1); chi=np.where(tau43<=cut,chi,1e9); res[c]=chi; i=np.unravel_index(np.argmin(chi),chi.shape); m=chi.min()
        for lab,dc in (('1sig(4par)',4.72),('2sig',9.70)):
            sel=chi<=m+dc; iT,iN_,iH,iO=np.where(sel)
            tauH=TAUH[iT,iN_,iH,1]; tauO=TAUO[iT,iN_,iO,1]
            rows.append(dict(cut=cutname,region=c,level=lab,chi2min=m,best_Tkin=TK[i[0]],best_logn=LN[i[1]],best_logNH_dv=LNC[i[2]]-np.log10(DV),best_logNO_dv=LNC[i[3]]-np.log10(DV),
                best_logX=LNC[i[2]]-LNC[i[3]],best_tau32_HCN=TAUH[i[0],i[1],i[2],1],best_tau32_HCO=TAUO[i[0],i[1],i[3],1],
                Tkin_min=TK[iT].min(),Tkin_max=TK[iT].max(),logn_min=LN[iN_].min(),logn_max=LN[iN_].max(),logX_min=(LNC[iH]-LNC[iO]).min(),logX_max=(LNC[iH]-LNC[iO]).max(),
                tau32_HCN_min=tauH.min(),tau32_HCN_max=tauH.max(),tau32_HCO_max=tauO.max(),npts=int(sel.sum())))
        print(cutname,c,'best',dict(T=TK[i[0]],logn=LN[i[1]],logNH=LNC[i[2]],logNO=LNC[i[3]]),'chi2 %.2f'%m,'pred',np.round(P[i],3))
    def minchi_vs_X(chi):
        out={}
        for x in np.round(np.arange(-3,3.01,0.25),2):
            mask=np.isclose(np.round(logX,2),x)
            if mask.any(): ih,io=np.where(mask); out[x]=chi[:,:,ih,io].min()
        return out
    cN,cS=minchi_vs_X(res['N']),minchi_vs_X(res['S']); xs=sorted(set(cN)&set(cS)); tot={x:cN[x]+cS[x] for x in xs}; xb=min(tot,key=tot.get); free=res['N'].min()+res['S'].min()
    EX.append(dict(cut=cutname,common_logX=xb,chi2_equalX=tot[xb],chi2_free=free,delta_chi2=tot[xb]-free,chiN_min_at_X=cN[xb],chiS_min_at_X=cS[xb]))
    print(cutname,'equal-X: logX',xb,'chi2',round(tot[xb],2),'free',round(free,2),'delta',round(tot[xb]-free,2))
R=pd.DataFrame(rows); R.to_csv('/workspace/out/q1/q1_radex_fits.csv',index=False); pd.set_option('display.width',250,'display.max_columns',50); print(R.round(2).to_string())
pd.DataFrame(EX).to_csv('/workspace/out/q1/q1_radex_equalX_test.csv',index=False)
# what R32 contrast can equal-X models produce while matching the ladders (chi2_ladders only, 4 ladder data): allowed sets per region, compare predicted R32
lad_chi={c:(((P[...,:4]-OBS[c][:4,0])/OBS[c][:4,1])**2).sum(-1) for c in OBS}
rr=[]
for x in (-1.0,-0.5,-0.25,0.0,0.25,0.5,1.0):
    mask=np.isclose(np.round(logX,2),x); ih,io=np.where(mask)
    if not mask.any(): continue
    q={}
    for c in OBS:
        ch=lad_chi[c][:,:,ih,io]; j=np.unravel_index(np.argmin(ch),ch.shape); q[c]=(ch.min(),P[j[0],j[1],ih[j[2]],io[j[2]],4])
        # range of R32 among ladder-acceptable models (delta chi2<4 on 4 data)
        ok=ch<=ch.min()+4; ii=np.where(ok); r32=P[ii[0],ii[1],ih[ii[2]],io[ii[2]],4]; q[c]+=(r32.min(),r32.max())
    rr.append(dict(logX=x,chi2_ladders_N=q['N'][0],R32_pred_N=q['N'][1],R32_range_N=f"{q['N'][2]:.2f}-{q['N'][3]:.2f}",chi2_ladders_S=q['S'][0],R32_pred_S=q['S'][1],R32_range_S=f"{q['S'][2]:.2f}-{q['S'][3]:.2f}"))
RR=pd.DataFrame(rr); RR.to_csv('/workspace/out/q1/q1_radex_equalX.csv',index=False); print(RR.round(3).to_string())
# ---- H13CN opacity bound on tau(HCN 4-3): intensity ratio HCN/H13CN > Rlim  <=> (1-e^-t)/(1-e^-t/A) > Rlim  (A = 12C/13C abundance ratio, assumed 30,50,90)
tb=[]
for A in (30,50,90):
    for Rl in (2.3,9.5):
        t=np.logspace(-2,2,4000); rr_=(1-np.exp(-t))/(1-np.exp(-t/A)); ok=t[rr_>Rl]
        tb.append(dict(A_12C13C=A,HCN_over_H13CN_lower=Rl,tau_HCN43_max=ok.max() if len(ok) else 0.))
TB=pd.DataFrame(tb); TB.to_csv('/workspace/out/q1/q1_h13cn_tau_bound.csv',index=False); print(TB.round(2).to_string())
