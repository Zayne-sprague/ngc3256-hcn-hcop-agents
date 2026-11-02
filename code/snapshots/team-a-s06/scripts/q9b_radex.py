# Q9b fix 4: RADEX (LAMDA hcn.dat/hco+.dat, H2 collider, Tbg 2.73 K, FWHM 150 km/s; columns scale as N/dv) grids for HCN and HCO+ 3-2 and 4-3,
# beam-filling = 1 for both lines (same region, matched beams). Observables: HCN/HCO+ 3-2 Tb ratio at the N core (0.30"), N wings, S, plus Tb peaks (ceiling check).
# python3 -I q9b_radex.py
import sys,subprocess,os,itertools,re,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
TK=np.array([20,30,50,100,150,300.]); LN=np.arange(3,8.01,0.5); LNC=np.arange(12.5,18.51,0.25); DV=150.
XS=[0.316,1.0,3.16,10.0]   # N(HCN)/N(HCO+)
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(a):
    mol,T,ln,lN=a; wd=f'/tmp/rq9_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n200 400\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True); out={}
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out[int(p[0])]=(fl(p[7]),fl(p[6]),fl(p[5]))   # up -> (T_R, tau, Tex)
    return [out.get(3,(np.nan,)*3),out.get(4,(np.nan,)*3)]
def table(mol):
    path=f'/workspace/out/q9b/radex_{mol}.npz'
    if os.path.exists(path): d=np.load(path); return d['A']
    jobs=list(itertools.product([mol],TK,LN,LNC))
    with ProcessPoolExecutor(12) as ex: r=list(ex.map(run,jobs,chunksize=40))
    A=np.array(r).reshape(len(TK),len(LN),len(LNC),2,3); np.savez(path,A=A); return A
if __name__=='__main__':
    AH=table('hcn'); AO=table('hco+')    # [T,n,N,(3-2,4-3),(TR,tau,Tex)]
    from q9b_lib import TBF
    OBS={'N core 0.30" (r<0.25")':(0.982*TBF,0.013*TBF),'N core 1.05" (r<0.5")':(0.769*TBF,0.004*TBF),
         'N wings 100-260 stacked 0.30" (r<0.25")':(1.051*TBF,0.131*TBF),'N wings 100-260 stacked 1.05" (Xb4d r<0.5")':(1.186*TBF,0.098*TBF),
         'S core 0.30" (r<0.25")':(0.414*TBF,0.014*TBF),'S wings 100-260 stacked 0.30" (r<0.25")':(0.443*TBF,0.045*TBF)}
    rows=[]
    dN=LNC[1]-LNC[0]
    for X in XS:
        sh=int(round(np.log10(X)/dN))   # N(HCO+)=N(HCN)/X -> index shift
        for it,T in enumerate(TK):
            for iN in range(len(LNC)):
                jN=iN-sh
                if jN<0 or jN>=len(LNC): continue
                for i_n,ln in enumerate(LN):
                    h=AH[it,i_n,iN]; o=AO[it,i_n,jN]
                    rows.append(dict(X=X,T=T,logn=ln,logNHCN=LNC[iN],logNHCO=LNC[jN],R32=h[0,0]/o[0,0],R43=h[1,0]/o[1,0],lad_H=h[1,0]/h[0,0],lad_O=o[1,0]/o[0,0],
                                     tauH32=h[0,1],tauO32=o[0,1],TexH32=h[0,2],TexO32=o[0,2],TRH32=h[0,0],TRO32=o[0,0]))
    G=pd.DataFrame(rows); G.to_csv('/workspace/out/q9b/radex_grid_ratios.csv',index=False)
    print(len(G),'grid points')
    # reliability cuts: RADEX is not trustworthy for tau>~1e3 (negative/jumping tau,T_R above logN/dv~17.5): keep logN(HCN,HCO+)<=17.5 (per 150 km/s) and |tau|<1000
    G=G[(G.logNHCN<=17.5)&(G.logNHCO<=17.5)&(G.tauH32<1000)&(G.tauO32<1000)&(G.tauH32>0)&(G.tauO32>0)]
    # native-beam peak Tb (Q4 entry, Xb4b 0.25") : N HCN 2.63 HCO+ 2.82 K; S 1.34 / 2.94 K -> beam-filling eta = Tb/T_R must be <=1, i.e. T_R >= Tb_obs
    TBO={'N':(2.63,2.82),'S':(1.34,2.94)}
    OBSREG={'N core 0.30" (r<0.25")':'N','N core 1.05" (r<0.5")':'N','N wings 100-260 stacked 0.30" (r<0.25")':'N','N wings 100-260 stacked 1.05" (Xb4d r<0.5")':'N','S core 0.30" (r<0.25")':'S','S wings 100-260 stacked 0.30" (r<0.25")':'S'}
    tabs=[]
    for name,(R,e) in OBS.items():
        tH,tO=TBO[OBSREG[name]]
        s=G[(np.abs(G.R32-R)<max(e,0.05*R))&(G.TRH32>=tH)&(G.TRO32>=tO)]
        print(f'\n{name}: Tb ratio {R:.3f}+-{max(e,0.05*R):.3f} (model tolerance: max(err,5%)); allowed grid points (T_R>=native Tb, tau<1000, logN<=17.5) {len(s)}')
        for X,g in s.groupby('X'):
            thick=((g.tauH32>1)&(g.tauO32>1)).mean(); etaO=tO/g.TRO32
            print(f'  X={X}: n={len(g)}; log n {g.logn.min()}-{g.logn.max()}; T {g["T"].min():.0f}-{g["T"].max():.0f} K; tau(HCN32) {g.tauH32.min():.2f}-{g.tauH32.max():.1f} (med {g.tauH32.median():.1f}); tau(HCO+32) {g.tauO32.min():.2f}-{g.tauO32.max():.1f} (med {g.tauO32.median():.1f}); frac both tau>1 {thick:.2f}; filling eta=Tb/T_R(HCO+) {etaO.min():.3f}-{etaO.max():.2f}; ladder HCN {g.lad_H.median():.2f} HCO+ {g.lad_O.median():.2f}')
            tabs.append(dict(obs=name,X=X,n=len(g),logn_lo=g.logn.min(),logn_hi=g.logn.max(),T_lo=g['T'].min(),T_hi=g['T'].max(),tauH_med=g.tauH32.median(),tauO_med=g.tauO32.median(),frac_thick=thick,eta_lo=etaO.min(),eta_hi=etaO.max()))
    pd.DataFrame(tabs).to_csv('/workspace/out/q9b/radex_allowed.csv',index=False)
    # thick+thermal limit: tau>>1 and Tex equal -> ratio -> 1
    thick=G[(G.tauH32>5)&(G.tauO32>5)]; print('\nboth tau>5: R32 range %.2f-%.2f (median %.2f); with Tex equal within 10%%: %.2f-%.2f'%(thick.R32.min(),thick.R32.max(),thick.R32.median(),thick[np.abs(thick.TexH32/thick.TexO32-1)<.1].R32.min(),thick[np.abs(thick.TexH32/thick.TexO32-1)<.1].R32.max()))
    low=G[(G.R32<0.5)&(G.R32>0.3)]; print('R32 0.3-0.5 (S-like): fraction with both tau>3: %.2f; HCN tau>3 & HCO+ tau>3 & |Tex diff|<10%%: %d of %d'%(((low.tauH32>3)&(low.tauO32>3)).mean(),((low.tauH32>3)&(low.tauO32>3)&(np.abs(low.TexH32/low.TexO32-1)<.1)).sum(),len(low)))
    print(low.groupby('X').apply(lambda g:pd.Series(dict(n=len(g),frac_both_thick=((g.tauH32>3)&(g.tauO32>3)).mean(),tauH_med=g.tauH32.median(),tauO_med=g.tauO32.median()))).round(2).to_string())
    hi=G[(G.R32>0.9)&(G.R32<1.2)]; print('R32 0.9-1.2 (N-like):'); print(hi.groupby('X').apply(lambda g:pd.Series(dict(n=len(g),frac_both_thick=((g.tauH32>3)&(g.tauO32>3)).mean(),tauH_med=g.tauH32.median(),tauO_med=g.tauO32.median()))).round(2).to_string())
    # ratio vs column at fixed n,T to show how R32 depends on N (opacity) per X : e.g. T=100 K, n=1e6
    print('\nR32 vs logN(HCN) at T=100 K, log n=6 for each X:')
    sub=G[(G['T']==100)&(G.logn==6)]; print(sub.pivot(index='logNHCN',columns='X',values='R32').round(2).to_string())
    print('\nR32 vs logN(HCN) at T=100 K, log n=5:'); sub=G[(G['T']==100)&(G.logn==5)]; print(sub.pivot(index='logNHCN',columns='X',values='R32').round(2).to_string())
