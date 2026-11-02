# Q4 round 2, fix 2: brightness-temperature (T_b) ratio as the primary opacity discriminator.
# Observed (s43, Xb4b native beam 0.245x0.237", r<0.25" peak Gaussian+baseline, RJ, observed frame; flux-scale systematic ~10% common, 5% on a ratio):
#   N core: T_b(HCN 3-2)=1.461 K, T_b(HCO+ 3-2)=1.556 K (ratio 0.94);  S core: 0.792 K, 1.869 K (ratio 0.42);  width ratios 1.08(N, tentative) / 1.00 (S).
# Model: one RADEX LVG zone per line, shared Tkin, n(H2), dv=150 km/s FWHM, and shared beam-filling eta:  T_b = eta * T_R(RADEX).  HCN and HCO+ columns independent (A=N(HCN)/N(HCO+)).
# Allowed region = grid cells with T_R,HCN/T_R,HCO+ within the observed ratio +-tolerance AND eta=T_b(HCO+)/T_R(HCO+) <= 1 (and eta_HCN = T_b(HCN)/T_R(HCN) <= 1).  Branches by tau(3-2) at line centre.
# Also: the Gaussian-profile opacity-broadening factor f(tau)=FWHM_obs/FWHM_intr to predict the HCN/HCO+ width ratio per cell.
import subprocess,os,json,itertools,numpy as np,pandas as pd
os.makedirs('/workspace/out/radex_q4b',exist_ok=True); os.chdir('/workspace/out/radex_q4b')
LO={'hcn':{3:265.8864,4:354.5055},'hco+':{3:267.5576,4:356.7342}}
def radex(mol,Tk,nH,N,dV=150.):
    inp=f"/opt/lamda/{mol}.dat\nout_{mol}.txt\n170 370\n{Tk}\n1\nH2\n{nH}\n2.73\n{N}\n{dV}\n0\n"
    subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True); o={}
    for l in open(f'out_{mol}.txt'):
        p=l.split()
        if len(p)>=12 and p[1]=='--':
            for J,f in LO[mol].items():
                if abs(float(p[4])-f)<2e-3: o[J]=(float(p[6]),float(p[7]),float(p[8]))   # Tex,tau,T_R
    return o
Tks=[20,30,45,60,90,150]; ns=[1e4,1e5,3e5,1e6,3e6,1e7,1e8]; lN=np.arange(13,19.01,0.25)
cf='/workspace/out/q4r2/s48_grid.csv'
if not os.path.exists(cf):
    rows=[]
    for Tk in Tks:
        for n in ns:
            H={l:radex('hcn',Tk,n,10**l) for l in lN}; O={l:radex('hco+',Tk,n,10**l) for l in lN}
            for lh in lN:
                for lo in lN:
                    A=10**(lh-lo)
                    if A<0.2 or A>30 or 3 not in H[lh] or 3 not in O[lo]: continue
                    rows.append(dict(Tk=Tk,n=n,lNH=lh,lNO=lo,A=A,TexH=H[lh][3][0],tauH=H[lh][3][1],TRH=H[lh][3][2],TexO=O[lo][3][0],tauO=O[lo][3][1],TRO=O[lo][3][2],
                                     TRH4=H[lh][4][2] if 4 in H[lh] else np.nan,TRO4=O[lo][4][2] if 4 in O[lo] else np.nan,tauH4=H[lh][4][1] if 4 in H[lh] else np.nan))
    pd.DataFrame(rows).to_csv(cf,index=False)
df=pd.read_csv(cf); df=df[(df.TRH>0)&(df.TRO>0)].copy()
# opacity broadening factor of a Gaussian profile: FWHM(tau0)/FWHM(thin)
x=np.linspace(0,8,40001)
def fb(t):
    if t<1e-3: return 1.0
    y=1-np.exp(-t*np.exp(-0.5*x**2)); y/=y[0]; i=np.argmax(y<0.5); xh=x[i]; return xh/1.17741
ts=np.logspace(-3,3,300); fv=np.array([fb(t) for t in ts])
f=lambda t:np.interp(np.log10(np.clip(t,1e-3,1e3)),np.log10(ts),fv)
print('opacity broadening factor FWHM/FWHM_thin at tau0 = 0.1,1,3,10,30,100,300:',[round(float(f(t)),3) for t in (.1,1,3,10,30,100,300)])
df['wr_pred']=f(df.tauH)/f(df.tauO); df['dR']=df.TRH/df.TRO
obs={'N core':(1.461,1.556,0.10),'S core':(0.792,1.869,0.10)}
def branch(r):
    if r.tauH<1 and r.tauO<1: return 'both thin'
    if r.tauH>=3 and r.tauO>=3: return 'both thick (tau>=3)'
    if r.tauO>=3 and r.tauH<3: return 'HCO+ thick, HCN thin/mid'
    if r.tauH>=3 and r.tauO<3: return 'HCN thick, HCO+ thin/mid'
    return 'intermediate'
out={}
for nm,(tH,tO,tol) in obs.items():
    ro=tH/tO; m=df[(abs(df.dR/ro-1)<tol+0.05)].copy()   # +-15%: 10% calibration-free scatter + 5% flux-scale
    m['eta']=tO/m.TRO; m['etaH']=tH/m.TRH
    m=m[(m.eta<=1)&(m.etaH<=1.1)&(m.eta>0.003)]; m['br']=m.apply(branch,axis=1)
    print(f'\n===== {nm}: T_b ratio HCN/HCO+ = {ro:.3f} (+-15% window); allowed cells {len(m)} of {len(df)} ({100*len(m)/len(df):.1f}%)')
    g=m.groupby('br')
    summ=g.agg(cells=('A','size'),A_min=('A','min'),A_med=('A','median'),A_max=('A','max'),eta_min=('eta','min'),eta_med=('eta','median'),eta_max=('eta','max'),TexO_min=('TexO','min'),TexO_med=('TexO','median'),TexO_max=('TexO','max'),
              tauH_min=('tauH','min'),tauH_max=('tauH','max'),tauO_min=('tauO','min'),tauO_max=('tauO','max'),wr_min=('wr_pred','min'),wr_med=('wr_pred','median'),wr_max=('wr_pred','max'))
    pd.set_option('display.width',300,'display.max_columns',40); print(summ.round(2).T.to_string())
    # column density N(HCN)/dv and N(HCN) ranges per branch
    m['NH_dv']=10**m.lNH/150.; print('log N(HCN)[cm-2 per 150 km/s] range per branch:'); print(g.lNH.agg(['min','median','max']).round(2).to_string())
    print('Tkin / n of allowed cells per branch (min-max):'); print(g.agg(Tk_min=('Tk','min'),Tk_max=('Tk','max'),n_min=('n','min'),n_max=('n','max')).to_string())
    out[nm]=dict(ratio=ro,cells=len(m),frac=len(m)/len(df),summary=summ.round(3).reset_index().to_dict('records'),lNH=g.lNH.agg(['min','median','max']).round(2).reset_index().to_dict('records'))
    m.to_csv(f'/workspace/out/q4r2/s48_allowed_{nm.split()[0]}.csv',index=False)
# N: how contrived is the thin branch? thin cells need A (and eta) -- list the 'both thin' cells with the R_T ratio inside window as a function of n
# Compare: A needed for the N core thin branch vs A needed for S/outer N (T_b ratio 0.42) with the same grid
for nm in out: print(nm,'allowed fraction of grid',round(out[nm]['frac'],3))
tn=pd.read_csv('/workspace/out/q4r2/s48_allowed_N.csv'); ts_=pd.read_csv('/workspace/out/q4r2/s48_allowed_S.csv')
print('\nThin-branch A: N core median %.2f (%.2f-%.2f) vs S median %.2f (%.2f-%.2f)'%tuple([tn[tn.br=='both thin'].A.median(),tn[tn.br=='both thin'].A.min(),tn[tn.br=='both thin'].A.max(),ts_[ts_.br=='both thin'].A.median(),ts_[ts_.br=='both thin'].A.min(),ts_[ts_.br=='both thin'].A.max()]) if (tn.br=='both thin').any() and (ts_.br=='both thin').any() else ['nan']*6)
# 'minimum filling': for each branch the minimum eta; and Tb(HCN)/Tb(HCO+) = 1 requires tau_H>=? in a single-Tex slab:
print('\nSingle-Tex slab: T_b ratio HCN/HCO+ vs tau(HCO+)=tauO for tauH=r*tauO  (same Tex, Delta v) -> ratio = (1-exp(-tauH))/(1-exp(-tauO)); gives 0.94 for tauH/tauO(thin)=r:')
for tO in (0.1,0.5,1,2,3,5,10,30):
    for r in (1,2,3):
        print(f'  tauO={tO:5.1f} tauH={r*tO:5.1f}: ratio {(1-np.exp(-r*tO))/(1-np.exp(-tO)):.3f}',end=';')
    print()
json.dump(out,open('/workspace/out/q4r2/s48_tb_discriminator.json','w'),indent=1,default=float)
# extra: excitation-temperature contrast per branch (S: what the 0.42 ratio needs when both lines are thick) and fraction of allowed cells per branch
for nm in ('N','S'):
    t=pd.read_csv(f'/workspace/out/q4r2/s48_allowed_{nm}.csv'); t['TexRatio']=t.TexH/t.TexO
    print(f'\n[{nm} core] fraction of allowed cells per branch:',(t.br.value_counts(normalize=True)).round(2).to_dict())
    for b,g in t.groupby('br'): print(f'  {b:28s} TexH/TexO median {g.TexRatio.median():.2f} (16-84%: {g.TexRatio.quantile(.16):.2f}-{g.TexRatio.quantile(.84):.2f}); TexH median {g.TexH.median():.0f} K; A<1 fraction {np.mean(g.A<1):.2f}; A>=3 fraction {np.mean(g.A>=3):.2f}')
