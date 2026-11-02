# Q4a-2: collapse the wing-model fits of q4_h13cn.py into per-region HCN/H13CN 4-3 lower limits (Tb ratio = flux ratio*(nu_H13CN/nu_HCN)^2=*0.949) and the optical-depth bound for A=12C/13C=30,50,90.
# "adequate" = chi2/(ndata-npar-1) <= 1.5 ; conservative = smallest 3-sigma lower limit among adequate models (both fit ranges); best-model = the lowest-chi2 model.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
OUT='/workspace/out/q4/'
d=pd.read_csv(OUT+'h13cn_common_beam.csv'); d['chi2red']=d.chi2/(d.ndata-d.npar-1); K=(345.339760/354.505478)**2
rows=[]; tt=[]
def taumax(R,A):
    t=np.logspace(-2,2,6000); r=(1-np.exp(-t))/(1-np.exp(-t/A)); ok=t[r>R]; return ok.max() if len(ok) else 0.
def tau_of(R,A):   # tau(HCN) giving observed Tb ratio R (decreasing in tau); nan if R>=A
    if R>=A: return 0.
    t=np.logspace(-2,2,6000); r=(1-np.exp(-t))/(1-np.exp(-t/A)); return t[np.argmin(abs(r-R))]
for reg,g in d.groupby('region',sort=False):
    ad=g[(g.chi2red<=1.5)&g.alpha_3s_upper.notna()&(g.alpha_3s_upper>0)]; 
    if ad.empty: continue
    Rcons=(ad.HCN_over_H13CN_flux_lower.min())*K; b=ad.loc[ad.chi2.idxmin()]; Rb=b.HCN_over_H13CN_flux_lower*K
    # range of best-fit alphas among adequate models (flux ratio HCN/H13CN)
    a_lo,a_hi=ad.alpha_best.min(),ad.alpha_best.max()
    rows.append(dict(region=reg,n_adequate=len(ad),models=';'.join(sorted(set(ad.model+'@'+ad.fit_GHz))),R_Tb_lower_conservative=Rcons,R_Tb_lower_bestmodel=Rb,best_model=b.model+'@'+b.fit_GHz,
                     alpha_best_min=a_lo,alpha_best_max=a_hi,R_flux_if_detected_max=1/a_lo if a_lo>0 else np.nan,R_flux_if_detected_min=1/a_hi if a_hi>0 else np.nan))
    for A in (30,50,90):
        tt.append(dict(region=reg,A_12C13C=A,tau43_max_conservative=taumax(Rcons,A),tau43_max_bestmodel=taumax(Rb,A),
                       tau43_if_H13CN_real_at_max_alpha=tau_of(1/a_hi*K,A) if a_hi>0 else np.nan,tau43_if_H13CN_real_at_min_alpha=tau_of(1/a_lo*K,A) if a_lo>0 else np.nan))
R=pd.DataFrame(rows); R.round(3).to_csv(OUT+'h13cn_Rlimits.csv',index=False); T=pd.DataFrame(tt); T.round(2).to_csv(OUT+'h13cn_tau.csv',index=False)
pd.set_option('display.width',250); print(R.round(3).to_string()); print(T.round(2).to_string())
