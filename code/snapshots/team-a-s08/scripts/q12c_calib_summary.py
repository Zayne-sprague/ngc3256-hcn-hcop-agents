# Q12c calibration summary: per truth/setting, distribution of dlnZ(2cT-1z) and of the one-zone profile chi2; false-positive rate P(dlnZ>=3 | truth one-zone [or mis-specified one-zone]) with Clopper-Pearson 95% bounds,
# power P(dlnZ>=3 | truth two-comp), and p-value of the observed N data (dlnZ_obs, chi2_obs from q12c_summary lnZ table) under each truth. Output out/q12c/calib_summary.csv
import pandas as pd, numpy as np, glob, os
from scipy.stats import beta
D='/workspace/out/q12c/'; rows=[]
fs=glob.glob(D+'runs/cal_*_sim*.csv'); R=pd.concat([pd.read_csv(f) for f in fs],ignore_index=True)
R['sim']=[int(os.path.basename(f).split('_sim')[1].split('_')[0]) for f in fs]
obs=pd.read_csv(D+'lnZ_table.csv').set_index(['region','tag']); 
def cp(k,n): return (beta.ppf(.025,k,n-k+1) if k>0 else 0.,beta.ppf(.975,k+1,n-k) if k<n else 1.)
for tag,gg in R.groupby('tag'):
    P=gg.pivot_table(index='sim',columns='model',values='logZ'); C=gg[gg.model=='1z'].set_index('sim').chi2_prof
    P=P.dropna(); d=(P['2cT']-P['1z']); c=C.loc[P.index]
    ts,st=tag.split('_')[1],tag.split('_')[2]; n=len(d); k3=int((d>=3).sum()); k1=int((d>=1).sum())
    ob=obs.loc[('N','base' if st=='base' else 'geo')]; dobs=ob['2cT']-ob['1z']
    lo,hi=cp(k3,n)
    rows.append(dict(truth=ts,fit=st,nsim=n,dlnZ_med=d.median(),dlnZ_16=d.quantile(.16),dlnZ_84=d.quantile(.84),P_dlnZ_ge3=k3/n,P_ge3_lo95=lo,P_ge3_hi95=hi,P_dlnZ_ge1=k1/n,dlnZ_obs=dobs,P_dlnZ_ge_obs=(d>=dobs).mean(),
        chi2_1z_med=c.median(),P_chi2_ge_9p1=(c>=9.13).mean(),P_chi2_ge_3p2=(c>=3.21).mean()))
S=pd.DataFrame(rows).sort_values(['truth','fit']); S.round(3).to_csv(D+'calib_summary.csv',index=False); pd.set_option('display.width',250,'display.max_columns',30); print(S.round(3).to_string())
