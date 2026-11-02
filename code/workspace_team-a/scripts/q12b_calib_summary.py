# Calibration of chi2_prof and Delta lnZ (2cT - 1z) on simulated N ratio data (truth = best-fit 1z or 2cT; 40 sims each) -> out/q12b/calib.csv, p-values of the observed values.
import pandas as pd, numpy as np, glob
R=pd.concat([pd.read_csv(f).assign(f=f) for f in glob.glob('/workspace/out/q12b/runs/cal_*.csv')])
R['truth']=R.f.str.extract(r'truth_(\w+)\.csv')[0]; R['sim']=R.f.str.extract(r'_sim(\d+)_')[0].astype(int)
P=R.pivot_table(index=['truth','sim'],columns='model',values=['chi2_prof','logZ'])
P['dlnZ']=P[('logZ','2cT')]-P[('logZ','1z')]; P.columns=['_'.join(c).strip('_') for c in P.columns]
P.to_csv('/workspace/out/q12b/calib.csv')
obs_chi,obs_dz=9.24,3.25
for t in ('1z','2cT'):
    x=P.loc[t]; print(t,'n',len(x),'chi2_prof(1z fit) median %.2f 16-84 [%.2f,%.2f]'%tuple(np.percentile(x.chi2_prof_1z,[50,16,84])),'frac>=9.24: %.3f'%(x.chi2_prof_1z>=obs_chi).mean(),
      '| dlnZ median %.2f [%.2f,%.2f]'%tuple(np.percentile(x.dlnZ,[50,16,84])),'frac>=3.25: %.3f'%(x.dlnZ>=obs_dz).mean(),'frac>=3: %.3f'%(x.dlnZ>=3).mean())
