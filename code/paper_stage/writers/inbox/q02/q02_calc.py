import json,numpy as np
from scipy.optimize import brentq
d=json.load(open('/workspace/out/q47/ratios10.json'))
rng=np.random.default_rng(1)
for mol in ['HCN','HCO+']:
  for est in 'TAP':
    n=d[f'{mol}|N|{est}'];s=d[f'{mol}|S|{est}']
    a=rng.normal(n['R'],n['err'],200000);b=rng.normal(s['R'],s['err'],200000)
    ok=b>0; dr=a[ok]/b[ok]
    print(mol,est,'N %.4f+-%.4f S %.4f+-%.4f'%(n['R'],n['err'],s['R'],s['err']),'double ratio nominal %.2f'%(n['R']/s['R']),'MC median %.2f 16-84 %.1f-%.1f 2.5%%=%.2f  P(<=1)=%.3f  (S<=0: %.2f)'%((np.median(dr),)+tuple(np.percentile(dr,[16,84,2.5]))+((dr<=1).sum()/len(a),(~ok).mean())))
print()
for R in (0.022,0.028):
  print('R',R,'thin-both X=1/R=%.0f (lower limit if 12C thick)'%(1/R))
  for X in (40,60,66,90,135):
    f=lambda t:t/(1-np.exp(-t))-X*R
    t=brentq(f,1e-6,200) if X*R>1 else 0
    print('  X=%d: tau12=%.2f tau13=%.3f'%(X,t,t/X))
