"""r10(a): HCN/HCO+ 1-0 N/S from Harada+2018 (ApJ 855, 49; 2018ApJ...855...49H) Table 13, 1.7" cubes, +-150 km/s about 2775 km/s.
Inputs are TABULATED values transcribed from the paper (python3 /task/tools/lit.py read 2018ApJ...855...49H --grep "HCN 1-0 88.634"); nothing is re-measured here.
Columns: I dV (K km/s), err; Ipeak (mK). Errors propagated in quadrature."""
import numpy as np
T={'N':dict(HCN=(197.23,0.96,1111.22),HCO=(205.25,0.95,1173.21)),
   'S':dict(HCN=(56.68,0.96,510.57),HCO=(87.87,0.95,773.12))}
R={}
for k in 'NS':
    h,o=T[k]['HCN'],T[k]['HCO']; r=h[0]/o[0]; e=r*np.hypot(h[1]/h[0],o[1]/o[0]); rp=h[2]/o[2]
    R[k]=(r,e,rp); print(f'{k}: HCN/HCO+ integrated = {r:.3f} +- {e:.3f}; peak-intensity ratio = {rp:.3f}')
ns=R['N'][0]/R['S'][0]; ens=ns*np.hypot(R['N'][1]/R['N'][0],R['S'][1]/R['S'][0])
print(f'N/S integrated = {ns:.3f} +- {ens:.3f}; from 2-decimal ratios 0.96/0.65 = {0.96/0.65:.3f}; N/S peak = {R["N"][2]/R["S"][2]:.3f}')
# our own matched-beam measurement (2016.1.00965 X94a spw27, out/q55/r10_core.csv, scripts/q55f_r10_core.py), beams 1.4 and 1.9 bracket 1.7"
import csv
d={(float(r['beam']),r['reg']):(float(r['R10_Tb']),float(r['err'])) for r in csv.DictReader(open('/workspace/out/q55/r10_core.csv'))}
for b in (1.0,1.4,1.9,2.5):
    n,s=d[(b,'N')],d[(b,'S')]; x=n[0]/s[0]; print(f'team-a beam {b}": N {n[0]:.3f} S {s[0]:.3f} N/S {x:.3f} +- {x*np.hypot(n[1]/n[0],s[1]/s[0]):.3f}')
