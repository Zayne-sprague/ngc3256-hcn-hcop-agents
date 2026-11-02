"""Q55f: HCN/HCO+ 1-0 at sub-arcsec to 2.5" from ONE cube (2016.1.00965 X94a spw27, 0.76x0.63"; both lines in the same spw, same uv coverage)
vs 2015.1.00993 (2.1x1.7", complete); fixed window +-350 km/s (q54_chem conventions); Tb ratio = Jy ratio * 1.0126. Output out/q55/r10_core.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q54_chem as c
A='/archive/files/2016.1.00965.S/member.uid___A001_X879_X94a.NGC_3256_sci.spw27.cube.I.pbcor.fits'
c.L['HCN1-0_a']=(A,88.631602); c.L['HCO+1-0_a']=(A,89.188526)
rows=[]
for bm in (1.0,1.4,1.9,2.5):
    for reg in 'NS':
        h=c.meas('HCN1-0_a',bm,reg); o=c.meas('HCO+1-0_a',bm,reg)
        r=h[0]/o[0]*1.0126; e=r*np.hypot(h[1]/h[0],o[1]/o[0])
        rows.append(dict(beam=bm,reg=reg,IHCN=h[0],eHCN=h[1],IHCO=o[0],eHCO=o[1],R10_Tb=r,err=e)); print(bm,reg,round(h[0],1),round(o[0],1),round(r,3),round(e,3),flush=True)
pd.DataFrame(rows).round(4).to_csv('/workspace/out/q55/r10_core.csv',index=False)
