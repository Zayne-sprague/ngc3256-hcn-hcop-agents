# list curated-catalogue lines falling inside each spw at vsys=2780 km/s (radio, LSRK)
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import pandas as pd, numpy as np
from curated import curated
C=299792.458; V=2780.
SPW={'B5 spw22':(163.244,165.110),'B5 spw20':(164.536,166.402),'B5 spw16':(174.607,176.469),'B5 spw18':(176.366,178.227),
'B7 spw22':(338.871,340.668),'B7 spw20':(340.572,342.396),'B7 spw16':(350.622,352.423),'B7 spw18':(352.379,354.180),
'B6 spw31/22':(248.981,250.845),'B6 spw29/20':(250.738,252.602),'B6 spw27/18':(262.808,264.672),'B6 spw25/16':(264.565,266.429)}
df=curated(250)
rows=[]
for k,(lo,hi) in SPW.items():
    f0lo,f0hi=lo/(1-V/C),hi/(1-V/C)
    s=df[(df.f0>=f0lo)&(df.f0<=f0hi)].copy(); s['spw']=k; s['fsky']=s.f0*(1-V/C); rows.append(s)
al=pd.concat(rows)
al.to_csv('/workspace/out/inband_lines_curated.csv',index=False)
pd.set_option('display.width',200,'display.max_rows',500)
for k in SPW:
    s=al[(al.spw==k)&(al.logA>-4.3)&(al.Eu<200)].sort_values('f0')
    print('\n==',k); print(s[['sp','qn','f0','fsky','Eu','logA']].to_string(index=False))
