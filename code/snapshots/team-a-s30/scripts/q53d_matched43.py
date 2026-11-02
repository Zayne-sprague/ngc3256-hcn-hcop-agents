"""Q53 r3: HCN/HCO+ 4-3 at N and S from archive cubes (2018.1.00493.S X1df1 spw29 HCN; 2018.1.00223.S X858 spw27 HCO+), channel-wise
convolution to common circular beams, T_B (RJ) moment-0 from Gaussian fit (same method as q53b). 1"=188 pc. Output out/q53/matched43.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4_lib
from q4_lib import spectrum,kb,vel
from q14a_beamseq_funcs import bright
P={'HCN':'/archive/files/2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw29.cube.I.pbcor.fits',
   'HCO+':'/archive/files/2018.1.00223.S/member.uid___A001_X133d_X858.NGC3256_sci.spw27.cube.I.pbcor.fits'}
F0={'HCN':354.505477,'HCO+':356.734223}
cur=[None]; q4_lib.f=lambda *a,**k: P[cur[0]]
rows=[]
for bm in (1.0,1.4,1.9,2.5):
  for reg in ('N','S'):
    I={}
    for ln in P:
        cur[0]=ln
        fr,s,beam,nb=spectrum('x','x',reg,tb=(bm,bm,0.))
        v=vel(fr,F0[ln]); y=s*kb(fr.mean(),bm,bm)
        try: p,eA,sg=bright(v,y); I[ln]=(p[0],eA*1.3,p[1],p[2])
        except Exception as e: print('fail',ln,bm,reg,e); break
    if len(I)<2: continue
    r=I['HCN'][0]/I['HCO+'][0]; e=r*np.hypot(I['HCN'][1]/I['HCN'][0],I['HCO+'][1]/I['HCO+'][0])
    rows.append(dict(beam_as=bm,pc=bm*188,reg=reg,IHCN=I['HCN'][0],eHCN=I['HCN'][1],IHCO=I['HCO+'][0],eHCO=I['HCO+'][1],ratio=r,err=e,v0=I['HCN'][2],fwHCN=I['HCN'][3],v0O=I['HCO+'][2],fwHCO=I['HCO+'][3]))
    print({k:(round(x,3) if not isinstance(x,str) else x) for k,x in rows[-1].items()})
pd.DataFrame(rows).to_csv('/workspace/out/q53/matched43.csv',index=False)
