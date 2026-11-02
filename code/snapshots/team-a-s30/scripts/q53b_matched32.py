"""Q53 r2: HCN/HCO+ 3-2 at N and S in common circular beams (channel-wise image-plane convolution, q4_lib.spectrum),
Xb4b (0.25" native) for 0.3-1.4", Xb4d (0.93x0.74) for >=1.1" as uv-consistency check. 1" = 188 pc. Gaussian fit to each line at the fixed N/S position.
Output out/q53/matched32.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from q4_lib import spectrum,kb,vel
from q14a_beamseq_funcs import bright
rows=[]
LN={'HCO+':(267.557626,{'Xb4b':'spw25','Xb4d':'spw25'}),'HCN':(265.886434,{'Xb4b':'spw27','Xb4d':'spw27'})}
for ous,beams in (('Xb4b',[0.30,0.40,0.55,0.75,1.0,1.4,1.9]),('Xb4d',[1.1,1.4,1.9,2.5])):
    for bm in beams:
        for reg in ('N','S'):
            I={}
            for ln,(f0,sp) in LN.items():
                tb=(bm,bm,0.)
                try: fr,s,beam,nb=spectrum(ous,sp[ous],reg,tb=tb)
                except ValueError: break
                v=vel(fr,f0); y=s*kb(fr.mean(),bm,bm)
                p,eA,sg=bright(v,y); I[ln]=(p[0],eA*1.3,p[1],p[2])
            if len(I)<2: continue
            r=I['HCN'][0]/I['HCO+'][0]; e=r*np.hypot(I['HCN'][1]/I['HCN'][0],I['HCO+'][1]/I['HCO+'][0])
            rows.append(dict(ous=ous,beam_as=bm,pc=bm*188,reg=reg,IHCN=I['HCN'][0],eHCN=I['HCN'][1],IHCO=I['HCO+'][0],eHCO=I['HCO+'][1],ratio=r,err=e,v0=I['HCN'][2],fw=I['HCN'][3]))
            print(rows[-1])
pd.DataFrame(rows).to_csv('/workspace/out/q53/matched32.csv',index=False)
