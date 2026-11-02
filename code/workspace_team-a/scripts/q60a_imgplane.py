"""Q60 step 1 (image plane): HCN & HCO+ 4-3 integrals (fixed window VC+-350 km/s, linear baseline 450-950 km/s) at N and S from every 4-3 archive cube
(Xb47 ACA single tuning; 2018.1.00493 X1df1 HCN 0.7"; 2018.1.00223 X858 HCO+ 0.3") convolved to common circular beams (T_B, K km/s).
Same line, same beam, different array = flux-scale / missing-flux test. Output out/q60/imgplane.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4_lib
from q4_lib import spectrum,kb,vel
VC={'N':2785.,'S':2750.}; HW=350.; BL=(450.,950.)
def win(v,y,vc):
    b=(np.abs(v-vc)>BL[0])&(np.abs(v-vc)<BL[1])&np.isfinite(y); c=np.polyfit(v[b]-vc,y[b],1)
    r=y[b]-np.polyval(c,v[b]-vc); sg=1.4826*np.median(np.abs(r-np.median(r))); dv=abs(np.median(np.diff(v)))
    w=(np.abs(v-vc)<HW); n=w.sum(); I=np.sum(y[w]-np.polyval(c,v[w]-vc))*dv
    return I,np.hypot(sg*dv*np.sqrt(n),sg*dv*n/np.sqrt(b.sum()))
SRC={('HCN','Xb47'):('Xb47','spw16',None),('HCO+','Xb47'):('Xb47','spw18',None),
 ('HCN','X1df1'):('x','x','/archive/files/2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw29.cube.I.pbcor.fits'),
 ('HCO+','X858'):('x','x','/archive/files/2018.1.00223.S/member.uid___A001_X133d_X858.NGC3256_sci.spw27.cube.I.pbcor.fits')}
F0={'HCN':354.505477,'HCO+':356.734223}; orig=q4_lib.f; cur=[None]
rows=[]
for bm in [4.5,5.5,7.0]:
  for (ln,tag),(ous,spw,path) in SRC.items():
    for reg in 'NS':
      if path: q4_lib.f=lambda *a,**k: path
      else: q4_lib.f=orig
      try: fr,s,beam,nb=spectrum(ous,spw,reg,tb=(bm,bm,0.))
      except ValueError: continue
      v=vel(fr,F0[ln]); y=s*kb(fr.mean(),bm,bm); I,e=win(v,y,VC[reg])
      rows.append(dict(beam=bm,line=ln,src=tag,reg=reg,I=I,e=e)); print(rows[-1],flush=True)
d=pd.DataFrame(rows); d.round(3).to_csv('/workspace/out/q60/imgplane.csv',index=False)
p=d.pivot_table(index=['beam','reg'],columns=['line','src'],values='I'); print(p.round(1))
