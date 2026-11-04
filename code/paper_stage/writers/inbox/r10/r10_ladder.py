"""r10: HCN & HCO+ 4-3 vs 3-2 T_b (K km/s, fixed window VC+-350, linear baseline 450-950) at N and S at matched circular beams 1.0" and 1.5".
Reuses q53g_fixedwin_all.py machinery. 3-2: Xb4b (and Xb4d at 1.5"); 4-3: 2018.1.00493 X1df1 (HCN) + 2018.1.00223 X858 (HCO+).
Missing flux: 4-3 high-res/Xb47 ratio at 4.5" from out/q60/imgplane.csv. Run: python3 -I r10_ladder.py"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4_lib
from q4_lib import spectrum,kb,vel
VC={'N':2785.,'S':2750.}; HW=350.; BL=(450.,950.)
P43={'HCN':'/archive/files/2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw29.cube.I.pbcor.fits',
     'HCO+':'/archive/files/2018.1.00223.S/member.uid___A001_X133d_X858.NGC3256_sci.spw27.cube.I.pbcor.fits'}
F43={'HCN':354.505477,'HCO+':356.734223}; F32={'HCO+':267.557626,'HCN':265.886434}; SP={'HCO+':'spw25','HCN':'spw27'}
orig=q4_lib.f; cur=[None]
def win(v,y,vc):
    b=(np.abs(v-vc)>BL[0])&(np.abs(v-vc)<BL[1])&np.isfinite(y); c=np.polyfit(v[b]-vc,y[b],1)
    r=y[b]-np.polyval(c,v[b]-vc); sg=1.4826*np.median(np.abs(r-np.median(r))); dv=abs(np.median(np.diff(v)))
    w=np.abs(v-vc)<HW; n=w.sum(); I=np.sum(y[w]-np.polyval(c,v[w]-vc))*dv
    return I,np.hypot(sg*dv*np.sqrt(n),sg*dv*n/np.sqrt(b.sum()))
def run(J,ous,bm,reg,ln):
    if J==4: cur[0]=ln; q4_lib.f=lambda *a,**k: P43[cur[0]]; fr,s,_,_=spectrum('x','x',reg,tb=(bm,bm,0.)); f0=F43[ln]
    else: q4_lib.f=orig; fr,s,_,_=spectrum(ous,SP[ln],reg,tb=(bm,bm,0.)); f0=F32[ln]
    return win(vel(fr,f0),s*kb(fr.mean(),bm,bm),VC[reg])
rows=[]
for bm,ous32 in [(1.0,['Xb4b']),(1.5,['Xb4b','Xb4d'])]:
  for reg in 'NS':
    for ln in ('HCN','HCO+'):
      I4,e4=run(4,'x',bm,reg,ln)
      for o in ous32:
        I3,e3=run(3,o,bm,reg,ln); rows.append(dict(beam=bm,reg=reg,line=ln,ous32=o,I43=I4,e43=e4,I32=I3,e32=e3,r43_32=I4/I3,e=abs(I4/I3)*np.hypot(e4/I4,e3/I3)))
        print(rows[-1],flush=True)
pd.DataFrame(rows).round(4).to_csv('ladder.csv',index=False)
