"""Q53 r3b: HCN/HCO+ 4-3 at N, S from archive cubes; common fixed velocity window (same width both lines), line-free linear baseline,
noise from baseline channels; beams 1.0,1.9 (and 4.5" to compare with Xb47 4.1x2.7" same-project cubes = uv/flux-recovery check; 0.75 min usable);
position offsets (pointing) at 1.0". T_B in K km/s (RJ). Output out/q53/fixedwin43.csv, spectra npz out/q53/spec43.npz"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4_lib
from q4_lib import spectrum,kb,vel,REG
P={'HCN':'/archive/files/2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw29.cube.I.pbcor.fits',
   'HCO+':'/archive/files/2018.1.00223.S/member.uid___A001_X133d_X858.NGC3256_sci.spw27.cube.I.pbcor.fits',
   'HCN47':'/data/member.uid___A001_X1273_Xb47.ngc_3256_sci.spw16.cube.I.pbcor.fits',
   'HCO47':'/data/member.uid___A001_X1273_Xb47.ngc_3256_sci.spw18.cube.I.pbcor.fits'}
F0={'HCN':354.505477,'HCO+':356.734223}; F0['HCN47']=F0['HCN']; F0['HCO47']=F0['HCO+']
VC={'N':2785.,'S':2750.}; HW=350.; BL=(450.,950.)
cur=[None]; q4_lib.f=lambda *a,**k: P[cur[0]]
def integ(v,y,vc):
    w=np.abs(v-vc)<HW; b=(np.abs(v-vc)>BL[0])&(np.abs(v-vc)<BL[1])&np.isfinite(y)
    c=np.polyfit(v[b]-vc,y[b],1); r=y[b]-np.polyval(c,v[b]-vc); sg=1.4826*np.median(np.abs(r-np.median(r)))
    dv=abs(np.median(np.diff(v))); I=np.sum(y[w]-np.polyval(c,v[w]-vc))*dv; e=sg*dv*np.sqrt(w.sum())
    # baseline-fit uncertainty on the integral (linear fit covariance, approx)
    e=np.hypot(e,sg*dv*w.sum()/np.sqrt(max(b.sum(),1)))
    return I,e,sg,b.sum()
def meas(ln,reg,bm,pos=None):
    cur[0]=ln; tb=(bm,bm,0.) if bm else None
    fr,s,beam,nb=spectrum('x','x',reg,tb=tb,pos=pos); v=vel(fr,F0[ln]); k=kb(fr.mean(),*(beam[:2]))
    y=s*k; I,e,sg,nbl=integ(v,y,VC[reg]); return I,e,sg,nbl,v,y
rows=[]; spec={}
for bm,pairs in ((1.0,[('HCN','HCO+')]),(1.9,[('HCN','HCO+')]),(4.5,[('HCN','HCO+'),('HCN47','HCO47')])):
  for reg in ('N','S'):
    for a,b in pairs:
        A=meas(a,reg,bm); B=meas(b,reg,bm); spec[f'{a}_{reg}_{bm}']=(A[4],A[5]); spec[f'{b}_{reg}_{bm}']=(B[4],B[5])
        r=A[0]/B[0]; er=r*np.hypot(A[1]/A[0],B[1]/B[0])
        rows.append(dict(set=a+'/'+b,beam=bm,reg=reg,IHCN=A[0],eHCN=A[1],IHCO=B[0],eHCO=B[1],ratio=r,err=er,sgHCN=A[2],sgHCO=B[2]))
        print(rows[-1],flush=True)
pd.DataFrame(rows).round(4).to_csv('/workspace/out/q53/fixedwin43.csv',index=False)
np.savez('/workspace/out/q53/spec43.npz',**{k+'_v':v[0] for k,v in spec.items()},**{k+'_y':v[1] for k,v in spec.items()})
# pointing / peak position test at 1.0"
off=np.arange(-0.4,0.41,0.2); rows=[]
for reg in ('N','S'):
    ra0,de0=REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']
    for dx in off:
      for dy in off:
        pos=(ra0+dx/3600/np.cos(np.radians(de0)),de0+dy/3600)
        A=meas('HCN',reg,1.0,pos); B=meas('HCO+',reg,1.0,pos)
        rows.append(dict(reg=reg,dx=dx,dy=dy,IHCN=A[0],IHCO=B[0],ratio=A[0]/B[0]))
T=pd.DataFrame(rows); T.round(4).to_csv('/workspace/out/q53/posoff43.csv',index=False)
for reg in 'NS':
    t=T[T.reg==reg]; print(reg,'ratio range',t.ratio.min().round(3),t.ratio.max().round(3),'std',t.ratio.std().round(3),
      'HCN peak at',t.loc[t.IHCN.idxmax(),['dx','dy']].values,'ratio there',t.loc[t.IHCN.idxmax(),'ratio'].round(3),
      'HCO peak at',t.loc[t.IHCO.idxmax(),['dx','dy']].values,'ratio there',t.loc[t.IHCO.idxmax(),'ratio'].round(3))
