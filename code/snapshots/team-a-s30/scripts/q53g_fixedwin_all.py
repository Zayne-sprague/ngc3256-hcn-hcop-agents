"""Q53 r3b: fixed common velocity window (VC +-350 km/s; N 2785, S 2750; baseline 450-950 km/s linear) integrals of HCN & HCO+ in T_B, for 3-2 (Xb4b/Xb4d) and 4-3
(X1df1/X858) at common circular beams; also split blue (<2780) / red (>2780) halves (spectra are double-horned: Gaussian is wrong). Output out/q53/fixedwin_all.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4_lib
from q4_lib import spectrum,kb,vel
VC={'N':2785.,'S':2750.}; HW=350.; BL=(450.,950.); VS=2780.
P43={'HCN':'/archive/files/2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw29.cube.I.pbcor.fits',
     'HCO+':'/archive/files/2018.1.00223.S/member.uid___A001_X133d_X858.NGC3256_sci.spw27.cube.I.pbcor.fits'}
F43={'HCN':354.505477,'HCO+':356.734223}; F32={'HCO+':267.557626,'HCN':265.886434}; SP={'HCO+':'spw25','HCN':'spw27'}
orig=q4_lib.f; cur=[None]
def win(v,y,vc,lo,hi):
    b=(np.abs(v-vc)>BL[0])&(np.abs(v-vc)<BL[1])&np.isfinite(y); c=np.polyfit(v[b]-vc,y[b],1)
    r=y[b]-np.polyval(c,v[b]-vc); sg=1.4826*np.median(np.abs(r-np.median(r))); dv=abs(np.median(np.diff(v)))
    w=(v>=lo)&(v<hi); n=w.sum(); I=np.sum(y[w]-np.polyval(c,v[w]-vc))*dv
    return I,np.hypot(sg*dv*np.sqrt(n),sg*dv*n/np.sqrt(b.sum()))
def run(J,ous,bm,reg,ln):
    if J==4: cur[0]=ln; q4_lib.f=lambda *a,**k: P43[cur[0]]; fr,s,beam,nb=spectrum('x','x',reg,tb=(bm,bm,0.)); f0=F43[ln]
    else: q4_lib.f=orig; fr,s,beam,nb=spectrum(ous,SP[ln],reg,tb=(bm,bm,0.)); f0=F32[ln]
    v=vel(fr,f0); y=s*kb(fr.mean(),bm,bm); vc=VC[reg]
    return {'all':win(v,y,vc,vc-HW,vc+HW),'blue':win(v,y,vc,vc-HW,VS),'red':win(v,y,vc,VS,vc+HW)}
rows=[]
CFG=[(3,'Xb4b',[0.3,0.4,0.55,0.75,1.0,1.4,1.9]),(3,'Xb4d',[1.1,1.4,1.9,2.5]),(4,'x',[1.0,1.4,1.9,2.5])]
for J,ous,beams in CFG:
  for bm in beams:
    for reg in 'NS':
      try: H=run(J,ous,bm,reg,'HCN'); O=run(J,ous,bm,reg,'HCO+')
      except ValueError: continue
      for k in ('all','blue','red'):
        r=H[k][0]/O[k][0]; e=abs(r)*np.hypot(H[k][1]/H[k][0],O[k][1]/O[k][0])
        rows.append(dict(J=J,ous=ous,beam=bm,pc=bm*188,reg=reg,win=k,IHCN=H[k][0],eHCN=H[k][1],IHCO=O[k][0],eHCO=O[k][1],ratio=r,err=e))
      print(J,ous,bm,reg,*[round(rows[-3+i]['ratio'],3) for i in range(3)],flush=True)
pd.DataFrame(rows).round(4).to_csv('/workspace/out/q53/fixedwin_all.csv',index=False)
