"""Q54: chemistry ratios (HNC/HCN, CN/HCN) at N and S in matched circular beams from archive cubes.
Fixed velocity window (N 2785, S 2750 km/s, +-350 km/s) after a linear baseline (450-950 km/s from VC), Tb integrals (K km/s).
Usage: python3 q54_chem.py ; output out/q54/chem.csv. Lines whose cube is missing are skipped."""
import sys,os,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q4_lib
from q4_lib import spectrum,kb,vel
A='/archive/files/'
VC={'N':2785.,'S':2750.}; HW=350.; BL=(450.,950.)
# line: (file, rest GHz)
L={'HCO+1-0_00993':(A+'2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits',89.188526),
   'HCO+3-2_Xb4d':('/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw25.cube.I.pbcor.fits',267.557626),
   'HCN3-2_Xb4d':('/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw27.cube.I.pbcor.fits',265.886434),
   'HCN1-0_00965':(A+'2016.1.00965.S/member.uid___A001_X879_X94a.NGC_3256_sci.spw27.cube.I.pbcor.fits',88.631602),
   'HNC1-0_00965':(A+'2016.1.00965.S/member.uid___A001_X879_X94e.NGC_3256_sci.spw25.cube.I.pbcor.fits',90.663568),
   'HCN1-0_00993':(A+'2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits',88.631602),
   'HNC1-0_00993':(A+'2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw1_89515MHz.12m.cube.I.pbcor.fits',90.663568),
   'HCN4-3_X1df1':(A+'2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw29.cube.I.pbcor.fits',354.505477),
   'CN3-2_X1df1':(A+'2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw25.cube.I.pbcor.fits',340.2477),  # N=3-2 J=7/2-5/2 group; see windows
   'HCN1-0_00412':(A+'2015.1.00412.S/member.uid___A001_X33e_X145.ari_l.NGC_3256_sci.spw1_88144MHz.12m.cube.I.pbcor.fits',88.631602),
   'CN1-0_X4d':(A+'2018.1.00223.S/member.uid___A001_X133f_X4d.NGC3256_sci.spw27.cube.I.pbcor.fits',113.490982),
   'HNC1-0_00412':(A+'2015.1.00412.S/member.uid___A001_X33e_X149.NGC_3256_sci.spw25.cube.I.pbcor.fits',90.663568),
   'HNC3-2_00412':(A+'2015.1.00412.S/member.uid___A001_X33e_X169.ari_l.ngc_3256_sci.spw3_268777MHz.12m.cube.I.pbcor.fits',271.981142)}
def win(v,y,vc,lo,hi):
    b=(np.abs(v-vc)>BL[0])&(np.abs(v-vc)<BL[1])&np.isfinite(y); c=np.polyfit(v[b]-vc,y[b],1)
    r=y[b]-np.polyval(c,v[b]-vc); sg=1.4826*np.median(np.abs(r-np.median(r))); dv=abs(np.median(np.diff(v)))
    w=(v>=lo)&(v<hi); n=w.sum(); I=np.sum(y[w]-np.polyval(c,v[w]-vc))*dv
    return I,np.hypot(sg*dv*np.sqrt(n),sg*dv*n/np.sqrt(b.sum()))
def meas(key,bm,reg,hw=HW):
    path,f0=L[key]; q4_lib.f=lambda *a,**k: path
    fr,s,beam,nb=spectrum('x','x',reg,tb=(bm,bm,0.)); v=vel(fr,f0); y=s*kb(fr.mean(),bm,bm); vc=VC[reg]
    I,e=win(v,y,vc,vc-hw,vc+hw); sg=1.4826*np.median(np.abs(np.diff(y)))/np.sqrt(2)
    return I,e,nb,sg,v.min(),v.max()
if __name__=='__main__':
    rows=[]
    for key in L:
        if not os.path.exists(L[key][0]): print('missing',key); continue
        for bm in [0.7,1.0,1.4,1.9,2.5]:
            for reg in 'NS':
                try: I,e,nb,sg,v0,v1=meas(key,bm,reg)
                except ValueError: continue
                rows.append(dict(line=key,beam=bm,reg=reg,I=I,e=e,noise_chan_K=sg,nbeam=f'{nb[0]:.2f}x{nb[1]:.2f}',vmin=v0,vmax=v1))
        print(key,flush=True)
    pd.DataFrame(rows).round(3).to_csv('/workspace/out/q54/chem.csv',index=False)
