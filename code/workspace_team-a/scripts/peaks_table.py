# Peaks of each tracer's moment-0 map (vsys+-350 km/s, pbcor) within the primary-beam FWHM; S/N from per-pixel line-free rms
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from lineset import RA0,DE0
from moments import mom, LINES
from scipy import ndimage as ndi
import pandas as pd
RMAX={'Xb55':17.,'Xb47':9.,'Xb4f':11.,'Xb4d':11.,'Xb4b':11.}
rows=[]
for ous,spw,name,f0 in LINES:
    if ous=='Xb4b':
        dd=np.load(f'/workspace/out/m0_Xb4b_{spw}.npz'); m0=dd['m0_core']; e0=dd['sig']*dd['dv']*np.sqrt(dd['n_core']); x0,y0=int(dd['x0']),int(dd['y0'])
        hd=fits.getheader(f(ous,spw))
    else:
        hd,m0,e0,dv,n=mom(ous,spw,name,f0); x0=y0=0
    w=WCS(hd).celestial; pix=abs(hd['CDELT2'])*3600; bmaj=hd['BMAJ']*3600
    ny,nx=m0.shape; Y,X=np.mgrid[:ny,:nx]; ra,de=w.pixel_to_world_values(X+x0,Y+y0)
    dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600; r=np.hypot(dx,dy)
    e1=np.where(np.isfinite(e0),e0,np.nanmedian(e0)); snr=np.nan_to_num(m0)/ndi.uniform_filter(e1,max(3,int(bmaj/pix)))
    fp=max(3,int(hd['BMIN']*3600/pix)) if ous!='Xb4b' else int(0.6/pix)
    lm=(m0==ndi.maximum_filter(np.nan_to_num(m0),size=fp))&(snr>8)&(r<RMAX[ous])
    ys,xs=np.where(lm)
    for y,x in zip(ys,xs):
        K=1.222e6/((f0*(1-2780/C))**2*bmaj*hd['BMIN']*3600)
        rows.append(dict(ous=ous,spw=spw,line=name,snr=float(snr[y,x]),I_m0_Jybkms=float(m0[y,x]),I_m0_Kkms=float(m0[y,x]*K),dRA=float(dx[y,x]),dDec=float(dy[y,x]),ra=float(ra[y,x]),dec=float(de[y,x])))
df=pd.DataFrame(rows).sort_values(['ous','line','snr'],ascending=[True,True,False]); df.to_csv('/workspace/out/tracer_peaks.csv',index=False)
pd.set_option('display.width',200,'display.max_rows',300); print(df.round(2).to_string())
