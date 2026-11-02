# Continuum table: N / S peak (mJy/beam, pbcor), r=2.5" aperture flux densities, additional >5 sigma sources, mfs vs line-free cont comparison
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
SRC={**NUC,'E1':(156.96608,-43.90350),'SE':(156.96685,-43.90635)}   # E1, SE = additional 0.24"-resolution continuum sources (s1_continuum)
rows=[]
for fn in files('*.pbcor.fits')+files('*tapered*cont*'):
    if 'cube' in fn or 'repBW' in fn: continue
    d,h=load2d(fn); w=WCS(h).celestial; ps=abs(h['CDELT1'])*3600; ny,nx=d.shape
    bmaj,bmin=h['BMAJ']*3600,h['BMIN']*3600; bpix=np.pi*bmaj*bmin/(4*np.log(2))/ps**2
    yy,xx=np.mgrid[:ny,:nx]
    mid=w.world_to_pixel_values(0.5*(NUC['N'][0]+NUC['S'][0]),0.5*(NUC['N'][1]+NUC['S'][1]))
    ri=np.hypot(xx-mid[0],yy-mid[1])*ps; fin=np.isfinite(d)
    ann=fin&(ri>max(6,1.5*bmaj))&(ri<max(12,3*bmaj)); v=d[ann]; rms=1.4826*np.median(np.abs(v-np.median(v)))   # off-nuclear MAD (includes faint sources; slightly conservative)
    r={'file':short(fn),'nu_GHz':(h.get('CRVAL3',np.nan)/1e9),'beam':f'{bmaj:.2f}x{bmin:.2f}','rms_mJy_b':rms*1e3}
    for k,(ra,dec) in SRC.items():
        x,y=w.world_to_pixel_values(ra,dec); r_pix=max(1.0,0.6*bmin/ps)
        m=np.hypot(xx-x,yy-y)<=max(r_pix,1)
        sub=np.where(m,d,-np.inf); iy,ix=np.unravel_index(np.argmax(sub),d.shape)
        r[k+'_peak_mJy_b']=d[iy,ix]*1e3
        if ps<0.5 and k in ('N','S','E1','SE'):
            ap=np.hypot(xx-x,yy-y)<=2.5/ps; r[k+'_S2.5_mJy']=np.nansum(d[ap])/bpix*1e3
    ap=np.hypot(xx-mid[0],yy-mid[1])<=8/ps; r['tot_r8_mJy']=np.nansum(d[ap&fin])/bpix*1e3
    rows.append(r)
df=pd.DataFrame(rows); df.to_csv('/workspace/out/s6_continuum.csv',index=False)
pd.set_option('display.width',250,'display.max_columns',30,'display.max_rows',100)
print(df.round(3).to_string())
