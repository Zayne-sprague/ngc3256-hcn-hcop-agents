import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
L=pd.read_csv('/workspace/out/line_table_v2.csv')
f='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw25.cube.I.pbcor.fits'
h=fits.open(f,memmap=True); d=h[0].data[0] if h[0].data.ndim==4 else h[0].data
w=WCS(h[0].header).celestial; H=h[0].header
nu=H['CRVAL3']+(np.arange(d.shape[0])-(H['CRPIX3']-1))*H['CDELT3']
bm=h[1].data if len(h)>1 else None
for reg in ['N','S','Bridge']:
    r=L[(L.ous=='Xb4d')&(L.region==reg)&(L.kind=='peakpix')&(L.line=='HCN 3-2')].iloc[0]
    x,y=w.world_to_pixel_values(r.ra,r.dec); x=int(round(float(x)));y=int(round(float(y)))
    out={}
    for ln,f0,z in [('HCN',265.886434,0),('HCO+',267.557626,0)]:
        for zz in [0.00935]:
            fs=f0*1e9/(1+zz)
            ch=np.argmin(abs(nu-fs)) 
            sp=d[:,y,x]
            sel=abs(nu-fs)<fs*500/3e5   # +-500 km/s
            k=1.222e6/(f0**2*r.bmaj*r.bmin)
            out[ln]=np.nanmax(sp[sel])*1e3*k*1e-3
    print(reg,x,y,out,out['HCN']/out['HCO+'])
