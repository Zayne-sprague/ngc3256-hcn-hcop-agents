import numpy as np,pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
D='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.%s.cube.I.pbcor.fits'
reg=pd.read_csv('/workspace/out/regions.csv'); print(reg[['region','dRA_arcsec','dDec_arcsec']].to_string())
c=299792.458
def mom(spw,f0,x,y,hw=350):
    h=fits.open(D%spw,memmap=True);hd=h[0].header
    n=hd['NAXIS3'];nu=hd['CRVAL3']+(np.arange(n)-hd['CRPIX3']+1)*hd['CDELT3']
    v=c*(1-nu/1e9/f0);dv=abs(v[1]-v[0])
    w=WCS(hd).celestial
    px,py=w.world_to_pixel_values(156.96346+x/3600/np.cos(np.deg2rad(-43.90461)),-43.90461+y/3600)
    px,py=int(round(float(px))),int(round(float(py)))
    a=np.nansum(np.array(h[0].data[0,:,py-1:py+2,px-1:px+2],float),axis=(1,2))/9  # 3x3 mean
    m=abs(v-2780)<hw; off=(abs(v-2780)>700)&(abs(v-2780)<1100)
    return (a[m]-np.median(a[off])).sum()*dv, np.std(a[off])*dv*np.sqrt(m.sum())
for r in reg.itertuples():
    if r.region not in('N','S','Bridge'):continue
    I1,e1=mom('spw27',265.886434,r.dRA_arcsec,r.dDec_arcsec);I2,e2=mom('spw25',267.557626,r.dRA_arcsec,r.dDec_arcsec)
    k=(267.557626/265.886434)**2  # Tb: nu^-2 conv, same beams (approx; different beams ~ few %)
    print(r.region,I1,I2,'Fnu ratio',I1/I2,'Tb ratio',I1/I2*k,'+-',I1/I2*k*np.hypot(e1/I1,e2/I2))
