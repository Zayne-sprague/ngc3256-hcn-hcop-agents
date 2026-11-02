import numpy as np, sys
from astropy.io import fits
f='/data/member.uid___A001_X1273_Xb55.ngc_3256_sci.spw%s.cube.I.pbcor.fits'
c=299792.458
def spec(spw,ra=156.96346,dec=-43.9046):
    h=fits.open(f%spw,memmap=True); d=h[0].data[0] if h[0].data.ndim==4 else h[0].data
    H=h[0].header
    nu=H['CRVAL3']+(np.arange(H['NAXIS3'])-H['CRPIX3']+1)*H['CDELT3']
    from astropy.wcs import WCS
    w=WCS(H).celestial; x,y=w.all_world2pix(ra,dec,0); x=int(round(float(x)));y=int(round(float(y)))
    return nu,d[:,y,x]*1e3,d[:,y-12:y+12,x-12:x+12]
for spw,rest,name in [('18',178.972,'HOC+ 2-1'),('22',165.07,'CH3OH J1-3')]:
    nu,s,_=spec(spw)
    # try both frequency conv: rest at z
    for v0 in [2750,2780]:
        fc=rest/(1+v0/c)*1.0
        fc=rest*(1-v0/c)
        m=np.abs((nu/1e9-fc)*1e3/ (fc*1e3)*c)<150
        dv=abs(nu[1]-nu[0])/nu[0]*c
        off=(np.abs((nu/1e9-fc)/fc*c)>400)&(np.abs((nu/1e9-fc)/fc*c)<1500)
        I=(s[m].sum())*dv/1e3; rms=s[off].std()
        print(name,v0,'nchan',m.sum(),'I(Jy/b km/s)=%.3f'%I,'rms/ch mJy %.2f'%rms,'sig %.1f'%(I*1e3/dv/(rms*np.sqrt(m.sum()))*1.0))
