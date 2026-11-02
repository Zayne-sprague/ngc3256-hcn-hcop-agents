import numpy as np
from astropy.io import fits
D='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.%s.cube.I.pbcor.fits'
c=299792.458
def load(spw):
    h=fits.open(D%spw,memmap=True)
    d=h[0].data[0] if h[0].data.ndim==4 else h[0].data
    H=h[0].header
    f=H['CRVAL3']+(np.arange(H['NAXIS3'])-H['CRPIX3']+1)*H['CDELT3']
    return d,f,H
def run(spw,f0,pos):
    d,f,H=load(spw)
    v=c*(f0/f-1)  # optical-ish; radio below
    vr=c*(1-f/f0*1.0)
    # sky freq of rest f0 at vsys
    res={}
    for name,(dx,dy) in pos.items():
        ra=H['CRVAL1'];dec=H['CRVAL2']
        x=int(round(H['CRPIX1']-1+dx/3600/H['CDELT1']/1)) if False else None
    return d,f,H
for spw,f0 in [('spw25',267.557626),('spw27',265.886434)]:
    d,f,H=load(spw)
    print(spw,d.shape,H['CRPIX1'],H['CDELT1']*3600,H['CRVAL1'],H['CRVAL2'])
    vel=c*(1-f/(f0))  # radio vs rest in GHz (f in Hz)
