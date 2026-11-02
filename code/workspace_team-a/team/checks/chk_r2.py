import numpy as np
from astropy.io import fits
f='/data/member.uid___A001_X1273_Xb55.ngc_3256_sci.spw22.cube.I.pbcor.fits'
h=fits.open(f);d=np.squeeze(h[0].data);hd=h[0].header
n=hd['NAXIS3'] if 'NAXIS3' in hd else 0
nu=hd['CRVAL3']+(np.arange(d.shape[0])-hd['CRPIX3']+1)*hd['CDELT3']
c=299792.458
# CH3OH Q J1-3 centre 165.0611 rest; v radio
v=c*(1-nu/1e9/165.0611*1.0)  # nu in Hz
v=c*(1-nu/165.0611e9)
cube=d
iy,ix=42,42; s=np.nansum(cube[(abs(v-2780)<300)],0); sub=s[34:50,34:50]; a,b=np.unravel_index(np.nanargmax(sub),sub.shape); iy,ix=34+a,34+b
sp=cube[:,iy,ix]
m=abs(v-2775)<250
off=(abs(v-2775)>600)&(abs(v-2775)<1500)
dv=abs(v[1]-v[0])
I=np.nansum(sp[m]-np.nanmedian(sp[off]))*dv
rms=np.nanstd(sp[off]);e=rms*dv*np.sqrt(m.sum())
print(iy,ix,I,e,I/e, 'chan rms',rms)
