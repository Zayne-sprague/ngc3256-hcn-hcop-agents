import numpy as np,warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
f='/archive/files/2016.1.00965.S/member.uid___A001_X879_X94a.NGC_3256_sci.spw29.cube.I.pbcor.fits'
h=fits.getheader(f); d=fits.getdata(f)[0]; print(d.shape)
nu=h['CRVAL3']+h['CDELT3']*(np.arange(d.shape[0]))
w=WCS(h).celestial
print('freq',nu[0]/1e9,nu[-1]/1e9)
RA0,DE0=156.96346,-43.90461
for nm,(dx,dy) in {'N':(0,2.63),'S':(-0.15,-2.64)}.items():
    ra=RA0+dx/3600/np.cos(np.deg2rad(DE0)); de=DE0+dy/3600
    x,y=w.world_to_pixel_values(ra,de); x=int(round(float(x)));y=int(round(float(y)))
    s=d[:,y-2:y+3,x-2:x+3].mean((1,2))*1e3
    v=299792.458*(99.02295e9/nu-1)  # optical? use radio below
    print(nm,x,y)
    for i in range(0,len(nu),6): print(f'{nu[i]/1e9:.3f} {s[i:i+6].round(2)}')
