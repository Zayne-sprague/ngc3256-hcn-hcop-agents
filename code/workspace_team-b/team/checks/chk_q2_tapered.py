import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
d='/data/member.uid___A001_X1273_Xb4d_ngc_3256_sci.spw27.cube.I.tapered.image.pbcor.fits'
h=fits.open(d); c=h[0].data.squeeze().astype(float); hd=h[0].header
print(c.shape, hd.get('BMAJ')*3600, hd.get('BMIN')*3600)
w=WCS(hd).celestial
f=(np.arange(c.shape[0])-(hd['CRPIX3']-1))*hd['CDELT3']+hd['CRVAL3']
rest=265.88618e9; z=2780/299792.458
v=(rest/(1+0)-f)/ (rest) *299792.458  # placeholder
f0=rest/(1+2780/299792.458) # radio approx handled below
v=(1-f/rest)*299792.458-2780
print(v.min(),v.max(),np.diff(v)[:2])
# S nucleus position from team regions
import json,re,sys
sys.path.insert(0,'/workspace/scripts')
from s72_lib import REGIONS
S=REGIONS['S'][:2]; print(S)
from astropy.coordinates import SkyCoord
def pix(ra,dec): return w.world_to_pixel_values(ra,dec)
ps=abs(hd['CDELT2'])*3600
def ap(ra,dec,r):
    x0,y0=pix(ra,dec); yy,xx=np.indices(c.shape[1:]); return np.hypot(xx-x0,yy-y0)*ps<r
Sra,Sdec=S
dec_l=Sdec-1.88/3600; ra_l=Sra+0.0
bm=hd['BMAJ']*3600; bn=hd['BMIN']*3600; nb=1.1331*bm*bn/ps**2
dv=abs(np.diff(v).mean())
# baseline from v<-450 channels (few); report spectrum
for nm,(ra,dec) in dict(lobe=(Sra,dec_l),S=(Sra,Sdec),off=(Sra+12/3600/np.cos(np.radians(Sdec)),Sdec-6/3600)).items():
    m=ap(ra,dec,1.5)
    s=np.nanmean(c[:,m],axis=1)*m.sum()/nb
    print(nm,' '.join(f'{x*1e3:.1f}' for x in s))
    k=(v>=-405)&(v<=-150); print(nm,'wing flux Jy km/s',s[k].sum()*dv, 'chan',k.sum())
print(np.round(v,0))
