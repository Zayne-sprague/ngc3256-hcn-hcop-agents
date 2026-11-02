import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
D='/data/member.uid___A001_X1273_%s.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'
D2='/data/member.uid___A001_X1273_Xb4f.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits'
def ld(f):
    h=fits.open(f);return np.squeeze(h[0].data),h[0].header
tf,hf=ld(D2); tb=(hf['BMAJ']*3600,hf['BMIN']*3600)
pts={'N':(156.96346,-43.90390),'S':(156.96340,-43.90532)}
def pk(a,h,p):
    x,y=WCS(h).celestial.world_to_pixel_values(*p);x=int(round(float(x)));y=int(round(float(y)))
    return a[y-5:y+6,x-5:x+6].max() if h['CDELT1']<1e-4 else a[y-1:y+2,x-1:x+2].max()
for o in ('Xb4d','Xb4b'):
    a,h=ld(D%o); ps=abs(h['CDELT1'])*3600; b0=(h['BMAJ']*3600,h['BMIN']*3600)
    # circular-equivalent smoothing (PA ignored, geometric-mean approx)
    s0=np.sqrt(b0[0]*b0[1]);s1=np.sqrt(tb[0]*tb[1]);sg=np.sqrt(s1**2-s0**2)/2.3548/ps
    a=np.nan_to_num(a);b=gaussian_filter(a,sg)*s1**2/s0**2
    print(o,{k:round(float(pk(b,h,p)/pk(tf,hf,p)),2) for k,p in pts.items()}, 'sum all /Xb4f', )
