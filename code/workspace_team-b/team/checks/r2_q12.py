import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS
d='/data/member.uid___A001_X1273_'
def load(o,f):
    h=fits.open(d+o+'.ngc_3256_sci.'+f);im=np.squeeze(h[0].data);hd=h[0].header
    return im,hd
ra,de=156.96346,-43.90461
res={}
for o in['Xb4f','Xb4d']:
    im,hd=load(o,'spw25_27_29_31.cont.I.pbcor.fits' if o=='Xb4d' else 'spw16_18_20_22.cont.I.pbcor.fits')
    w=WCS(hd).celestial;x,y=w.world_to_pixel_values(ra,de);pix=abs(hd['CDELT2'])*3600
    bm=np.sqrt(hd['BMAJ']*hd['BMIN'])*3600
    im=np.nan_to_num(im)
    # smooth to 10" beam; flux in Jy/10"beam
    s=np.sqrt(max(10**2-bm**2,0))/2.355/pix
    sm=gaussian_filter(im,s)
    yy,xx=np.indices(im.shape);r=np.hypot(xx-x,yy-y)*pix
    # integrated flux: sum(image Jy/beam)/beam_area_pix for r<R, original image
    ba=1.1331*hd['BMAJ']*hd['BMIN']*3600**2/pix**2
    res[o]=[im[r<R].sum()/ba*1e3 for R in(4,8,12)]
    print(o,bm,res[o])
print([a/b for a,b in zip(res['Xb4d'],res['Xb4f'])])
