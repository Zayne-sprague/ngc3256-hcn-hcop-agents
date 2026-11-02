import numpy as np,glob,sys
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
D='/data/'
f0=267.55763e9; V=2780.
NUC={'N':(156.96346,-43.90388),'S':(156.96340,-43.90533)}
def m0(ous,spw,win,base=True):
    fn=glob.glob(D+f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]
    h=fits.open(fn,memmap=True); hd=h[0].header; w=WCS(hd)
    n=hd['NAXIS3']; nu=hd['CRVAL3']+(np.arange(n)+1-hd['CRPIX3'])*hd['CDELT3']
    v=2.99792458e5*(f0-nu)/f0
    cube=np.nan_to_num(h[0].data[0].astype(np.float32))
    ref=(v-V<-600)|((v-V>450)&(v-V<600)); core=abs(v-V)<win
    c=cube-(np.median(cube[ref],axis=0) if base else 0)
    dv=abs(np.median(np.diff(v))); m=c[core].sum(0)*dv
    ps=abs(hd['CDELT1'])*3600
    return m,ps,(hd['BMAJ']*3600,hd['BMIN']*3600),w.celestial
def smooth(m,ps,b0,B):
    # circular target; isotropic kernel width via quadrature using geometric mean beam (approx)
    bg=np.sqrt(b0[0]*b0[1]); s=np.sqrt(max(B**2-bg**2,0))/2.3548/ps
    return gaussian_filter(m,s)*(B*B)/(bg*bg)
for win in (250,350,500):
  for base in (True,False):
    r={}
    for ous,spw in (('Xb4d','spw25'),('Xb4f','spw16')):
        m,ps,b0,wc=m0(ous,spw,win,base); s=smooth(m,ps,b0,10.)
        for p,(ra,dec) in NUC.items():
            x,y=wc.world_to_pixel_values(ra,dec); r[(ous,p)]=s[int(round(float(y))),int(round(float(x)))]
    print(win,base,{p:round(r[('Xb4d',p)]/r[('Xb4f',p)],3) for p in NUC})
