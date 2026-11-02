import glob, os, warnings, numpy as np
warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
D='/data/'
C=299792.458
VSYS=2780.0   # km/s radio-LSRK systemic adopted from dense-gas line centroids (s3b_systemic.py: 2780+-5); nominal task value 2750
def vrad(nu_obs,nu0): return C*(1-np.asarray(nu_obs)/nu0)          # radio-definition LSRK velocity [km/s] (what the cubes use)
def rest_of(nu_obs,v=None): return nu_obs/(1-(VSYS if v is None else v)/C)  # rest freq assuming systemic radio velocity v
ZSYS=1/(1-VSYS/C)-1   # equivalent redshift for radio-def 2750 km/s (0.00926); note optical z=0.0094 <-> v_radio=2792 km/s
OUS={'Xb55':'B5 12m?','Xb47':'B7','Xb4f':'B6 low-res','Xb4d':'B6 ~1"','Xb4b':'B6 0.24"'}
def files(pat): return sorted(glob.glob(D+pat))
def load2d(fn):
    h=fits.open(fn,memmap=True); d=np.squeeze(h[0].data); return d,h[0].header
def short(fn): 
    b=os.path.basename(fn).replace('member.uid___A001_X1273_','')
    return b.replace('.ngc_3256_sci','').replace('.I','').replace('.fits','')

NUC={'N':(156.96346,-43.90388),'S':(156.96340,-43.90533)}  # from s1_continuum (0.24" mfs peaks)
def cube_info(fn):
    h=fits.open(fn,memmap=True)
    hd=h[0].header; w=WCS(hd)
    n=hd['NAXIS3']; nu=(hd['CRVAL3']+(np.arange(n)+1-hd['CRPIX3'])*hd['CDELT3'])
    return h,hd,w,nu
def pix(hd,ra,dec):
    w=WCS(hd).celestial; x,y=w.world_to_pixel_values(ra,dec); return float(x),float(y)
