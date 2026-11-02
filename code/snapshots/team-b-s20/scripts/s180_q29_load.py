"""Q29 helper: load Band 3 (2016.1.00965.S) pbcor cubes, cut out, convert to RJ T_b, convolve to a circular beam.
Shared by s181/s182. Conventions as s170 (flux-conserving unit-area kernel on T_b; radio LSRK velocity vs systemic 2780 km/s)."""
import sys, numpy as np, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/scripts')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
C=299792.458; V0=2780.
A='/archive/files/2016.1.00965.S/member.uid___A001_X879_%s.NGC_3256_sci.%s.cube.I.pbcor.fits'
CUBES={'X94a_spw25':A%('X94a','spw25'),'X94a_spw27':A%('X94a','spw27'),'X94e_spw25':A%('X94e','spw25')}
CEN=(156.9638,-43.9046)
def tbk(nu,a,b): return 1.222e6/(nu**2*a*b)      # Jy/beam -> K, nu GHz, beam arcsec
def cov(maj,mn,pa):
    s1,s2=maj/2.3548,mn/2.3548; pa=np.deg2rad(pa); u=np.array([-np.sin(pa),np.cos(pa)]); w=np.array([np.cos(pa),np.sin(pa)])
    return s1**2*np.outer(u,u)+s2**2*np.outer(w,w)
def kernel(Ct,Cn,pix):
    Ck=(Ct-Cn)/pix**2; ev=np.linalg.eigvalsh(Ck); assert ev.min()>-1e-9,ev
    sg=np.sqrt(np.sort(ev)[::-1]); R=int(np.ceil(5*sg[0]))+1; y,x=np.mgrid[-R:R+1,-R:R+1]
    Ci=np.linalg.inv(Ck+1e-9*np.eye(2)); q=Ci[0,0]*x*x+2*Ci[0,1]*x*y+Ci[1,1]*y*y; k=np.exp(-0.5*q); return k/k.sum()
def load(key,half=20.0,fref=None):
    path=CUBES.get(key,key)
    """returns Jy/beam cutout (nchan,ny,nx) float32, nu (GHz), pix('), beam(maj,min,pa), pixfn(ra,dec)->(x,y) in cutout"""
    h=fits.open(path,memmap=True); hd=h[0].header; d=h[0].data
    while d.ndim>3: d=d[0]
    nu=(hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3'])/1e9
    pix=abs(hd['CDELT2'])*3600; w=WCS(hd).celestial; cx,cy=[float(t) for t in w.world_to_pixel_values(*CEN)]; n=int(half/pix)
    ix,iy=int(round(cx)),int(round(cy)); y0=max(iy-n,0); x0=max(ix-n,0)
    cut=np.array(d[:,y0:iy+n+1,x0:ix+n+1],dtype=np.float32)
    bm=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
    def pxy(pos):
        a,b=w.world_to_pixel_values(*pos); return float(a)-x0,float(b)-y0
    return cut,nu,pix,bm,pxy
def to_tb_conv(cut,nu,pix,bm,theta,fnu):
    """Jy/beam -> RJ T_b [K] using native beam and frequency fnu (GHz), then convolve to circular beam theta (arcsec)"""
    Tb=cut*np.float32(tbk(fnu,bm[0],bm[1]))
    if theta is None: return Tb
    k=kernel(cov(theta,theta,0),cov(*bm),pix)
    out=np.empty_like(Tb)
    for i in range(Tb.shape[0]): out[i]=fftconvolve(np.nan_to_num(Tb[i]),k,'same')
    return out

X='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.%s.cube.I.pbcor.fits'
XB4D={'HCO+32':(X%'spw25',267.5576),'HCN32':(X%'spw27',265.8864)}
