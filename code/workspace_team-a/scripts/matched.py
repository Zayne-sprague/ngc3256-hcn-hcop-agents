# Convolve cube (pbcor, Jy/beam) to a target beam in the image plane (elliptical Gaussian kernel from covariance difference), return spectrum at (ra,dec).
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
FW=2*np.sqrt(2*np.log(2))
def cov(maj,mn,pa):
    p=np.deg2rad(pa); u=np.array([np.sin(p),np.cos(p)]); w=np.array([np.cos(p),-np.sin(p)])
    return (maj/FW)**2*np.outer(u,u)+(mn/FW)**2*np.outer(w,w)
def kernel(hd,tb,nb):
    """tb,nb=(maj,min,pa) arcsec,arcsec,deg. returns normalised kernel or None if target smaller than native"""
    S=cov(*tb)-cov(*nb)
    ev=np.linalg.eigvalsh(S)
    if ev.min()<-1e-3*ev.max(): return None
    S=S+1e-6*np.eye(2)*0
    pix=abs(hd['CDELT2'])*3600; n=int(np.ceil(5*np.sqrt(max(ev.max(),1e-6))/pix))+1
    j,i=np.mgrid[-n:n+1,-n:n+1]; E=i*hd['CDELT1']*3600; N=j*hd['CDELT2']*3600
    X=np.stack([E,N],-1); Si=np.linalg.pinv(S)
    k=np.exp(-0.5*np.einsum('...i,ij,...j->...',X,Si,X)); return k/k.sum()
def smoothed_spectrum(path,tb,ra,dec,chan=None,radius_pix=None):
    h=fits.open(path,memmap=True); hd=h[0].header; cube=h[0].data[0] if h[0].data.ndim==4 else h[0].data
    nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
    x,y=[int(round(float(a))) for a in WCS(hd).celestial.world_to_pixel_values(ra,dec)]
    k=kernel(hd,tb,nb)
    if k is None: raise ValueError('target beam smaller than native in some direction')
    n=k.shape[0]//2; pad=n+1
    sl=slice(None) if chan is None else slice(*chan)
    sub=np.array(cube[sl,max(0,y-pad):y+pad+1,max(0,x-pad):x+pad+1],float)
    sub=np.nan_to_num(sub)
    om=lambda b: np.pi*b[0]*b[1]/(4*np.log(2))
    out=np.array([fftconvolve(s,k,mode='same') for s in sub])
    yc,xc=y-max(0,y-pad),x-max(0,x-pad)
    return out[:,yc,xc]*om(tb)/om(nb), sub[:,yc,xc], nb
