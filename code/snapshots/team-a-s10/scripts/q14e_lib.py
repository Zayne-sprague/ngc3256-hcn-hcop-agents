# Q14e library: smoothed cube cutouts (K, RJ) in a common beam, region sampling, matched moment-0 estimator with local baseline and empirical null.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
from scipy.ndimage import map_coordinates
from scipy.optimize import curve_fit
from common import f,freqaxis,C
from matched import kernel
from q4_lib import REG,kb,OM,vel,robust_sigma
CEN=(156.9635,-43.9047)   # N-S midpoint
def cut(ous,spw,tb=None,half=22.,cen=CEN,kind='cube.I.pbcor'):
    """cutout (+-half arcsec about cen) of cube smoothed channel-wise to tb=(maj,min,pa) -> returns freq GHz, cube [K] (nchan,ny,nx), wcs(cutout), pix arcsec, beam used, native beam"""
    h=fits.open(f(ous,spw,kind),memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9
    nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); w=WCS(hd).celestial; pix=abs(hd['CDELT2'])*3600
    px,py=[float(a) for a in w.world_to_pixel_values(*cen)]; n=int(half/pix)
    k=kernel(hd,tb,nb) if tb is not None else None
    if tb is not None and k is None: raise ValueError('target beam smaller than native')
    pad=(k.shape[0]//2+2) if k is not None else 0
    x0=max(0,int(px)-n-pad); x1=int(px)+n+pad+1; y0=max(0,int(py)-n-pad); y1=int(py)+n+pad+1
    d=h[0].data; d=d[0] if d.ndim==4 else d
    a=np.nan_to_num(np.array(d[:,y0:y1,x0:x1],np.float32))
    beam=tb if tb is not None else nb
    if k is not None:
        sc=OM(tb)/OM(nb); a=np.array([fftconvolve(c,k.astype(np.float32),mode='same')*sc for c in a],np.float32)
        a=a[:,pad:-pad,pad:-pad]; x0+=pad; y0+=pad
    kk=kb(fr.mean(),beam[0],beam[1]); a=a*kk
    cw=w.slice((slice(y0,y0+a.shape[1]),slice(x0,x0+a.shape[2])))
    return fr,a,cw,pix,beam,nb
def sample(a,cw,ra,dec):
    px,py=[float(z) for z in cw.world_to_pixel_values(ra,dec)]
    return np.array([map_coordinates(c,[[py],[px]],order=1,mode='nearest')[0] for c in a])
def gauss(v,A,v0,fw,a,b): return A/(fw/2.3548*2.5066)*np.exp(-0.5*((v-v0)/(fw/2.3548))**2)+a+b*(v-v0)/100.
def hco_fit(v,y,v00=2780.,fw0=150.):
    m=np.abs(v-v00)<700
    p,cv=curve_fit(gauss,v[m],y[m],p0=[y[m].max()*fw0*1.06*0.9,v00,fw0,0,0],maxfev=20000); return p
def m0(v,y,c,W,excl=None,bw=2.0):
    """moment-0 (K km/s) in window c+-W/2 with a linear baseline fitted (3-sigma clipped) to channels in c+-(W/2+bw*W) outside the window and outside excl=(lo,hi) [real line]"""
    dv=abs(np.median(np.diff(v))); inw=np.abs(v-c)<W/2
    side=(np.abs(v-c)<W/2+bw*W)&~inw
    if excl is not None: side&=~((v>excl[0])&(v<excl[1]))
    if side.sum()<8: return np.nan,np.nan
    use=side.copy(); 
    for _ in range(4):
        p=np.polyfit(v[use]-c,y[use],1); r=y-np.polyval(p,v-c); sg=robust_sigma(r[use]); use=side&(np.abs(r)<3*sg)
    return float(np.sum(r[inw])*dv),sg*dv*np.sqrt(inw.sum())
def m0map(v,a,c,W,excl=None,bw=2.0):
    """per-pixel moment-0 map (units of a x km/s) in window c+-W/2 with a per-pixel linear baseline from the side bands (same as m0())"""
    dv=abs(np.median(np.diff(v))); inw=np.abs(v-c)<W/2; side=(np.abs(v-c)<W/2+bw*W)&~inw
    if excl is not None: side&=~((v>excl[0])&(v<excl[1]))
    x=v[side]-c; A=np.vstack([np.ones_like(x),x]).T; P=np.linalg.pinv(A)       # (2,ns)
    ys=a[side].reshape(side.sum(),-1); co=P@ys; base=(co[0][None,:]+np.outer(v[inw]-c,co[1]))
    mp=(a[inw].reshape(inw.sum(),-1)-base).sum(0)*dv
    return mp.reshape(a.shape[1:])
def pos_pix(cw,ra,dec): return [float(z) for z in cw.world_to_pixel_values(ra,dec)]
def val(mp,cw,ra,dec):
    px,py=pos_pix(cw,ra,dec); return float(map_coordinates(mp,[[py],[px]],order=1,mode='nearest')[0])
