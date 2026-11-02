# Q7 library: cached cubes -> baseline-subtracted moment-0, matched-beam smoothing, annulus sums with null-window noise.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
from scipy.signal import fftconvolve
from matched import kernel
C=299792.458; VSYS=2780.; RA0,DE0=156.96346,-43.90461
OM=lambda b: np.pi*b[0]*b[1]/(4*np.log(2))
BASE=((-900,-520),(520,2000))     # baseline channels (km/s from VSYS): outside HCN v2=1f (+401 in HCO+), HOC+ (-1000)
class Cube:
    def __init__(s,ous,ln):
        z=np.load(f'/workspace/out/q7/cube_{ous}_{ln}.npz'); s.cube=z['cube']; s.v=z['v']-VSYS; s.pix=float(z['pix']); s.beam=tuple(z['beam']); s.ous=ous; s.ln=ln
        s.hd={'CDELT1':float(z['cd1']),'CDELT2':float(z['cd2'])}; s.ra=z['ra']; s.de=z['de']
        s.dv=abs(np.median(np.diff(s.v))); s.ny,s.nx=s.cube.shape[1:]
        m=np.zeros(len(s.v),bool)
        for a,b in BASE: m|=(s.v>=a)&(s.v<=b)
        s.bm=m; X=np.vstack([np.ones(m.sum()),s.v[m]/500.]).T; P=np.linalg.pinv(X)   # linear baseline per pixel
        c=np.tensordot(P,s.cube[m],axes=(1,0)); s.cube=s.cube-(np.tensordot(np.vstack([np.ones(len(s.v)),s.v/500.]).T,c,axes=(1,0))).astype(np.float32)
        s.ch_rms=1.4826*np.median(np.abs(s.cube[m]-np.median(s.cube[m],axis=0)),axis=0)   # per-pixel baseline rms (Jy/beam)
        xx=(s.ra-RA0)*np.cos(np.deg2rad(DE0))*3600; yy=(s.de-DE0)*3600
        s.X,s.Y=np.meshgrid(xx,yy)
    def m0(s,lo,hi,arr=None):
        m=(s.v>=lo)&(s.v<=hi); return (s.cube if arr is None else arr)[m].sum(0)*s.dv
def smooth(img,c,tb):
    """convolve a native-beam map (Jy/beam km/s) to target beam tb=(maj,min,pa); returns map in Jy/beam(tb) km/s"""
    k=kernel(c.hd,tb,c.beam)
    if k is None: raise ValueError('target smaller than native')
    return fftconvolve(np.nan_to_num(img),k,mode='same')*OM(tb)/OM(c.beam)
def peak_near(c,img,x0,y0,rad=0.6):
    m=np.hypot(c.X-x0,c.Y-y0)<rad; j=np.argmax(np.where(m,img,-1e9)); iy,ix=np.unravel_index(j,img.shape); return c.X[iy,ix],c.Y[iy,ix]
def annulus_flux(c,img,x0,y0,r1,r2,tb):
    """integrated flux (Jy km/s) in annulus: sum(img)/beam area in pixels"""
    m=(np.hypot(c.X-x0,c.Y-y0)>=r1)&(np.hypot(c.X-x0,c.Y-y0)<r2); return img[m].sum()/(OM(tb)/c.pix**2), m.sum()
