# Q4 library: region spectra (Jy/beam) from pbcor cubes convolved channel-wise to a common (elliptical) target beam; region positions from out/regions.csv.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
from scipy.ndimage import map_coordinates
from common import f,freqaxis,C
from matched import kernel
OUT='/workspace/out/q4/'
OM=lambda b: np.pi*b[0]*b[1]/(4*np.log(2))
REG=pd.read_csv('/workspace/out/regions.csv').set_index('region')
def kb(freq_ghz,bmaj,bmin): return 1.222e6/(freq_ghz**2*bmaj*bmin)     # K per Jy/beam (RJ)
def spectrum(ous,spw,region,tb=None,kind='cube.I.pbcor',sep='.',pad_as=14.,pos=None):
    """spectrum (freq GHz, Jy/beam in the target beam, native beam, null rms estimate) at a region position (bilinear). tb=(maj,min,pa) arcsec/deg or None (native)."""
    h=fits.open(f(ous,spw,kind,sep),memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9
    nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); w=WCS(hd).celestial
    ra,de=(REG.loc[region,'RA_deg'],REG.loc[region,'Dec_deg']) if pos is None else pos
    px,py=[float(a) for a in w.world_to_pixel_values(ra,de)]
    pix=abs(hd['CDELT2'])*3600
    k=kernel(hd,tb,nb) if tb is not None else None
    if tb is not None and k is None: raise ValueError('target beam smaller than native')
    n=(k.shape[0]//2+2) if k is not None else 2
    x0=max(0,int(px)-n); y0=max(0,int(py)-n); x1=int(px)+n+1; y1=int(py)+n+1
    a=np.nan_to_num(np.array(h[0].data[0,:,y0:y1,x0:x1] if h[0].data.ndim==4 else h[0].data[:,y0:y1,x0:x1],np.float32))
    if k is not None: sc=OM(tb)/OM(nb); a=np.array([fftconvolve(c,k.astype(np.float32),mode='same')*sc for c in a],np.float32)
    s=np.array([map_coordinates(c,[[py-y0],[px-x0]],order=1,mode='nearest')[0] for c in a])
    return fr,s,(tb if tb is not None else nb),nb
def vel(fr,f0): return (1-fr/f0)*C
def robust_sigma(x): x=x[np.isfinite(x)]; return 1.4826*np.median(np.abs(x-np.median(x)))

def fit_target(v,s,v0,fw,off,extras=(),tmin=40,tmax=300,vmax=950):
    """joint linear fit near a strong line (centre v0, width fw): baseline(1)+A_w*exp(-|v-v0|/tau) [core |v-v0|<1.5fw masked] + Gaussian at v0+off (+ Gaussians at v0+extras, same sigma=fw/2.355). Returns amplitude-integrals (mJy/beam km/s or input units) with errors (stat+tau spread)."""
    m=(np.abs(v-v0)>1.5*fw)&(np.abs(v-v0)<vmax); x=v[m]-v0; y=s[m]; sig=fw/2.355; best=None
    for tau in np.linspace(tmin,tmax,27):
        E=np.exp(-np.abs(x)/tau); G=np.exp(-0.5*((x-off)/sig)**2)
        for it in range(3):
            A=np.vstack([np.ones_like(x),x,E,G]+[np.exp(-0.5*((x-e)/sig)**2) for e in extras]).T
            if it==0: w=np.ones_like(y)
            c,*_=np.linalg.lstsq(A*w[:,None],y*w,rcond=None); r=y-A@c; sd=robust_sigma(r); w=np.ones_like(y)/sd
        chi=np.sum((r/sd)**2); cov=np.linalg.inv((A*w[:,None]).T@(A*w[:,None]))
        if best is None or chi<best[0]: best=(chi,tau,c,cov,sd,A,w)
        res=(chi,tau,c,cov)
    chi,tau,c,cov,sd,A,w=best
    # tau systematics: spread of the vib amplitude among tau values with delta chi2<4
    amps=[];
    for t2 in np.linspace(tmin,tmax,27):
        E=np.exp(-np.abs(x)/t2); G=np.exp(-0.5*((x-off)/sig)**2); A2=np.vstack([np.ones_like(x),x,E,G]+[np.exp(-0.5*((x-e)/sig)**2) for e in extras]).T
        c2,*_=np.linalg.lstsq(A2*w[:,None],y*w,rcond=None); chi2=np.sum(((y-A2@c2)*w)**2)
        if chi2-chi<4: amps.append(c2[3])
    dv=abs(np.median(np.diff(v))); norm=sig*np.sqrt(2*np.pi)/dv   # channels
    Iv=c[3]*norm*dv; ev=np.hypot(np.sqrt(cov[3,3]),np.std(amps) if len(amps)>1 else 0)*norm*dv
    # wing level in the vib window predicted by the exponential + linear terms alone (mJy/beam km/s), width 2 FWHM
    xx=np.arange(off-fw,off+fw,dv); Iw=np.sum(c[0]+c[1]*xx+c[2]*np.exp(-np.abs(xx)/tau))*dv; ew=0.
    return Iv,ev,Iw,ew,len(x),tau
