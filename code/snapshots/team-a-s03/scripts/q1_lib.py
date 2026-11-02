# Q1 library: moment-0 maps from pbcor cubes (radio velocity, LSRK), smoothing to a common circular beam, empirical noise from line-free windows.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from matched import kernel
from scipy.signal import fftconvolve
RA0,DE0=156.96346,-43.90461
VSYS=2780.; HW=350.   # window centre / half width (km/s): +-350 keeps HCN v2=1f (+401 km/s redshifted from HCO+ 3-2; measured effect on the +-500 window ratio <=0.5%, q1_velratio.py) out of the HCO+ windows
OM=lambda b: np.pi*b[0]*b[1]/(4*np.log(2))
# (OUS, spw, rest GHz) of the lines used
LINES={
 'HCN 2-1':('Xb55','spw16',177.261111),'HCO+ 2-1':('Xb55','spw18',178.375056),
 'HCN 4-3':('Xb47','spw16',354.505478),'HCO+ 4-3':('Xb47','spw18',356.734223),'CS 7-6':('Xb47','spw22',342.883000),
 'HCN 3-2@Xb4f':('Xb4f','spw18',265.886434),'HCO+ 3-2@Xb4f':('Xb4f','spw16',267.557626),'SO 6_5-5_4@Xb4f':('Xb4f','spw22',251.825770),
 'HCN 3-2@Xb4d':('Xb4d','spw27',265.886434),'HCO+ 3-2@Xb4d':('Xb4d','spw25',267.557626),'SO 6_5-5_4@Xb4d':('Xb4d','spw31',251.825770),
 'HCN 3-2@Xb4b':('Xb4b','spw27',265.886434),'HCO+ 3-2@Xb4b':('Xb4b','spw25',267.557626),'SO 6_5-5_4@Xb4b':('Xb4b','spw31',251.825770)}
def offsets(hd,ny,nx,x0=0,y0=0):
    w=WCS(hd).celestial; Y,X=np.mgrid[:ny,:nx]; ra,de=w.pixel_to_world_values(X+x0,Y+y0)
    return (ra-RA0)*np.cos(np.deg2rad(DE0))*3600,(de-DE0)*3600
def load(key,crop=None,nullspec=True,window=None,gap=None):
    """return dict: m0 (Jy/beam km/s, native), nulls (list of null-window m0 maps), beam, pix, dx,dy, hd.
    crop: half-size in arcsec of the central region kept (needed for the 1.3 GB Xb4b cubes)."""
    ous,spw,f0=LINES[key]; h=fits.open(f(ous,spw),memmap=True); hd=h[0].header
    fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C;   # radio velocity (km/s) w.r.t. rest frequency f0
    dv=abs(np.median(np.diff(v))); pix=abs(hd['CDELT2'])*3600
    ny,nx=hd['NAXIS2'],hd['NAXIS1']; cx,cy=[float(a) for a in WCS(hd).celestial.world_to_pixel_values(RA0,DE0)]
    if crop: hp=int(crop/pix); x0=max(0,int(cx)-hp); x1=min(nx,int(cx)+hp+1); y0=max(0,int(cy)-hp); y1=min(ny,int(cy)+hp+1)
    else: x0,x1,y0,y1=0,nx,0,ny
    def win(lo,hi):
        m=np.where((v>=lo)&(v<=hi))[0]
        a=np.array(h[0].data[0,m[0]:m[-1]+1,y0:y1,x0:x1],float); return np.nansum(a,axis=0)*dv,len(m)
    wlo,whi=window if window else (VSYS-HW,VSYS+HW); hw=(whi-wlo)/2.
    m0,n=win(wlo,whi)
    # null windows: half the line-window width (350 km/s, rescaled by sqrt(2) to the 700 km/s window), >= 600 km/s from VSYS, >= 60 km/s inside the band edge
    nulls=[]; g=gap if gap else (450. if LINES[key][0]=='Xb47' else 600.)   # Band 7 spws are narrower in velocity: smaller guard
    if nullspec:
        for hh in (hw,hw/2.,hw/4.):   # shrink the null window if the band has too little line-free velocity range; noise rescaled by sqrt(width ratio)
            for lo in np.arange(v.min()+60,v.max()-hh-60,hh/2.):
                if lo+hh<VSYS-g or lo>VSYS+g: nulls.append((lo,win(lo,lo+hh)[0]*np.sqrt(2.*hw/hh)))
            if len(nulls)>=2: break
    dx,dy=offsets(hd,y1-y0,x1-x0,x0,y0)
    return dict(x0=x0,y0=y0,m0=m0,nulls=nulls,beam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),pix=pix,dx=dx,dy=dy,hd=hd,f0=f0,nch=n,dv=dv,vr=(v.min(),v.max()))
def smooth(d,fw,arr=None):
    """convolve to a circular beam of FWHM fw (arcsec); result in Jy/beam(new) km/s. Returns None if native beam larger."""
    nb=d['beam']; tb=(fw,fw,0.)
    k=kernel(d['hd'],tb,nb)
    if k is None: return None
    a=np.nan_to_num(d['m0'] if arr is None else arr)
    return fftconvolve(a,k,mode='same')*OM(tb)/OM(nb)
def noise_profile(d,fw,rbins=(0,4,8,12,16,20,26)):
    """empirical sigma (Jy/beam km/s) of smoothed null maps vs radius; robust (1.4826 MAD) per annulus, over all null windows"""
    r=np.hypot(d['dx'],d['dy']); out=[]
    S=[smooth(d,fw,n) if fw else n for _,n in d['nulls']]
    for a,b in zip(rbins[:-1],rbins[1:]):
        m=(r>=a)&(r<b); sg=[];md=[];nn=0
        for s in S:   # per-window robust sigma; the quietest quartile of windows is the noise (other windows contain real lines, e.g. HOC+ at -1000 km/s of HCO+)
            vals=s[m]; vals=vals[np.isfinite(vals)&(vals!=0)]
            if len(vals)>20: sg.append(1.4826*np.median(np.abs(vals-np.median(vals)))); md.append(np.median(vals)); nn=len(vals)
        out.append((a,b,np.percentile(sg,25) if sg else np.nan,np.median(md) if md else np.nan,nn))
    return out
def sigma_at(prof,r):
    for a,b,s,_,_ in prof:
        if a<=r<b and np.isfinite(s): return s
    ss=[s for *_,s,_,_ in [(0,0,p[2],0,0) for p in prof] if np.isfinite(s)]
    return ss[-1] if ss else np.nan
def sample(img,d,x,y):
    """bilinear sample at offset (dRA,dDec) arcsec (exact WCS pixel position, cropped-array offset removed)"""
    from scipy.ndimage import map_coordinates
    w=WCS(d['hd']).celestial
    px,py=[float(a) for a in w.world_to_pixel_values(RA0+x/3600/np.cos(np.deg2rad(DE0)),DE0+y/3600)]
    return float(map_coordinates(img,[[py-d['y0']],[px-d['x0']]],order=1,mode='nearest')[0])
def Tb_factor(freq_ghz,fw):
    """K per (Jy/beam) for a circular Gaussian beam of FWHM fw arcsec, Rayleigh-Jeans: 1.222e3/(nu^2 theta^2)"""
    return 1.222e6/(freq_ghz**2*fw*fw)
