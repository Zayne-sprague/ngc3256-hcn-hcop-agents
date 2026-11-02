# Q23: convolve the 0.82x0.68" CO(2-1) cube to circular beams and measure peak T_b (RJ, Planck) at the N/S core positions
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
C=299792.458; H=6.62607e-27; KB=1.380649e-16; CL=2.99792458e10
P='/archive/files/2015.1.00902.S/member.uid___A001_X2fb_X881.ari_l.ngc_3256_sci.spw2_228401MHz.12m.cube.I.pbcor.fits'
RA0,DE0=156.96346,-43.90461
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
h=fits.open(P,memmap=True); hd=h[0].header; cube=np.nan_to_num(h[0].data[0].astype(float))
nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/230.538e9)*C
bmaj,bmin,bpa=hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']; pix=abs(hd['CDELT2'])*3600
w=WCS(hd).celestial
ny,nx=cube.shape[1:]; Y,X=np.mgrid[:ny,:nx]; ra,de=w.pixel_to_world_values(X,Y)
dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600
off=np.abs(v-2780)>550
def TRJ(S,nu0,a,b): return 1.222e3*S*1e3/((nu0/1e9)**2*a*b)
def TPl(S,nu0,a,b):
    I=S/(1.1331*a*b*(np.pi/180/3600)**2)*1e-23; return H*nu0/KB/np.log1p(2*H*nu0**3/CL**2/I)
out={}
# circular target beams; the kernel is the quadrature difference (circularised using the geometric-mean beam: a small approximation, stated)
bgeo=np.sqrt(bmaj*bmin)
for tgt in [None,1.0,1.7,2.5]:
    if tgt is None: c=cube; a,b=bmaj,bmin; lab='native'
    else:
        k=np.sqrt(tgt**2-bgeo**2)/pix/2.3548
        c=np.stack([gaussian_filter(ch,k) for ch in cube])*(tgt**2/bgeo**2)   # Jy/beam scales with beam area for conserved Jy/pixel
        a=b=tgt; lab=f'{tgt}"'
    for nm,(x0,y0) in CORES.items():
        r=np.hypot(dx-x0,dy-y0); iy,ix=np.unravel_index(np.argmin(r),r.shape)
        sp=c[:,iy,ix]; sig=np.std(sp[off]); k_=np.argmax(sp)
        # integrated-intensity-weighted: also report T_b of channel-averaged 3-ch smoothed peak for noise
        out[f'{nm}_{lab}']=dict(peak_Jy=float(sp[k_]),noise_Jy=float(sig),v=float(v[k_]),Tb_RJ=float(TRJ(sp[k_],nu[k_],a,b)),Tb_Planck=float(TPl(sp[k_],nu[k_],a,b)),e_Tb_RJ=float(TRJ(sig,nu[k_],a,b)),
                                M0_Kkms=float(TRJ(np.sum(sp[~off])*0+np.sum(sp[(v>2780-550)&(v<2780+550)]),nu[k_],a,b)*abs(np.median(np.diff(v)))))
        print(lab,nm,{kk:round(vv,3) for kk,vv in out[f'{nm}_{lab}'].items()})
json.dump(out,open('/workspace/out/q23/co_conv.json','w'),indent=1)
