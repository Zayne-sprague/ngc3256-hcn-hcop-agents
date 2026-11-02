# Q30: matched-beam 12CO(2-1)/13CO(2-1) T_b at the S and N cores (2015.1.00902.S vs 2015.1.00412.S spw0)
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.convolution import convolve_fft, Gaussian2DKernel
from radio_beam import Beam
C=299792.458
P12='/archive/files/2015.1.00902.S/member.uid___A001_X2fb_X881.ari_l.ngc_3256_sci.spw2_228401MHz.12m.cube.I.pbcor.fits'
P13='/archive/files/2015.1.00412.S/member.uid___A001_X33e_X155.ari_l.ngc_3256_sci.spw0_218561MHz.12m.cube.I.pbcor.fits'
RA0,DE0=156.96346,-43.90461
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
W=pd=None
import csv
hcn={r['src']:r for r in csv.DictReader(open('/workspace/out/q5/q5e_width.csv')) if r['line']=='HCN' and r['aperture_r']=='0.0'}
def load(P,F0):
    h=fits.open(P,memmap=True); hd=h[0].header; cube=np.array(h[0].data[0],float)
    nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/F0)*C
    return hd,cube,nu,v
def conv(hd,cube,tb):
    b0=Beam(hd['BMAJ']*u.deg,hd['BMIN']*u.deg,hd['BPA']*u.deg)
    k=tb.deconvolve(b0); pix=abs(hd['CDELT2'])*u.deg
    ker=k.as_kernel(pix); fac=tb.sr.value/b0.sr.value
    out=np.empty_like(cube)
    for i in range(cube.shape[0]):
        out[i]=convolve_fft(np.nan_to_num(cube[i]),ker,boundary='fill',allow_huge=True,normalize_kernel=False)*fac
    return out
def sample(hd,cube,x,y):
    w=WCS(hd).celestial
    ra=RA0+x/3600/np.cos(np.deg2rad(DE0)); de=DE0+y/3600
    px,py=w.world_to_pixel_values(ra,de); ix,iy=int(round(float(px))),int(round(float(py)))
    return cube[:,iy,ix]
def TRJ(S,nu,bmaj,bmin): return 1.222e3*S*1e3/((nu/1e9)**2*bmaj*bmin)
res={}
for bname,bs in [('native',None),('1.0',1.0),('1.5',1.5),('2.5',2.5)]:
    out={}
    for line,P,F0 in [('12',P12,230.538e9),('13',P13,220.3987e9)]:
        hd,cube,nu,v=load(P,F0)
        if bs is None: bm=(hd['BMAJ']*3600,hd['BMIN']*3600); cc=cube
        else:
            tb=Beam(bs*u.arcsec,bs*u.arcsec,0*u.deg); cc=conv(hd,cube,tb); bm=(bs,bs)
        for c,(x,y) in CORES.items():
            sp=sample(hd,cc,x,y)
            # velocity offset set: also off-position nulls
            nulls=[sample(hd,cc,x+dx,y+dy) for dx,dy in [(0.4,0),(-0.4,0),(0,0.4),(0,-0.4)]] if bs else []
            out[(line,c)]=dict(v=v,nu=nu,sp=sp,bm=bm,nulls=nulls,dv=abs(np.median(np.diff(v))))
    # common velocity grid: 12CO channels (5.08 km/s); interpolate 13CO (1.95 MHz=2.66 km/s) by averaging
    res[bname]=out
    print(bname,'beams 12',out[('12','S')]['bm'],'13',out[('13','S')]['bm'])
import pickle; pickle.dump(res,open('/workspace/out/q30/spectra.pkl','wb'))
