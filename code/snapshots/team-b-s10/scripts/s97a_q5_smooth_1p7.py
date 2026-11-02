# Q5 round 2: Xb4d (0.95x0.75") Band 6 cubes smoothed per channel to a common circular 1.7" beam (Harada+2018 resolution), spectra at fixed positions (peak pixel of
# the smoothed cube, pbcor Jy/beam in the 1.7" beam) saved in the s90 npz format as spec_Xb4s_spwXX.npz. Convolution kernel per channel = deconvolution of the channel beam
# (BEAMS table) from the target beam with radio_beam; flux rescaled by the beam-area ratio (Jy/beam of the target beam, as CASA imsmooth does).
# usage: python3 -I s97a_q5_smooth_1p7.py [Xb4s|Xb4t]
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits; from astropy.wcs import WCS; import astropy.units as u
from astropy.convolution import convolve_fft; from radio_beam import Beam
from common import files
from regions import REGIONS
TAG=sys.argv[1] if len(sys.argv)>1 else 'Xb4s'   # Xb4s: circular 1.7" (Harada); Xb4t: 4.31x2.79" PA 86.7 (the Xb47 CS 7-6 beam)
TB=Beam(1.7*u.arcsec,1.7*u.arcsec,0*u.deg) if TAG=='Xb4s' else Beam(4.312*u.arcsec,2.788*u.arcsec,86.7*u.deg); REG=['N','S','Bridge','Sblue','W','E1','E2']
for sp in ('25','27','29','31'):
    fn=files(f'*Xb4d*spw{sp}.cube.I.pbcor.fits')[0]; h=fits.open(fn); hd=h[0].header; dat=h[0].data[0]; nz=dat.shape[0]
    ps=abs(hd['CDELT1'])*3600
    nu=hd['CRVAL3']+(np.arange(nz)-(hd['CRPIX3']-1))*hd['CDELT3']; wc=WCS(hd).celestial
    pos={k:[int(round(float(t))) for t in wc.world_to_pixel_values(REGIONS[k][0],REGIONS[k][1])] for k in REG}
    y0,y1,x0,x1=[int(v) for v in (min(p[1] for p in pos.values())-40,max(p[1] for p in pos.values())+40,min(p[0] for p in pos.values())-40,max(p[0] for p in pos.values())+40)]
    out={'nu':nu,'bmaj':TB.major.value,'bmin':TB.minor.value,'apr_arcsec':0.}; cache={}
    for k,(px,py) in pos.items(): out[k+'_ap']=np.full(nz,np.nan); out[k+'_pix']=np.full(nz,np.nan)
    for i in range(nz):
        im=np.nan_to_num(dat[i,y0:y1,x0:x1]); 
        if not np.any(im): continue
        b=Beam(hd['BMAJ']*u.deg,hd['BMIN']*u.deg,hd['BPA']*u.deg)  # pbcor cubes carry one (median) beam in the header
        key=(round(b.major.value,3),round(b.minor.value,3),round(b.pa.value,0))
        if key not in cache:
            try: kb=TB.deconvolve(b); ker=kb.as_kernel(ps*u.arcsec,x_size=61,y_size=61).array; fac=(TB.sr/b.sr).value; cache[key]=(ker/ker.sum(),fac)
            except Exception: cache[key]=None
        if cache[key] is None: continue
        ker,fac=cache[key]; sm=convolve_fft(im,ker,boundary='fill',allow_huge=True)*fac
        for k,(px,py) in pos.items(): v=sm[py-y0,px-x0]; out[k+'_ap'][i]=v; out[k+'_pix'][i]=v
    np.savez(f'/workspace/out/q5/spec_{TAG}_spw{sp}.npz',**out); print('spw',sp,'done',len(cache),'beams',flush=True)
