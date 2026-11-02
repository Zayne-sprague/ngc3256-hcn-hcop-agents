# Q30: CO(1-0) (2018.1.00223.S X4d, 1.92x1.68") and CO(3-2) (X858, 0.39x0.35") at the cores, matched to 2.5" and compared with CO(2-1) (same convolution as q30_13co.py)
import numpy as np, json, warnings, sys; warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/scripts')
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.convolution import convolve_fft
from radio_beam import Beam
C=299792.458; RA0,DE0=156.96346,-43.90461
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
PS={'10':('/archive/files/2018.1.00223.S/member.uid___A001_X133f_X4d.NGC3256_sci.spw29.cube.I.pbcor.fits',115.2712018e9),
    '32':('/archive/files/2018.1.00223.S/member.uid___A001_X133d_X858.NGC3256_sci.spw29.cube.I.pbcor.fits',345.7959899e9)}
BS=2.5; TB=Beam(BS*u.arcsec,BS*u.arcsec,0*u.deg)  # CO(1-0) native 1.92x1.68"; common beam 2.5"
def TRJ(S,nu,bm): return 1.222e3*S*1e3/((nu/1e9)**2*bm[0]*bm[1])
out={}
for k,(P,F0) in PS.items():
    h=fits.open(P,memmap=True); hd=h[0].header; nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/F0)*C
    w=WCS(hd).celestial; b0=Beam(hd['BMAJ']*u.deg,hd['BMIN']*u.deg,hd['BPA']*u.deg)
    ker=TB.deconvolve(b0).as_kernel(abs(hd['CDELT2'])*u.deg); fac=TB.sr.value/b0.sr.value
    sel=np.where(np.abs(v-2780)<1500)[0] if k=='10' else np.arange(len(v))
    print(k,'v range',v.min(),v.max(),'nsel',len(sel),'beam',hd['BMAJ']*3600,hd['BMIN']*3600,flush=True)
    d={'v':v[sel],'nu':nu[sel]}
    pos={}
    for c,(x,y) in CORES.items():
        for tag,(dx,dy) in {'':(0,0),'_n1':(0.4,0),'_n2':(-0.4,0),'_n3':(0,0.4),'_n4':(0,-0.4)}.items():
            ra=RA0+(x+dx)/3600/np.cos(np.deg2rad(DE0)); de=DE0+(y+dy)/3600
            px,py=w.world_to_pixel_values(ra,de); pos[c+tag]=(int(round(float(px))),int(round(float(py))))
    nat={p:[] for p in pos}; cv={p:[] for p in pos}
    for i in sel:
        im=np.nan_to_num(np.array(h[0].data[0,i],float))
        # sample native: peak within 1 native pixel box? use nearest pixel
        c1=convolve_fft(im,ker,boundary='fill',allow_huge=True,normalize_kernel=False)*fac
        for p,(ix,iy) in pos.items(): nat[p].append(im[iy,ix]); cv[p].append(c1[iy,ix])
    d['nat']={p:np.array(a) for p,a in nat.items()}; d['conv']={p:np.array(a) for p,a in cv.items()}
    d['bm']=(hd['BMAJ']*3600,hd['BMIN']*3600)
    out[k]=d
import pickle; pickle.dump(out,open('/workspace/out/q30/spectra_other_2p5.pkl','wb'))
