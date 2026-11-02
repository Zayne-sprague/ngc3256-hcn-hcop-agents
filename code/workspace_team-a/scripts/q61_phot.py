# Q61 step 1: matched circular-beam continuum photometry at N and S in Bands 3/5/6/7 (archive pbcor continuum images).
# python3 -I q61_phot.py -> out/q61/phot.json, phot.log.  Peak within 0.5" of the nominal position of the smoothed image (Jy/beam of the common beam).
import json,numpy as np,warnings; warnings.filterwarnings('ignore')
from astropy.io import fits; from astropy.wcs import WCS; from scipy.ndimage import map_coordinates
D='/data/member.uid___A001_X1273_'
F={'B3':('/archive/files/2016.1.00965.S/member.uid___A001_X879_X94a.NGC_3256_sci.spw29.mfs.I.pbcor.fits',98.05),
 'B5':(D+'Xb55.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits',170.7),
 'B6d':(D+'Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits',257.8),
 'B6f':(D+'Xb4f.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits',257.7),
 'B7':(D+'Xb47.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits',346.5)}
POS={'N':(156.96341,-43.90392),'S':(156.96341,-43.90533)}
def cov(maj,mn,pa):
    s1,s2=maj/2.3548,mn/2.3548; t=np.radians(pa); R=np.array([[np.sin(t),np.cos(t)],[np.cos(t),-np.sin(t)]]); return R@np.diag([s1**2,s2**2])@R.T
def load(k):
    h=fits.open(F[k][0]); hd=h[0].header; return np.squeeze(h[0].data).astype(float),hd,WCS(hd).celestial
def smooth_to(im,hd,fwhm):
    px=abs(hd['CDELT2'])*3600; bm,bn,pa=hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']
    K=np.eye(2)*(fwhm/2.3548)**2-cov(bm,bn,pa)
    if np.linalg.eigvalsh(K).min()<-1e-9: return None
    ny,nx=im.shape; KX,KY=np.meshgrid(np.fft.fftfreq(nx,px),np.fft.fftfreq(ny,px))
    arg=2*np.pi**2*(K[0,0]*KX**2-2*K[0,1]*KX*KY+K[1,1]*KY**2)
    return np.fft.ifft2(np.fft.fft2(np.nan_to_num(im))*np.exp(-arg)).real*(fwhm**2/(bm*bn))
BEAMS={1.1:['B3','B6d'],1.5:['B3','B6d'],2.5:['B3','B6d'],4.5:['B3','B6d','B7'],6.5:['B3','B6d','B6f','B7'],10.5:['B3','B5','B6d','B6f','B7']}
res={}; log=[]
IM={k:load(k) for k in F}
for fw,ks in BEAMS.items():
    for k in ks:
        im,hd,w=IM[k]; sm=smooth_to(im,hd,fw)
        if sm is None: log.append(f'{fw}" {k}: native beam larger than target'); continue
        px=abs(hd['CDELT2'])*3600
        # noise: MAD in an emission-poor annulus (12-25" from the N-S midpoint) of the smoothed image (upper bound: residual extended emission)
        ra0,de0=156.96341,-43.90462; yy,xx=np.mgrid[0:im.shape[0],0:im.shape[1]]
        x0,y0=[float(a) for a in w.world_to_pixel_values(ra0,de0)]; r=np.hypot(xx-x0,yy-y0)*px
        m=(r>12)&(r<25)&np.isfinite(im)&(im!=0)
        sig=1.4826*np.nanmedian(np.abs(sm[m]-np.nanmedian(sm[m]))) if m.sum()>50 else np.nan
        for s,(ra,de) in POS.items():
            x,y=[float(a) for a in w.world_to_pixel_values(ra,de)]
            best=(-1e9,0,0)
            for dy in np.arange(-0.5,0.51,px/2):
                for dx in np.arange(-0.5,0.51,px/2):
                    if np.hypot(dx,dy)>0.5: continue
                    v=map_coordinates(sm,[[y+dy/px],[x+dx/px]],order=1)[0]
                    if v>best[0]: best=(v,dx,dy)
            res[f'{fw}|{k}|{s}']=dict(nu=F[k][1],mJy_beam=best[0]*1e3,rms_mJy=sig*1e3,off=(best[1],best[2]))
            log.append(f'{fw:5.1f}" {k:4s} {F[k][1]:6.1f} GHz {s}: {best[0]*1e3:8.3f} mJy/beam (rms~{sig*1e3:.3f}; peak offset {best[1]:+.2f},{best[2]:+.2f}")')
json.dump(res,open('/workspace/out/q61/phot.json','w'),indent=1); open('/workspace/out/q61/phot.log','w').write('\n'.join(log)); print('\n'.join(log))
