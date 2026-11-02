# Q36 step 1: nuclear continuum brightness temperatures at MATCHED circular resolution from the archive continuum images,
# and the dust/free-free decomposition. Output out/q36/q36_cont.json + log.
import glob,json,numpy as np,warnings; warnings.filterwarnings('ignore')
from astropy.io import fits; from astropy.wcs import WCS
D='/data/member.uid___A001_X1273_'
F={'B5':('Xb55.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits',170.7),'B6f':('Xb4f.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits',257.7),
   'B7':('Xb47.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits',346.5),'B6d':('Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits',257.8),'B6b':('Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits',257.8)}
POS={'N':(156.96341,-43.90392),'S':(156.96341,-43.90533)}
def cov(maj,min_,pa):   # covariance (in arcsec^2, x=East? use sky pixel axes): sigma^2 from FWHM; PA east of north
    s1,s2=maj/2.3548,min_/2.3548; t=np.radians(pa); R=np.array([[np.sin(t),np.cos(t)],[np.cos(t),-np.sin(t)]])  # columns: major axis dir (RA,Dec) -> x=RA(east) y=Dec
    return R@np.diag([s1**2,s2**2])@R.T
def load(k):
    h=fits.open(D+F[k][0]); hd=h[0].header; im=np.squeeze(h[0].data).astype(float); w=WCS(hd).celestial
    return im,hd,w
def smooth_to(im,hd,fwhm):
    """convolve to a circular beam of FWHM `fwhm`'' (analytic Gaussian kernel = target - native); returns image in Jy/beam of the new beam"""
    px=abs(hd['CDELT2'])*3600; bm,bn,pa=hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']
    C0=cov(bm,bn,pa); Ct=np.eye(2)*(fwhm/2.3548)**2; K=Ct-C0
    ev=np.linalg.eigvalsh(K)
    if ev.min()<-1e-9: return None
    ny,nx=im.shape; ky=np.fft.fftfreq(ny,px); kx=np.fft.fftfreq(nx,px)   # cycles/arcsec
    KX,KY=np.meshgrid(kx,ky); q=np.array([KX,KY]); # x = RA axis (increasing to the east is -x in pixels: symmetric kernel so sign irrelevant)
    # FT of Gaussian with covariance K (kernel orientation: pixel x = -RA, y=Dec -> flip sign of xy term)
    Kp=np.array([[K[0,0],-K[0,1]],[-K[0,1],K[1,1]]])
    arg=2*np.pi**2*(Kp[0,0]*KX**2+2*Kp[0,1]*KX*KY+Kp[1,1]*KY**2)
    out=np.fft.ifft2(np.fft.fft2(np.nan_to_num(im))*np.exp(-arg)).real
    return out*(fwhm**2/(bm*bn))     # Jy/beam rescale: beam area ratio
def tb(Jy,nu,fwhm_maj,fwhm_min): return 1.222e6*Jy/(nu**2*fwhm_maj*fwhm_min)      # RJ brightness temp [K], nu in GHz
res={}; log=[]
for target in (10.5,6.5):
    for k in ('B5','B6f','B7'):
        im,hd,w=load(k); sm=smooth_to(im,hd,target)
        if sm is None: continue
        for r,(ra,dec) in POS.items():
            x,y=w.world_to_pixel_values(ra,dec); x,y=int(round(float(x))),int(round(float(y)))
            # peak within 1 pixel (bilinear)
            from scipy.ndimage import map_coordinates
            v=float(map_coordinates(sm,[[float(w.world_to_pixel_values(ra,dec)[1])],[float(w.world_to_pixel_values(ra,dec)[0])]],order=1)[0])
            s=1.4826*np.nanmedian(np.abs(sm[np.isfinite(im)]-np.nanmedian(sm[np.isfinite(im)])))
            res[f'{k}_{target}_{r}']=dict(nu=F[k][1],Jy_beam=v,mJy=v*1e3,Tb_K=tb(v,F[k][1],target,target),rms_mJy_rough=float(s*1e3))
            log.append(f'{target:5.1f}" {k:4s} {F[k][1]:6.1f} GHz {r}: {v*1e3:7.2f} mJy/beam  Tb_RJ={tb(v,F[k][1],target,target)*1e3:8.2f} mK')
# native-resolution compact core values
for k in ('B6b','B6d'):
    im,hd,w=load(k)
    for r,(ra,dec) in POS.items():
        x,y=w.world_to_pixel_values(ra,dec); x,y=int(round(float(x))),int(round(float(y)))
        cut=im[y-3:y+4,x-3:x+4]; v=float(np.nanmax(cut)); bm,bn=hd['BMAJ']*3600,hd['BMIN']*3600
        res[f'{k}_native_{r}']=dict(nu=F[k][1],Jy_beam=v,mJy=v*1e3,Tb_K=tb(v,F[k][1],bm,bn),beam=[bm,bn])
        log.append(f'native {k} {bm:.2f}x{bn:.2f}" {r}: peak {v*1e3:.3f} mJy/beam Tb_RJ={tb(v,F[k][1],bm,bn):.3f} K')
json.dump(res,open('/workspace/out/q36/q36_cont.json','w'),indent=1); open('/workspace/out/q36/q36_cont.log','w').write('\n'.join(log)); print('\n'.join(log))
