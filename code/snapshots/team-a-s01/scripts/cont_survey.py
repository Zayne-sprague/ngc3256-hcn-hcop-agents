# peak finding in all 2D continuum / mfs images
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy import ndimage as ndi
import glob
rows=[]
for p in sorted(glob.glob(D+'*cont*pbcor*.fits')+glob.glob(D+'*mfs.I.pbcor*.fits')):
    h=fits.open(p); hd=h[0].header; im=np.squeeze(h[0].data)
    w=WCS(hd).celestial
    # robust noise from MAD of outer region
    ok=np.isfinite(im)
    sig=1.4826*np.nanmedian(np.abs(im[ok]-np.nanmedian(im[ok])))
    mx=np.nanmax(im); iy,ix=np.unravel_index(np.nanargmax(im),im.shape)
    ra,dec=w.pixel_to_world_values(ix,iy)
    # local maxima > 5 sigma
    lm=(im==ndi.maximum_filter(np.nan_to_num(im),size=int(max(3,2.5/ (abs(hd['CDELT2'])*3600)))))&(im>5*sig)
    ys,xs=np.where(lm)
    pk=sorted([(im[y,x],*[float(v) for v in w.pixel_to_world_values(x,y)]) for y,x in zip(ys,xs)],reverse=True)[:6]
    print(os.path.basename(p).replace('member.uid___A001_X1273_',''),f'beam={hd["BMAJ"]*3600:.2f}x{hd["BMIN"]*3600:.2f} rms~{sig*1e3:.3f}mJy max={mx*1e3:.2f}mJy at {ra:.5f},{dec:.5f} npk>5s={len(pk)}')
    for a in pk[:5]: print('     ',f'{a[0]*1e3:.2f} mJy/b  {a[1]:.5f} {a[2]:.5f}  S/N~{a[0]/sig:.0f}')
