# continuum / mfs survey: peak, position, noise, beam for every 2D pbcor image
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy import ndimage
rows=[]
for fn in files('*.pbcor.fits')+files('*tapered*cont*'):
    if 'cube' in fn or 'repBW' in fn: continue
    if fn.endswith('pbcor.fits') is False: continue
    d,h=load2d(fn); w=WCS(h).celestial
    ny,nx=d.shape
    # robust noise from outer-edge-free region: MAD of whole image minus brightest? use corners (20%)
    fin=d[np.isfinite(d)]; fy,fx=np.where(np.isfinite(d))
    # noise: MAD of pixels in the outer 25% annulus of the *finite* (pb-limited) area
    rr=np.hypot(fx-nx/2,fy-ny/2); sel=fin[rr>np.percentile(rr,75)]
    rms=1.4826*np.median(np.abs(sel-np.median(sel)))
    iy,ix=np.unravel_index(np.nanargmax(d),d.shape)
    ra,dec=w.pixel_to_world_values(ix,iy)
    # local maxima >5 sigma
    mx=ndimage.maximum_filter(np.nan_to_num(d),size=max(3,int(1.5/ (abs(h['CDELT1'])*3600))))
    pk=np.argwhere((np.nan_to_num(d)==mx)&(d>5*rms))
    pk=sorted(pk,key=lambda p:-d[p[0],p[1]])[:6]
    pks=[(round(float(w.pixel_to_world_values(p[1],p[0])[0]),5),round(float(w.pixel_to_world_values(p[1],p[0])[1]),5),round(float(d[p[0],p[1]])*1e3,3)) for p in pk]
    print(short(fn),'beam %.2fx%.2f'%(h['BMAJ']*3600,h['BMIN']*3600),'nu=%.2f'%(h.get('CRVAL3',0)/1e9),'peak %.3f mJy/b @%.5f %.5f'%(d[iy,ix]*1e3,ra,dec),'rms %.3f mJy'%(rms*1e3),'SNR %.0f'%(d[iy,ix]/rms),'n>5s',len(pk),'nfinite',fin.size,'/',d.size)
    print('    peaks(ra,dec,mJy):',pks)
