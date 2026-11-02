# Spectra for every region in every measurement cube: peak pixel (refined within 0.6" of nominal, max of |v|<350 integrated map where line exists, else mean map),
# beam-aperture mean (diameter = geometric-mean beam FWHM), saved per cube.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
os.makedirs('/workspace/out/regspec',exist_ok=True)
def run(fn):
    h,hd,w,nu=cube_info(fn); nz=len(nu); ps=abs(hd['CDELT1'])*3600; bmaj,bmin=hd['BMAJ']*3600,hd['BMIN']*3600
    wc=WCS(hd).celestial; ny,nx=hd['NAXIS2'],hd['NAXIS1']
    rad=0.6/ps; apr=0.5*np.sqrt(bmaj*bmin)/ps
    out={'nu':nu,'bmaj':bmaj,'bmin':bmin,'ps':ps}
    for k,(ra,dec,_) in REGIONS.items():
        x,y=wc.world_to_pixel_values(ra,dec); m=int(np.ceil(rad+apr+2))
        x0,x1=max(0,int(x)-m),min(nx,int(x)+m+1); y0,y1=max(0,int(y)-m),min(ny,int(y)+m+1)
        if x1-x0<3 or y1-y0<3 or not np.isfinite(h[0].data[0,nz//2,int(y),int(x)]) : continue
        patch=np.array(h[0].data[0,:,y0:y1,x0:x1],dtype=np.float32)        # (nz,py,px)
        yy,xx=np.mgrid[y0:y1,x0:x1]
        # strongest-line map: sum of the 20% highest channels' positive flux is unstable; use channel-integrated |signal| above 0 in the central 40% brightest channels
        tot=np.nansum(patch,axis=0)
        disk=np.hypot(xx-x,yy-y)<=rad
        iy,ix=np.unravel_index(np.argmax(np.where(disk,tot,-np.inf)),tot.shape)
        pk=patch[:,iy,ix]
        ap=np.hypot(xx-xx[iy,ix],yy-yy[iy,ix])<=apr
        out[k+'_pix']=pk; out[k+'_ap']=np.nanmean(patch[:,ap],axis=1); out[k+'_xy']=wc.pixel_to_world_values(xx[iy,ix],yy[iy,ix])
    np.savez('/workspace/out/regspec/'+short(fn)+'.npz',**out); print(short(fn),[k for k in REGIONS if k+'_pix' in out],flush=True)
if __name__=='__main__':
    fns=files('*Xb4b*spw*.cube.I.pbcor.fits')+files('*Xb4d*spw*.cube.I.pbcor.fits')+files('*Xb4f*.cube.I.pbcor.fits')+files('*Xb47*.cube.I.pbcor.fits')+files('*Xb55*.cube.I.pbcor.fits')
    for fn in fns: run(fn)
