# Extract spectra from every pbcor cube: nuclei (peak pixel & beam aperture), cube-peak pixel, r=2.5" sums, r=8" sum, image-plane noise
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
def run(fn):
    h,hd,w,nu=cube_info(fn); data=h[0].data  # (1,nchan,ny,nx) memmap
    nz=len(nu); ny,nx=data.shape[-2:]
    pixas=abs(hd['CDELT1'])*3600
    bmaj,bmin=hd['BMAJ']*3600,hd['BMIN']*3600
    beam_pix=np.pi*bmaj*bmin/(4*np.log(2))/pixas**2  # beam area in pixels
    yy,xx=np.mgrid[:ny,:nx]
    pos={k:pix(hd,*v) for k,v in NUC.items()}
    mid=pix(hd,0.5*(NUC['N'][0]+NUC['S'][0]),0.5*(NUC['N'][1]+NUC['S'][1]))
    def circ(c,r): return np.hypot(xx-c[0],yy-c[1])<=r/pixas
    ap_beam=2*np.sqrt(bmaj*bmin)/2  # radius = geometric-mean FWHM (diameter 2*FWHM)? -> see below
    ap_beam=0.5*np.sqrt(bmaj*bmin)  # radius = FWHM/2 -> diameter = geometric-mean beam FWHM
    masks={'N_ap':circ(pos['N'],ap_beam),'S_ap':circ(pos['S'],ap_beam),
           'N_r2.5':circ(pos['N'],2.5),'S_r2.5':circ(pos['S'],2.5),'tot_r8':circ(mid,8.0)}
    sp={k:np.full(nz,np.nan) for k in ['N_pix','S_pix']+list(masks)}
    nfin={k:0 for k in masks}
    avg=np.zeros((ny,nx)); cnt=np.zeros((ny,nx))
    step=16 if nx>500 else 64
    noise_img=np.full(nz,np.nan)
    ri=np.hypot(xx-mid[0],yy-mid[1])*pixas   # arcsec from nuclei midpoint
    for i0 in range(0,nz,step):
        c=np.array(data[0,i0:i0+step],dtype=np.float32)  # (m,ny,nx)
        m=c.shape[0]
        for k in ('N','S'):
            x,y=pos[k]; sp[k+'_pix'][i0:i0+m]=c[:,int(round(y)),int(round(x))]
        for k,mk in masks.items():
            v=c[:,mk]
            if k.endswith('r2.5') or k=='tot_r8':   # flux density in Jy: sum / beam area
                sp[k][i0:i0+m]=np.nansum(v,axis=1)/beam_pix
            else: sp[k][i0:i0+m]=np.nanmean(v,axis=1)
        avg+=np.nansum(c,axis=0); cnt+=np.isfinite(c).sum(0)
        # image-plane noise: MAD over finite pixels in annulus 6-12" (>=1.5-3 bmaj) from the nuclei midpoint (off-nuclear; diagnostic only)
        for j in range(m):
            f=np.isfinite(c[j])&(ri>=max(6.0,1.5*bmaj))&(ri<=max(12.0,3*bmaj)); vals=c[j][f]
            noise_img[i0+j]=1.4826*np.median(np.abs(vals-np.median(vals))) if vals.size>20 else np.nan
    mean_img=np.where(cnt>0,avg/np.maximum(cnt,1),np.nan)
    # restrict peak search to central 10" x 10" box around the midpoint (nuclei region) and finite
    box=circ(mid,10.0)&np.isfinite(mean_img)
    pk=np.where(box,mean_img,-np.inf); iy,ix=np.unravel_index(np.argmax(pk),pk.shape)
    sp['peak_pix']=np.array(data[0,:,iy,ix],dtype=np.float32)
    ra,dec=WCS(hd).celestial.pixel_to_world_values(ix,iy)
    out=os.path.join('/workspace/out/spectra',short(fn)+'.npz')
    np.savez(out,nu=nu,noise_img=noise_img,mean_img=mean_img,peak_xy=[ix,iy],peak_radec=[ra,dec],bmaj=bmaj,bmin=bmin,bpa=hd['BPA'],pixas=pixas,
             ap_diam=2*ap_beam,pos_N=pos['N'],pos_S=pos['S'],**sp)
    print(short(fn),'nchan',nz,'peakpix',ix,iy,'(%.5f,%.5f)'%(ra,dec),'imgnoise med %.2f mJy'%(1e3*np.nanmedian(noise_img)),flush=True)
if __name__=='__main__':
    fns=sys.argv[1:] or (files('*.cube.I.pbcor.fits')+files('*.cube.I.manual.image*.fits')+files('*repBW.I.pbcor.fits')+files('*tapered*cube*')+files('*tapered.image.pbcor.fits'))
    for fn in sorted(set(fns)): 
        if os.path.exists('/workspace/out/spectra/'+short(fn)+'.npz') and not sys.argv[1:]: continue
        run(fn)
