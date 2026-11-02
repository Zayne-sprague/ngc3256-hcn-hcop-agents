# Moment maps (total, blue wing, red wing, velocity) of HCO+ 3-2 & HCN 3-2 (and 4-3, 2-1 low-res) from pbcor cubes, plus peak positions of emission regions
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import ndimage
LINES={'HCO+32':267.55763,'HCN32':265.88618,'HCO+43':356.73424,'HCN43':354.50547,'HCN21':177.26111,'HCO+21':178.37506}
JOBS=[('Xb4b.spw25','HCO+32'),('Xb4b.spw27','HCN32'),('Xb4d.spw25','HCO+32'),('Xb4d.spw27','HCN32'),('Xb4f.spw16','HCO+32'),('Xb4f.spw18','HCN32'),
      ('Xb47.spw18','HCO+43'),('Xb47.spw16','HCN43'),('Xb55.spw16','HCN21'),('Xb55.spw18','HCO+21')]
WIN=dict(core=(-350,350),blue=(-400,-150),red=(150,400))   # km/s relative to VSYS (radio)
def run(tag,line,fov=22.0):
    fn=files(f'*{tag.split(".")[0]}.ngc_3256_sci.{tag.split(".")[1]}.cube.I.pbcor.fits')[0]
    h,hd,w,nu=cube_info(fn); f0=LINES[line]*1e9; v=vrad(nu,f0)-VSYS
    pixas=abs(hd['CDELT1'])*3600; wc=WCS(hd).celestial
    cx,cy=wc.world_to_pixel_values(*[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)])
    r=int(fov/2/pixas); x0,x1=max(0,int(cx-r)),min(hd['NAXIS1'],int(cx+r)); y0,y1=max(0,int(cy-r)),min(hd['NAXIS2'],int(cy+r))
    sel=np.where(abs(v)<1200)[0]; i0,i1=sel[0],sel[-1]+1
    cube=np.array(h[0].data[0,i0:i1,y0:y1,x0:x1],dtype=np.float32); v=v[i0:i1]; dv=abs(np.median(np.diff(v)))
    off=abs(v)>650
    # per-pixel noise per channel from line-free channels of this cut (pbcor-aware)
    sig=np.nanstd(cube[off],axis=0)
    res={}
    for k,(a,b) in WIN.items():
        m=(v>=a)&(v<=b); n=m.sum()
        mom0=np.nansum(cube[m],axis=0)*dv                 # Jy/beam km/s
        err=sig*np.sqrt(n)*dv
        res[k]=(mom0,err)
    # moment 1 with 3 sigma clip per channel
    cl=np.where(cube>3*sig,cube,0)[(abs(v)<350)]; vv=v[abs(v)<350]
    mom1=np.nansum(cl*vv[:,None,None],axis=0)/np.where(np.nansum(cl,axis=0)>0,np.nansum(cl,axis=0),np.nan)
    ext=[(x0-cx)*pixas*-1,(x1-cx)*pixas*-1,(y0-cy)*pixas,(y1-cy)*pixas]   # offsets from midpoint, RA increases left -> "-" flip handled by x*-1... (E to left)
    fig,ax=plt.subplots(1,4,figsize=(18,4.6))
    for a,(k,(m0,e)) in zip(ax,res.items()):
        im=a.imshow(m0,origin='lower',extent=ext,cmap='magma'); a.set_title(f'{tag} {line} {k} [Jy/b km/s]',fontsize=8); plt.colorbar(im,ax=a)
        a.contour(np.flipud(m0/e)[::-1] if False else m0/e,levels=[5,10,20,50,100],colors='w',linewidths=.4,extent=ext,origin='lower')
    im=ax[3].imshow(mom1,origin='lower',extent=ext,cmap='RdBu_r',vmin=-250,vmax=250); ax[3].set_title('mom1 (3sig clip) km/s',fontsize=8); plt.colorbar(im,ax=ax[3])
    for a in ax:
        for k,(rn,dn) in NUC.items():
            xp,yp=wc.world_to_pixel_values(rn,dn); a.plot(-(xp-cx)*pixas,(yp-cy)*pixas,'c+',ms=10)
    plt.tight_layout(); plt.savefig(f'/workspace/out/mom_{tag}_{line}.png',dpi=55); plt.close()
    # peaks of core moment-0 (local maxima > 8 sigma, separated >= 1 beam)
    m0,e=res['core']; sn=m0/e; fp=max(3,int(abs(hd['BMAJ']*3600)/pixas))
    mx=ndimage.maximum_filter(np.nan_to_num(m0),size=2*fp+1); pk=np.argwhere((m0==mx)&(sn>8))
    pk=sorted(pk,key=lambda p:-m0[p[0],p[1]])[:8]
    print(f'## {tag} {line}: beam {hd["BMAJ"]*3600:.2f}x{hd["BMIN"]*3600:.2f}  dv={dv:.1f} km/s')
    for p in pk:
        ra,dec=wc.pixel_to_world_values(p[1]+x0,p[0]+y0); 
        print('   peak core mom0 %.3f Jy/b km/s (S/N %.0f) at %.5f %.5f ; blue %.3f red %.3f'%(m0[p[0],p[1]],sn[p[0],p[1]],ra,dec,res['blue'][0][p[0],p[1]],res['red'][0][p[0],p[1]]))
    np.savez(f'/workspace/out/mom_{tag}_{line}.npz',core=res['core'][0],core_err=res['core'][1],blue=res['blue'][0],blue_err=res['blue'][1],red=res['red'][0],red_err=res['red'][1],mom1=mom1,x0=x0,y0=y0)
if __name__=='__main__':
    for tag,line in JOBS: run(tag,line)
