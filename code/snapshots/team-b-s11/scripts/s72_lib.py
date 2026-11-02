# Q11 round-2 helpers: cube loader (all channels or v window, optional clean-mask cube), circular-beam convolution, aperture spectra, random-aperture noise.
import sys, os; sys.path.insert(0,'/workspace/scripts')
import numpy as np
import s30_matchedJ_ratios as S
from s30_matchedJ_ratios import REGIONS, F0, CEN, conv, VSYS, C
from common import files, cube_info, vrad
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import binary_dilation
KT=1.0126
HALF=24.0
def load_cube(ous,spw,f0,vmin=-1000,vmax=700,kind='pbcor',half=HALF):
    fn=files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.{kind}.fits')[0]; h_,hd,w,nu=cube_info(fn)
    v=vrad(nu,f0)-VSYS; ps=abs(hd['CDELT1'])*3600; wc=WCS(hd).celestial
    cx,cy=wc.world_to_pixel_values(*CEN); r=int(min(half,(min(hd["NAXIS1"],hd["NAXIS2"])//2-2)*ps)/ps); x0,y0=int(round(cx-r)),int(round(cy-r)); n=2*r
    ok=np.where((v>vmin)&(v<vmax))[0]; i0,i1=ok[0],ok[-1]+1
    cube=np.nan_to_num(np.array(h_[0].data[0,i0:i1,y0:y0+n,x0:x0+n],dtype=np.float32))
    beam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']) if 'BMAJ' in hd else None
    return cube,v[i0:i1],ps,wc,x0,y0,beam
def pxy_fn(wc,x0,y0):
    return lambda ra,dec:(wc.world_to_pixel_values(ra,dec)[0]-x0,wc.world_to_pixel_values(ra,dec)[1]-y0)
def sout_pos(ofs=(-0.2,-1.85)):
    r0=REGIONS['S'][:2]; return (r0[0]+ofs[0]/3600/np.cos(np.deg2rad(r0[1])), r0[1]+ofs[1]/3600)
def offset_pos(base,dx,dy):
    """position dx (arcsec W-E, +E) dy (+N) from base=(ra,dec)"""
    return (base[0]+dx/3600/np.cos(np.deg2rad(base[1])), base[1]+dy/3600)
def apmask(n,ps,xy,rad,rin=0.0):
    Y,X=np.mgrid[:n,:n]; d=np.hypot(X-xy[0],Y-xy[1])*ps; return (d<=rad)&(d>=rin)
def aspec(cube,ba,m): return cube[:,m].sum(1)/ba
def wsum(v,sp,a,b,dv): k=(v>=a)&(v<=b); return sp[k].sum()*dv
def random_pos(cube,v,ps,pxy,n,rad,rmin,rmax,rng,excl=()):
    nn=cube.shape[1]; dv=abs(np.median(np.diff(v)))
    m0=cube[np.abs(v)<350].sum(0)*dv; sg=np.median(np.abs(m0-np.median(m0)))*1.4826
    bright=binary_dilation(m0>4*sg,iterations=int(1.5/ps)); cx,cy=pxy(*CEN); out=[]
    while len(out)<n:
        x,y=rng.uniform(0,nn,2); d=np.hypot(x-cx,y-cy)*ps
        if rmin<d<rmax and not bright[int(y),int(x)] and all(np.hypot(x-e[0],y-e[1])*ps>2*rad for e in excl): out.append((x,y))
    return out
def prep_pair(ous,fwhm=1.05,spw=None,vmin=-1000,vmax=700):
    """HCN/HCO+ 3-2 cubes, pbcor, convolved to circular fwhm, per-pixel baseline = median(-900<v<-450). returns dict mol->R"""
    spw=spw or {'Xb4d':('spw27','spw25'),'Xb4b':('spw27','spw25')}[ous]; o={}
    for mol,sp in zip(('HCN','HCO+'),spw):
        c,v,ps,wc,x0,y0,b0=load_cube(ous,sp,F0[mol][3],vmin,vmax)
        c=conv(c,ps,b0,(fwhm,fwhm,0.0)) if fwhm else c
        c=c-np.median(c[(v>-900)&(v<-450)],axis=0); b=(fwhm,fwhm,0.0) if fwhm else b0
        o[mol]=dict(c=c,v=v,dv=abs(np.median(np.diff(v))),ps=ps,n=c.shape[1],ba=1.1331*b[0]*b[1]/ps**2,pxy=pxy_fn(wc,x0,y0),beam=b)
    return o
def ap(R,pos,rad=1.0): return aspec(R['c'],R['ba'],apmask(R['n'],R['ps'],R['pxy'](*pos),rad))
def W(R,sp,a,b): return wsum(R['v'],sp,a,b,R['dv'])
