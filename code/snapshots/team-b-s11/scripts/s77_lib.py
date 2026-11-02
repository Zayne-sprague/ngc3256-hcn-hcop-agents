# Q2 round 2 helpers: loaders for the standard, hand-made tapered (Xb4d) and manual (Xb4b) cubes onto one structure, all circularly convolved to a common beam.
import sys; sys.path.insert(0,'/workspace/scripts')
from s72_lib import *
import numpy as np
S0=REGIONS['S'][:2]; N0=REGIONS['N'][:2]
LOBE=offset_pos(S0,0.15,-1.85); RED=offset_pos(S0,0.15,3.3)
REST={'HCN':F0['HCN'][3],'HCO+':F0['HCO+'][3]}
SPW={'HCN':'spw27','HCO+':'spw25'}
def load_any(prod,mol,fwhm=1.05,half=12.0):
    """prod: std_Xb4d, std_Xb4b, tap_Xb4d, man_Xb4b (pbcor) or manflat_Xb4b. Returns dict(c,v,dv,ps,n,ba,pxy,edge) with c in Jy/beam of the circular beam (flux-rescaled), baseline NOT subtracted."""
    f0=REST[mol]; spw=SPW[mol]
    if prod.startswith('std'):
        c,v,ps,wc,x0,y0,b0=load_cube(prod[4:],spw,f0,-1000,700,half=half)
    else:
        if prod=='tap_Xb4d': fn=files(f'*Xb4d_ngc_3256_sci.{spw}.cube.I.tapered.image.pbcor.fits')[0]
        elif prod=='man_Xb4b': fn=files('*Xb4b*spw25.cube.I.manual.image.pbcor.fits')[0]
        elif prod=='manflat_Xb4b': fn=files('*Xb4b*spw25.cube.I.manual.image.fits')[0]
        h_,hd,w,nu=cube_info(fn); v=vrad(nu,f0)-VSYS; ps=abs(hd['CDELT1'])*3600; wc=WCS(hd).celestial
        cx,cy=wc.world_to_pixel_values(*CEN); r=int(min(half/ps,0.5*hd['NAXIS1']-abs(cx-0.5*hd['NAXIS1'])-2,0.5*hd['NAXIS2']-abs(cy-0.5*hd['NAXIS2'])-2)); x0,y0=int(round(cx-r)),int(round(cy-r))
        c=np.nan_to_num(np.array(h_[0].data[0,:,y0:y0+2*r,x0:x0+2*r],dtype=np.float32)); b0=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
    if fwhm: c=conv(c,ps,b0,(fwhm,fwhm,0.0)); b=(fwhm,fwhm,0.0)
    else: b=b0
    dv=abs(np.median(np.diff(v)))
    return dict(c=c,v=v,dv=dv,ps=ps,n=c.shape[1],ba=1.1331*b[0]*b[1]/ps**2,pxy=pxy_fn(wc,x0,y0),beam=b)
def field_pixels(R,rmin=3.5,rmax=None):
    n=R['n']; ps=R['ps']; m0=R['c'][np.abs(R['v'])<350].sum(0)*R['dv']; sg=np.median(np.abs(m0-np.median(m0)))*1.4826
    br=binary_dilation(m0>4*sg,iterations=int(1.5/ps)); Y,X=np.mgrid[:n,:n]; cx,cy=R['pxy'](*CEN); rc=np.hypot(X-cx,Y-cy)*ps
    return (~br)&(rc>rmin)&(rc<(rmax or n*ps*0.47))
def edge_window(v,dv,a,b):
    """channels whose centre lies in [a,b]; returns mask and the actual velocity edges"""
    k=(v>=a)&(v<=b); return k,(v[k].min()-dv/2,v[k].max()+dv/2)
