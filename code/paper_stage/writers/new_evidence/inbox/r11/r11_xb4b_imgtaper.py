"""r11: image-plane Gaussian uv-weighting (= convolution of the restored Xb4b 0.24" cube) of HCN/HCO+ 3-2 to circular 1.05" and 1.5" beams;
N/S ratios with the same estimators as r11_apertures.py.  Xb4b archive pbcor cubes (/data, spw27 HCN, spw25 HCO+), 8"x8" cut-outs around N/S.
Gaussian kernel FWHM = sqrt(B^2-b0^2); T_b is conserved by smoothing, so the kernel has unit area. Usage: python3 -I r11_xb4b_imgtaper.py"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve_fft, Gaussian2DKernel
from scipy.ndimage import map_coordinates
from regions import REGIONS
C=299792.458; V0=2780.; F0={'HCN32':265.88618,'HCOp32':267.55762}; SPW={'HCN32':'spw27','HCOp32':'spw25'}
D='/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.%s.cube.I.pbcor.fits'
def load(ln):
    h=fits.open(D%SPW[ln],memmap=True)[0]; H=h.header; w=WCS(H).celestial
    nu=(H['CRVAL3']+(np.arange(H['NAXIS3'])-(H['CRPIX3']-1))*H['CDELT3'])/1e9; v=C*(1-nu/F0[ln])-V0
    sel=np.where(abs(v)<1500)[0]
    cx,cy=w.world_to_pixel_values(156.9634,-43.9046); pix=abs(H['CDELT2'])*3600; hw=int(5.5/pix)
    x0,y0=int(cx)-hw,int(cy)-hw
    d=np.nan_to_num(np.array(h.data[0,sel[0]:sel[-1]+1,y0:y0+2*hw,x0:x0+2*hw],dtype=float))
    b0=np.sqrt(H['BMAJ']*H['BMIN'])*3600
    Tb=d*1.222e6/(F0[ln]**2*H['BMAJ']*3600*H['BMIN']*3600)
    return Tb,v[sel],pix,w,(x0,y0),b0,H['CDELT3']
def kw(v,dv): return np.clip(np.minimum(v+dv/2,350)-np.maximum(v-dv/2,-350),0,dv)/dv
res={}; lines=[]
for B in (1.05,1.5):
    D_={}
    for ln in F0:
        Tb,v,pix,w,o,b0,cd=load(ln); g=np.sqrt(B**2-b0**2)/pix/2.3548
        Tc=np.stack([convolve_fft(t,Gaussian2DKernel(g),boundary='fill',normalize_kernel=True,allow_huge=True) for t in Tb])
        lf=(v<-420)|(v>440); Tc=Tc-np.median(Tc[lf],axis=0); D_[ln]=(Tc,v,lf)
    Tc,v,lf=D_['HCOp32']; dv=abs(v[1]-v[0]); k=kw(v,dv); m0=np.tensordot(k,Tc,axes=(0,0))*dv; ny,nx=m0.shape; Y,X=np.mgrid[:ny,:nx]
    rows={}
    for rn in ('N','S'):
        px,py=[float(t)-o_ for t,o_ in zip(w.world_to_pixel_values(*REGIONS[rn][:2]),o)]
        mm=np.where(np.hypot(X-px,Y-py)*pix<0.4,m0,-np.inf); j=np.argmax(mm); cx,cy=float(X.flat[j]),float(Y.flat[j])
        gg=np.arange(-B*1.2,B*1.2+1e-9,0.05)/pix; gy,gx=np.meshgrid(gg,gg,indexing='ij'); rr=np.hypot(gy,gx)*pix
        for ln in F0:
            T,vl,lf_=D_[ln]; kl=kw(vl,abs(vl[1]-vl[0])); dvl=abs(vl[1]-vl[0]); cube=np.stack([map_coordinates(t,[gy+cy,gx+cx],order=3) for t in T])
            for lab,r in [('peak',None),('r0.25',.25),('r0.5',.5),('r0.75',.75),('rbeam',B)]:
                sp=T[:,int(cy),int(cx)] if r is None else cube[:,rr<=r].mean(axis=1)
                rows[(rn,ln,lab)]=((kl*sp).sum()*dvl, np.std(sp[lf_])*np.sqrt((kl**2).sum())*dvl)
    for lab in ('peak','r0.25','r0.5','r0.75','rbeam'):
        R={}
        for rn in 'NS':
            h,o_=rows[(rn,'HCN32',lab)],rows[(rn,'HCOp32',lab)]; R[rn]=(h[0]/o_[0],h[0]/o_[0]*np.hypot(h[1]/h[0],o_[1]/o_[0]))
        c=R['N'][0]/R['S'][0]; ce=c*np.hypot(R['N'][1]/R['N'][0],R['S'][1]/R['S'][0])
        res[f'{B}|{lab}']=dict(N=R['N'],S=R['S'],NS=(c,ce))
        lines.append(f'Xb4b conv {B}\" {lab:6s} N {R["N"][0]:.3f}+-{R["N"][1]:.3f} S {R["S"][0]:.3f}+-{R["S"][1]:.3f} N/S {c:.3f}+-{ce:.3f}')
open('r11_xb4b_imgtaper.txt','w').write('\n'.join(lines)); print('\n'.join(lines)); json.dump(res,open('r11_xb4b_imgtaper.json','w'))
