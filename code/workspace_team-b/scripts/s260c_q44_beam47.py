"""Q44(c): HCN/HCO+ 1-0, 3-2, 4-3 at ONE elliptical beam = the Xb47 4-3 native beam (4.14x2.70", PA 87), +-150 km/s box, aperture-mean r=0.6" (~1 pixel, 'beam-peak') and r=2.2" at the peaks
of HCO+ 1-0 (X961 2.2" map, within 0.4" of nominal; N,S,Bridge).  1-0: X961 spw0 (2015.1.00993, short-spacing complete) convolved image-plane; 3-2: Q34 12m+7m 'comb' 2.2" cubes (s191) convolved
(the 7m-alone Xb4f beam 5.7x3.5" is larger than the target, so uv-matched 12m+7m at 2.2" then smoothed); 4-3: Xb47 spw16 (HCN) and spw18 (HCO+) native (beam of the cube header, no convolution).
RJ T_b (K km/s), x(1+z) rest-frame as s193.  Baselines: linear fit to channels |v|>450 km/s from the line (and away from the other lines). Errors: line-free rms of the aperture-mean spectrum x box width.
usage: python3 -I s260c_q44_beam47.py -> out/q44/s260c.json, s260c.txt"""
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import *
from s182_q29_measure import REG,F0,vel,kern_w
from astropy.io import fits
from astropy.wcs import WCS
from regions import REGIONS
D='/archive/files/2015.1.00993.S/'
X961=D+'member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits'
B47='/data/member.uid___A001_X1273_Xb47.ngc_3256_sci.%s.cube.I.pbcor.fits'
TGT=(4.14,2.70,87.0); Z=0.00927; K=1+Z
log=[]
def P(s=''): print(s,flush=True); log.append(s)
# ---- 1-0 at the target beam
cut,nu,pix,bm,pxy=load(X961,half=22)
vO=vel(nu,88.6318*0+89.1885); vH=vel(nu,88.6318)
Tb2=to_tb_conv(cut,nu,pix,bm,2.2,np.mean(nu))              # 2.2" circular, for peak finding (as s185)
Ct=cov(*TGT); Tb10=np.empty_like(Tb2); kk=kernel(Ct,cov(*bm),pix); Tf=cut*np.float32(tbk(np.mean(nu),bm[0],bm[1]))
from scipy.signal import fftconvolve
for i in range(Tf.shape[0]): Tb10[i]=fftconvolve(np.nan_to_num(Tf[i]),kk,'same')
ny,nx=Tb2.shape[1:]; Y,X=np.mgrid[:ny,:nx]; dv=abs(vO[1]-vO[0])
m0=np.tensordot(kern_w(vO,dv,-150,150),Tb2,axes=(0,0))*dv
def peak(rn):
    px,py=pxy(REG[rn]); mm=np.where(np.hypot(X-px,Y-py)*pix<0.4,m0,-np.inf); j=np.argmax(mm)
    mp=(float(X.flat[j]),float(Y.flat[j])); dxs=(mp[0]-px)*pix; dys=(mp[1]-py)*pix
    return (REG[rn][0]-dxs/3600/np.cos(np.deg2rad(REG[rn][1])),REG[rn][1]+dys/3600)
POS={rn:peak(rn) for rn in ('N','S','Bridge')}
for rn,p in POS.items(): P(f'peak pos {rn}: RA {p[0]:.6f} Dec {p[1]:.6f}')
def boxint(v,sp,dvc,box=(-150,150),vfree=450.,vmax=1e9):
    o=np.argsort(v); v=v[o]; sp=sp[o]; fr=(np.abs(v)>vfree)&(np.abs(v)<vmax)&np.isfinite(sp); c=np.polyfit(v[fr],sp[fr],1); y=sp-np.polyval(c,v); sig=float(np.std(y[fr]))
    sel=(v>=box[0])&(v<=box[1]); return float(y[sel].sum()*dvc),float(sig*np.sqrt(sel.sum())*dvc)
def apm(shape,pix_,c,rad):
    Y_,X_=np.mgrid[:shape[0],:shape[1]]; return np.hypot(X_-c[0],Y_-c[1])*pix_<=rad
RES={}
def add(rn,rad,line,val): RES[f'{rn}|{rad}|{line}']=val; P(f' {rn:6s} r={rad}" {line:6s} {val[0]*K:7.2f} +- {val[1]*K:5.2f} K km/s')
# 1-0
for rn,p in POS.items():
    c=pxy(p)
    for rad in (0.6,2.2):
        m=apm((ny,nx),pix,c,rad); sp=Tb10[:,m].mean(axis=1)
        add(rn,rad,'HCN1',boxint(vH,sp,dv,vfree=550.,vmax=1500.)); add(rn,rad,'HCO1',boxint(vO,sp,dv,vfree=550.,vmax=1500.))
# 3-2 comb 2.2" -> target
C32='/workspace/out/q34/%s_2p2_comb.pbcor.fits'
for ln,f0,ky in (('HCN32',265.88618,'HCN3'),('HCOp32',267.55762,'HCO3')):
    h=fits.open(C32%ln); d=np.nan_to_num(h[0].data.astype(float)); hd=h[0].header; w=WCS(hd).celestial
    vv=(hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']); v=C*(1-vv/1e9/f0*1.0)-V0 if False else C*(1-(vv/1e9)/f0)-V0
    bm_=hd['BMAJ']*3600; p_=abs(hd['CDELT2'])*3600
    Tb=d*1.222e6/(f0**2*bm_**2); kk2=kernel(Ct,cov(bm_,bm_,0),p_); Tc=np.empty_like(Tb)
    for i in range(Tb.shape[0]): Tc[i]=fftconvolve(Tb[i],kk2,'same')
    dvc=abs(v[1]-v[0])
    for rn,p in POS.items():
        a,b=w.world_to_pixel_values(*p); c=(float(a),float(b))
        for rad in (0.6,2.2):
            m=apm(Tc.shape[1:],p_,c,rad); add(rn,rad,ky,boxint(v,Tc[:,m].mean(axis=1),dvc,vfree=420.))
# 4-3 native Xb47
for spw,f0,ky in (('spw16',354.50547,'HCN4'),('spw18',356.73424,'HCO4')):
    cut7,nu7,pix7,bm7,pxy7=load(B47%spw,half=22)
    Tb7=cut7*np.float32(tbk(np.mean(nu7),bm7[0],bm7[1])); v7=vel(nu7,f0); dv7=abs(v7[1]-v7[0])
    P(f'Xb47 {spw}: native beam {bm7[0]:.2f}x{bm7[1]:.2f} PA {bm7[2]:.0f}, channel {dv7:.1f} km/s')
    for rn,p in POS.items():
        c=pxy7(p)
        for rad in (0.6,2.2):
            m=apm(Tb7.shape[1:],pix7,c,rad); add(rn,rad,ky,boxint(v7,Tb7[:,m].mean(axis=1),dv7))
json.dump({k:[v[0]*K,v[1]*K] for k,v in RES.items()},open('/workspace/out/q44/s260c.json','w'),indent=1)
open('/workspace/out/q44/s260c.txt','w').write('\n'.join(log))
