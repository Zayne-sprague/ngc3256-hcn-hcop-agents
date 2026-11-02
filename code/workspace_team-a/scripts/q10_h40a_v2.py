# Q10 v2 (reviewer fixes 1-3): H40a free-free L/C expectation, matched/stacked H40a extraction, f_ff, 3-component continuum decomposition.
# python3 -I q10_h40a_v2.py  -> out/q10v2/h40a_v2.log, h40a_v2.json
import sys,json,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve
from q7_lib import Cube,smooth,peak_near,RA0,DE0
rng=np.random.default_rng(42)
log=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); log.append(s)
# ---------- 1. free-free line/continuum (Rohlfs & Wilson 2004 eq. 10.37 / Gordon & Sorochenko 2002 LTE)
def LC_ff(nu_ghz,Te,yp=0.08,integ=True):
    """dv*T_L/T_C [km/s]: 6.985e3 nu^1.1 Te^-1.15 /(1+y+) with dv the FWHM; x1.0645 gives the integral of a Gaussian line (km/s)"""
    v=6.985e3*nu_ghz**1.1*Te**-1.15/(1+yp); return v*(1.0645 if integ else 1.)
nuH=98.11
P(f'FREE-FREE H40a L/C at {nuH} GHz (LTE, optically thin, y+=0.08): FWHM*Tpeak/Tc | integral(Tl dv)/Tc [km/s]')
for Te in (5000,6000,8000,1e4,1.2e4): P(f'  Te={Te:6.0f} K: {LC_ff(nuH,Te,integ=False):5.1f} | {LC_ff(nuH,Te):5.1f}   (y+=0.0: {LC_ff(nuH,Te,0.0):5.1f}) ')
P(' Michiyama+2020 Table 6: R(N)=43+-3 km/s at Te=5900 K, R(S)=24+-2 at Te=10200 K -> formula gives', f'{LC_ff(nuH,5900,integ=False):.0f} and {LC_ff(nuH,10200,integ=False):.0f} (same convention).')
# ---------- 2. data
D='/archive/files/2016.1.00965.S/member.uid___A001_X879_X94a.NGC_3256_sci.spw29.'
hc=fits.getheader(D+'cube.I.pbcor.fits'); cub=np.nan_to_num(fits.getdata(D+'cube.I.pbcor.fits')[0].astype(np.float64))
hm=fits.getheader(D+'mfs.I.pbcor.fits'); mfs=np.nan_to_num(fits.getdata(D+'mfs.I.pbcor.fits')[0,0].astype(np.float64))
nu=hc['CRVAL3']+hc['CDELT3']*np.arange(cub.shape[0]); REST=99.02295e9
vel=299792.458*(REST/nu-1); VSYS=2780.; v=vel-VSYS; dv=abs(np.median(np.diff(v)))
if np.median(np.diff(v))<0: v=v[::-1].copy(); cub=cub[::-1].copy()   # np.interp needs ascending velocity
bm=(hc['BMAJ']*3600,hc['BMIN']*3600,hc['BPA']); pix=abs(hc['CDELT1'])*3600
w=WCS(hc).celestial; ny,nx=cub.shape[1:]
yy,xx=np.mgrid[0:ny,0:nx]; RAp,DEp=w.pixel_to_world_values(xx,yy)
X=(RAp-RA0)*np.cos(np.deg2rad(DE0))*3600; Y=(DEp-DE0)*3600
P(f'\nH40a cube: beam {bm[0]:.2f}x{bm[1]:.2f}, pix {pix:.3f}", dv={dv:.1f} km/s, v-VSYS {v.min():.0f}..{v.max():.0f}')
# baseline: linear from |v|>700 per pixel
base=np.abs(v)>700
def lin_sub(c,vv,bmask):
    A=np.vstack([np.ones(bmask.sum()),vv[bmask]/500.]).T; cc=np.tensordot(np.linalg.pinv(A),c[bmask],axes=(1,0))
    return c-np.tensordot(np.vstack([np.ones(len(vv)),vv/500.]).T,cc,axes=(1,0))
sub=lin_sub(cub,v,base)
# ---------- HCN velocity field (Xb4d 0.95") for kinematic registration
H=Cube('Xb4d','HCN'); O=Cube('Xb4d','HCO+')
tb=(1.05,1.05,0.)
cen=(np.abs(H.v)<350)
m0=smooth(H.m0(-350,350),H,tb); m1=smooth((H.cube[cen]*H.v[cen,None,None]).sum(0)*H.dv,H,tb)
vmap=np.where(m0>0.06,m1/np.where(m0>0.06,m0,1),np.nan)   # centroid relative to VSYS
# map to H40a grid
def to_H40(img):
    xs=H.X[0]; ys=H.Y[:,0]
    ix=np.interp(X,xs[::-1],np.arange(len(xs))[::-1]) if xs[0]>xs[-1] else np.interp(X,xs,np.arange(len(xs)))
    iy=np.interp(Y,ys[::-1],np.arange(len(ys))[::-1]) if ys[0]>ys[-1] else np.interp(Y,ys,np.arange(len(ys)))
    return map_coordinates(np.nan_to_num(img,nan=np.nan),[iy,ix],order=1,mode='nearest')
vfield=to_H40(np.where(np.isnan(vmap),-9999,vmap)); vfield=np.where(vfield<-5000,np.nan,vfield)
REG={'N':(0.0,2.63),'S':(-0.15,-2.64)}
PK={k:peak_near(O,smooth(O.m0(-350,350),O,tb),*c) for k,c in REG.items()}
for k in PK:
    r=np.hypot(X-PK[k][0],Y-PK[k][1]); P(f'{k}: HCO+ peak {PK[k][0]:+.2f},{PK[k][1]:+.2f}; HCN centroid velocity (rel. {VSYS:.0f}) core {np.nanmedian(vfield[r<0.5]):+.0f}, r<1.5 range {np.nanpercentile(vfield[r<1.5],5):+.0f}..{np.nanpercentile(vfield[r<1.5],95):+.0f} km/s')
# ---------- smooth cube+mfs to 1.05" circular and also keep native
def smooth_to(img,bmaj,bmin,bpa,target):
    sa=np.sqrt(max(target**2-bmaj**2,0))/2.3548/pix; sb=np.sqrt(max(target**2-bmin**2,0))/2.3548/pix
    if sa==0 and sb==0: return img
    n=int(5*max(sa,sb))+1; y_,x_=np.mgrid[-n:n+1,-n:n+1].astype(float); th=np.deg2rad(bpa)
    u=-x_*np.sin(th)+y_*np.cos(th); t=x_*np.cos(th)+y_*np.sin(th)
    k=np.exp(-0.5*((u/max(sa,1e-3))**2+(t/max(sb,1e-3))**2)); k/=k.sum()
    return fftconvolve(img,k,mode='same')*(target**2)/(bmaj*bmin)
def smooth_cube(c,target):
    return np.array([smooth_to(ch,*bm,target) for ch in c])
def kernel_cube(c,target):  # target Jy/beam(target)
    return smooth_cube(c,target)
CUBES={'native':(sub,mfs,bm[:2])}
sub105=smooth_cube(sub,1.05); mfs105=smooth_to(mfs,hm['BMAJ']*3600,hm['BMIN']*3600,hm['BPA'],1.05)
CUBES['1.05']=(sub105,mfs105,(1.05,1.05))
sub15=smooth_cube(sub,1.5); mfs15=smooth_to(mfs,hm['BMAJ']*3600,hm['BMIN']*3600,hm['BPA'],1.5)
CUBES['1.5']=(sub15,mfs15,(1.5,1.5))
# ---------- registered spectrum and matched-filter amplitude
def reg_spec(cube,mask,register=True,vref=0.):
    idx=np.argwhere(mask); out=np.zeros(len(v)); n=0
    for j,i in idx:
        s=cube[:,j,i]
        sh=vfield[j,i] if register and np.isfinite(vfield[j,i]) else 0.
        out+=np.interp(v+ (sh-vref),v,s,left=np.nan,right=np.nan) if False else np.interp(v,v-0,s) if False else np.interp(v+sh,v,s)
        n+=1
    return out/max(n,1)
# NB: shifting channel v+sh samples the spectrum at v+sh, i.e. moves a line centred at sh to 0.
def gauss(v,fw,c=0.): return np.exp(-0.5*((v-c)/(fw/2.3548))**2)
def mf_amp(spec,fw,lo=-1.5,hi=1.5):
    # amplitude of Gaussian template (fixed centre 0, fixed FWHM) by least squares in |v|<2.5 FWHM; integral = A*fw*1.0645
    k=np.abs(v)<2.5*fw; g=gauss(v,fw); ok=k&np.isfinite(spec)
    return (spec[ok]*g[ok]).sum()/(g[ok]**2).sum(), g[ok]
def amp_noise(cube,mask_fn,fw,n=60,register=True):
    """noise of the matched-filter integral from off-source apertures of identical geometry (centres on a ring/box in the line-free field)"""
    vals=[]
    tries=0
    while len(vals)<n and tries<2000:
        tries+=1
        cx,cy=rng.uniform(-9,9),rng.uniform(-9,9)
        if min(np.hypot(cx-PK['N'][0],cy-PK['N'][1]),np.hypot(cx-PK['S'][0],cy-PK['S'][1]))<4.5 or np.hypot(cx,cy)<3.5: continue
        m=mask_fn(cx,cy)
        if m.sum()<3: continue
        s=reg_spec(cube,m,register=False)   # no velocity field outside the nuclei
        a,_=mf_amp(s,fw); vals.append(a*fw*1.0645)
    return np.std(vals),np.mean(vals),len(vals)
def vel_null(spec,fw):
    # matched-filter amplitude (integral) with the template shifted to line-free velocities (same aperture & registration)
    vals=[]
    for c in list(np.arange(-3000,-1000,2.0*fw*0.8))+list(np.arange(1000,2400,2.0*fw*0.8)):
        k=np.abs(v-c)<2.5*fw; g=gauss(v,fw,c); ok=k&np.isfinite(spec)
        if ok.sum()<0.9*k.sum() or k.sum()<5: continue
        vals.append((spec[ok]*g[ok]).sum()/(g[ok]**2).sum()*fw*1.0645)
    return np.std(vals),len(vals)
def ap_mask(cx,cy,rin,rout):
    r=np.hypot(X-cx,Y-cy); return (r>=rin)&(r<rout)
# line widths from HCN registered at 1.05" in the same aperture (fit FWHM by Gaussian)
def hcn_fwhm(cx,cy,rin,rout):
    r=np.hypot(H.X-cx,H.Y-cy); m=(r>=rin)&(r<rout)
    vm=vmap[m]; idx=np.argwhere(m)
    sp=np.zeros(len(H.v)); nn=0
    for j,i in idx:
        sh=vmap[j,i] if np.isfinite(vmap[j,i]) else 0.; sp+=np.interp(H.v+sh,H.v,H.cube[:,j,i]); nn+=1
    sp/=nn
    from scipy.optimize import curve_fit
    k=np.abs(H.v)<500
    f=lambda x,a,c,s,b: a*np.exp(-0.5*((x-c)/s)**2)+b
    p,_=curve_fit(f,H.v[k],sp[k],p0=[sp.max(),0,70,0]); return abs(p[2])*2.3548,sp
RES={}
P('\n=== H40a matched/stacked extraction: velocity-registered on the HCN centroid field, Gaussian template with FWHM fixed to the registered HCN line in the same aperture')
P(' amplitude integral I in mJy/beam km/s (aperture-mean pbcor), noise from 60 off-source apertures of identical shape (incl. baseline ripple); C = mean 98 GHz continuum in the aperture (mJy/b); L/C in km/s')
for res_name in ('native','1.05','1.5'):
    cube,cont,bb=CUBES[res_name]
    for reg in ('N','S'):
        for (rin,rout) in ((0,0.5),(0,1.0),(0,1.5),(0,2.0),(0.5,1.0),(1.0,1.5),(1.5,2.5)):
            cx,cy=PK[reg]
            m=ap_mask(cx,cy,rin,rout)
            # restrict to pixels nearer to this nucleus
            ox,oy=PK['S' if reg=='N' else 'N']; m&=(np.hypot(X-cx,Y-cy)<np.hypot(X-ox,Y-oy))
            fw,hsp=hcn_fwhm(cx,cy,rin,rout)
            sp_r=reg_spec(cube,m,True); sp_u=reg_spec(cube,m,False)
            a_r,_=mf_amp(sp_r,fw); a_u,_=mf_amp(sp_u,fw)
            mk=lambda a,b_: ap_mask(a,b_,rin,rout)
            en,mn,nn=amp_noise(cube,mk,fw)
            evn,nvn=vel_null(sp_r,fw); e_off=en*1e3; e_vel=evn*1e3
            I=a_r*fw*1.0645*1e3; e=max(e_off,e_vel); Iu=a_u*fw*1.0645*1e3
            C=cont[m].mean()*1e3; eC=0.063*max(1,1.0)/np.sqrt(max(1,m.sum()*pix**2/(1.1331*bb[0]*bb[1])))
            LC=I/C; eLC=e/C
            RES[f'{res_name}|{reg}|{rin}-{rout}']=dict(I=I,eI=e,Iu=Iu*1.0,C=C,LC=LC,eLC=eLC,fwhm=fw,npix=int(m.sum()),bias=mn*1e3)
            P(f' {res_name:6s} {reg} r={rin}-{rout}": FWHM_HCN={fw:5.0f} km/s npix={m.sum():4d} I={I:7.1f}+-{e:5.1f} (S/N {I/e:4.1f}; off-src {e_off:5.1f}, vel-null {e_vel:5.1f}) [unreg. {Iu:6.1f}] C={C:6.3f} L/C={LC:6.1f}+-{eLC:5.1f} km/s (off-source mean {mn*1e3:+.1f})')
# ---------- f_ff from L/C
P('\n=== f_ff(98 GHz) = (L/C)meas / (L/C)ff ; (L/C)ff integral convention; Te range 5000-12000 K (Te=1e4 central), y+=0.08, b_n*beta_n=1 assumed (<~10% effect at n=40)')
FF={}
for key,rr in RES.items():
    if not key.startswith(('1.05','native')): continue
    out=[]
    for Te in (5000,1e4,1.2e4):
        L=LC_ff(nuH,Te); out.append((rr['LC']/L,rr['eLC']/L))
    FF[key]=out
    P(f' {key:22s} f_ff(Te=1e4)={out[1][0]:5.2f}+-{out[1][1]:4.2f}   Te=5000: {out[0][0]:5.2f}+-{out[0][1]:4.2f}   Te=1.2e4: {out[2][0]:5.2f}+-{out[2][1]:4.2f}   3sigma UL(Te=1e4)={(max(out[1][0],0)+3*out[1][1]):.2f}')
# ---------- Michiyama+2020 comparison: flux inside large apertures (Jy km/s), 1.5" taper-like beam
P('\n=== comparison with Michiyama+2020 (Table 6: N 273+-14, S 143+-7 mJy km/s, beam 1.48x1.31", aperture 2.5x source size ~ >2")')
cube,cont,bb=CUBES['1.5']
for reg in ('N','S'):
    cx,cy=PK[reg]; ox,oy=PK['S' if reg=='N' else 'N']
    for rad in (1.5,2.0,2.5):
        m=ap_mask(cx,cy,0,rad)&(np.hypot(X-cx,Y-cy)<np.hypot(X-ox,Y-oy))
        fw,_=hcn_fwhm(cx,cy,0,rad); sp=reg_spec(cube,m,True); a,_=mf_amp(sp,fw)
        omega=1.1331*1.5**2/pix**2   # pixels per beam
        # total flux = sum over pixels / pix-per-beam ; matched filter amplitude applies to the aperture mean
        F=a*fw*1.0645*1e3*m.sum()/omega
        mk=lambda a_,b_,rad=rad: ap_mask(a_,b_,0,rad); en,_,_=amp_noise(cube,mk,fw,n=40); eF=en*fw/ (fw*1.0645)*0  # placeholder
        eF=en*1e3*m.sum()/omega
        Cflux=cont[m].sum()*1e3/omega
        P(f' {reg} r<{rad}" (1.5" beam): H40a flux {F:6.0f}+-{eF:4.0f} mJy km/s, 98 GHz continuum {Cflux:5.1f} mJy; L/C={F/Cflux:5.1f}+-{eF/Cflux:4.1f} km/s')
        RES[f'mich|{reg}|{rad}']=dict(F=F,eF=eF,Cflux=Cflux)
json.dump(dict(RES=RES,FF={k:v for k,v in FF.items()}),open('/workspace/out/q10v2/h40a_v2.json','w'),indent=1)
open('/workspace/out/q10v2/h40a_v2.log','w').write('\n'.join(log))
