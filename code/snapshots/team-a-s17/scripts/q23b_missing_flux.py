# Q23 round 2: missing-flux budget for CO(2-1) T_b at the S/N cores using the 12m 0.18" (2015.1.00714 spw1), 12m 0.82x0.68" (2015.1.00902),
# 7m mosaic (2015.1.00714 spw1) and SD (spw17) cubes. Image-plane "feather" correction: dT(v)=T_ref(core)-T_hi(convolved to ref beam)(core), added to hi-res T_b.
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve_fft
from astropy import units as u
from radio_beam import Beam
C=299792.458
RA0,DE0=156.96346,-43.90461
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
A='/archive/files/2015.1.00714.S/'; B='/archive/files/2015.1.00902.S/'
F={'hi':A+'member.uid___A001_X2d6_Xae.ari_l.NGC_3256_sci.spw1_229288MHz.12m.cube.I.pbcor.fits',
   'mid':B+'member.uid___A001_X2fb_X881.ari_l.ngc_3256_sci.spw2_228401MHz.12m.cube.I.pbcor.fits',
   '7m':A+'member.uid___A001_X2d6_Xb2.NGC_3256_sci.spw1_229289MHz.mos.7m.cube.I.pbcor.fits',
   'sd':A+'member.uid___A001_X2d6_Xb6.NGC_3256_sci.spw17.cube.I.sd.fits'}
VG=np.arange(2300,2845,5.08)  # common velocity grid for comparison (km/s), 5.08 = 0.82" cube channel
def TRJ(S,nu0,a,b): return 1.222e3*S*1e3/((nu0/1e9)**2*a*b)
def load(key,half):
    h=fits.open(F[key],memmap=True); hd=h[0].header; cube=h[0].data[0]
    nu=hd['CRVAL3']+(np.arange(hd['NAXIS3']))*hd['CDELT3']; v=(1-nu/230.538e9)*C
    k=np.where((v>VG.min()-15)&(v<VG.max()+15))[0]; k=k[np.argsort(v[k])]
    w=WCS(hd).celestial; pix=abs(hd['CDELT2'])*3600
    bm=Beam(hd['BMAJ']*u.deg,hd['BMIN']*u.deg,hd['BPA']*u.deg)
    return dict(hd=hd,cube=cube,v=v[k],nu=nu[k],k=k,w=w,pix=pix,beam=bm)
def cut(d,x0,y0,half):
    ra=RA0+x0/3600/np.cos(np.deg2rad(DE0)); de=DE0+y0/3600
    fx,fy=[float(a) for a in d['w'].world_to_pixel_values(ra,de)]; ix,iy=int(round(fx)),int(round(fy)); HW=int(round(half/d['pix']))
    ny,nx=d['cube'].shape[1:]
    y1,y2,x1,x2=max(iy-HW,0),min(iy+HW+1,ny),max(ix-HW,0),min(ix+HW+1,nx)
    sub=np.nan_to_num(np.array(d['cube'][d['k']][:,y1:y2,x1:x2],float))
    return sub,fx-x1,fy-y1
def interp_at(S,cx,cy):
    x0,y0=int(np.floor(cx)),int(np.floor(cy)); fx,fy=cx-x0,cy-y0
    return S[:,y0,x0]*(1-fx)*(1-fy)+S[:,y0,x0+1]*fx*(1-fy)+S[:,y0+1,x0]*(1-fx)*fy+S[:,y0+1,x0+1]*fx*fy
def conv(sub,src_beam,pix,tgt_beam):
    ker=tgt_beam.deconvolve(src_beam).as_kernel(pix*u.arcsec)
    fac=(tgt_beam.sr/src_beam.sr).decompose().value
    return np.array([convolve_fft(s,ker,boundary='fill',fill_value=0,allow_huge=True,normalize_kernel=True) for s in sub])*fac
def spec_T(S,d,bm):  # RJ K on VG
    T=TRJ(S,d['nu'],bm.major.to(u.arcsec).value,bm.minor.to(u.arcsec).value)
    return np.interp(VG,d['v'],T)
def boxavg(d,T,width=5.08):  # average a 2.54 km/s spectrum to width
    return T
res={}
D={k:load(k,0) for k in F}
for nm,(x0,y0) in CORES.items():
    r={}
    # ---- A: 0.18" convolved to the 0.82x0.68 PA-0 beam (proper elliptical) vs 0.82" cube, at the core
    sub,cx,cy=cut(D['hi'],x0,y0,4.0)
    Sc=conv(sub,D['hi']['beam'],D['hi']['pix'],D['mid']['beam'])
    # average the 2.54 km/s channels onto 5.08 km/s boxes (0.82 cube channel width) before comparing
    Tc=TRJ(interp_at(Sc,cx,cy),D['hi']['nu'],D['mid']['beam'].major.to(u.arcsec).value,D['mid']['beam'].minor.to(u.arcsec).value)
    vh=D['hi']['v']; Tbox=np.array([Tc[np.abs(vh-vv)<=2.54].mean() for vv in VG])
    subm,mx,my=cut(D['mid'],x0,y0,4.0)
    Tm=np.interp(VG,D['mid']['v'],TRJ(interp_at(subm,mx,my),D['mid']['nu'],D['mid']['beam'].major.to(u.arcsec).value,D['mid']['beam'].minor.to(u.arcsec).value))
    # core sits between pixels; also take best-aligned offset: compare peak and mean over line (v within +-100 of peak)
    kp=np.argmax(Tm); sel=np.abs(VG-VG[kp])<100
    r['A']={'T_mid_peak':float(Tm[kp]),'T_hi_conv_at_same_v':float(Tbox[kp]),'ratio_peak':float(Tbox[kp]/Tm[kp]),'ratio_line_mean_pm100':float(Tbox[sel].mean()/Tm[sel].mean()),'v_peak_mid':float(VG[kp]),
        'T_mid_peak_any':float(Tm.max()),'T_hi_conv_peak_any':float(Tbox.max()),'ratio_of_peaks':float(Tbox.max()/Tm.max())}
    # ---- B: 12m 0.82 cube convolved to 7m beam vs 7m cube
    sub,cx,cy=cut(D['mid'],x0,y0,25.0)
    b7=D['7m']['beam']
    S7=conv(sub,D['mid']['beam'],D['mid']['pix'],b7)
    T12_7=np.interp(VG,D['mid']['v'],TRJ(interp_at(S7,cx,cy),D['mid']['nu'],b7.major.to(u.arcsec).value,b7.minor.to(u.arcsec).value))
    sub7,qx,qy=cut(D['7m'],x0,y0,25.0)
    T7=np.interp(VG,D['7m']['v'],TRJ(interp_at(sub7,qx,qy),D['7m']['nu'],b7.major.to(u.arcsec).value,b7.minor.to(u.arcsec).value))
    sel=np.abs(VG-VG[np.argmax(T7)])<150
    r['B']={'T7m_peak':float(T7.max()),'T12m082_conv7m_peak':float(T12_7.max()),'ratio_peak':float(T12_7[np.argmax(T7)]/T7.max()),'ratio_line_sum_pm150':float(T12_7[sel].sum()/T7[sel].sum()),'v_peak':float(VG[np.argmax(T7)])}
    dT1=T7-T12_7
    # ---- C: 7m convolved to SD beam vs SD
    bsd=D['sd']['beam']
    sub7b,qx,qy=cut(D['7m'],x0,y0,45.0)
    Ssd=conv(sub7b,b7,D['7m']['pix'],bsd)
    T7_sd=np.interp(VG,D['7m']['v'],TRJ(interp_at(Ssd,qx,qy),D['7m']['nu'],bsd.major.to(u.arcsec).value,bsd.minor.to(u.arcsec).value))
    subs,sx,sy=cut(D['sd'],x0,y0,30.0)
    Tsd=np.interp(VG,D['sd']['v'],TRJ(interp_at(subs,sx,sy),D['sd']['nu'],bsd.major.to(u.arcsec).value,bsd.minor.to(u.arcsec).value))
    sel=np.abs(VG-VG[np.argmax(Tsd)])<250
    r['C']={'Tsd_peak':float(Tsd.max()),'T7m_conv_sd_peak':float(T7_sd.max()),'ratio_line_sum_pm250':float(T7_sd[sel].sum()/Tsd[sel].sum()),'v_peak':float(VG[np.argmax(Tsd)])}
    dT2=Tsd-T7_sd
    r['dT1_7m_minus_12m']=dT1.tolist(); r['dT2_sd_minus_7m']=dT2.tolist(); r['VG']=VG.tolist()
    r['T7']=T7.tolist(); r['T12_7']=T12_7.tolist(); r['Tsd']=Tsd.tolist(); r['T7_sd']=T7_sd.tolist()
    res[nm]=r
    print(nm,'A',{k:round(x,2) for k,x in r['A'].items()}); print(nm,'B',{k:round(x,2) for k,x in r['B'].items()}); print(nm,'C',{k:round(x,2) for k,x in r['C'].items()})
    print(nm,'dT1 (7m-12m@7m beam) at line peak range:',np.round(dT1[np.abs(VG-VG[np.argmax(T7)])<60],2),' dT2:',np.round(dT2[np.abs(VG-VG[np.argmax(T7)])<60],2))
json.dump(res,open('/workspace/out/q23b/missing_flux.json','w'),indent=1)
