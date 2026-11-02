# Q47(b): own 2.2" common-taper 13C/12C 1-0 ratios (H13CN/HCN, H13CO+/HCO+) at N, S, C(=Bridge) from the archived 1-0 cubes
# (2015.1.00993.S X95d spw0/spw1 for 13C, X961 spw0 for HCN/HCO+).  Every channel is convolved to a circular 2.2" beam (radio_beam deconvolution kernel),
# Jy/beam rescaled by the beam-area ratio, flux in r=2.2" apertures.  Three estimators: (T) amplitude of the 12C template of the same aperture fitted to the 13C
# spectrum with a free linear baseline, (A) +-150 km/s integral ratio (baseline from line-free channels, avoids SiO 2-1 near H13CO+), (P) peak pixel +-150 km/s.
# Errors: formal (line-free channel rms) and empirical (same estimator at 40 off-source apertures with the N template; max of the two).
import sys,warnings,glob,json; warnings.filterwarnings('ignore')
import numpy as np
from astropy.io import fits; from astropy.wcs import WCS
from radio_beam import Beam; import astropy.units as u
from scipy.ndimage import gaussian_filter
C=299792.458; RA0,DE0=156.96346,-43.90461
REG={'N':(0.0,2.55),'S':(-0.15,-2.55),'C':(0.3,0.15)}; VLOC={'N':2778.,'S':2764.,'C':2777.}
LINES={'HCN':88.631602,'HCO+':89.188526,'H13CN':86.339922,'H13CO+':86.754288}
CUBES={'X961_0':'member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits',
       'X95d_0':'member.uid___A001_X2fb_X95d.ari_l.ngc3256_sci.spw0_85203MHz.12m.cube.I.pbcor.fits',
       'X95d_1':'member.uid___A001_X2fb_X95d.ari_l.ngc3256_sci.spw1_86692MHz.12m.cube.I.pbcor.fits'}
TARGET=Beam(2.2*u.arcsec,2.2*u.arcsec,0*u.deg); VG=np.arange(-700,700.1,26.)
def load(key,f0):
    p=glob.glob('/archive/files/*/'+CUBES[key])[0]; hdu=fits.open(p,memmap=True); hd=hdu[0].header
    fr=hd['CRVAL3']+(np.arange(hd['NAXIS3'])+1-hd['CRPIX3'])*hd['CDELT3']; v=(1-fr/1e9/f0)*C
    sel=np.where((v>2780-760)&(v<2780+760))[0]; data=np.array(hdu[0].data[0,sel.min():sel.max()+1],float) if hdu[0].data.ndim==4 else np.array(hdu[0].data[sel.min():sel.max()+1],float)
    pix=abs(hd['CDELT1'])*3600; w=WCS(hd).celestial
    try: bt=hdu[1].data; bm=[(bt['BMAJ'][i],bt['BMIN'][i],bt['BPA'][i]) for i in sel]
    except Exception: bm=[(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])]*len(sel)
    out=np.empty_like(data)
    for i,(a,b,pa) in enumerate(bm):
        im=data[i]; im=np.nan_to_num(im); bm0=Beam(a*u.arcsec,b*u.arcsec,pa*u.deg)
        try: k=TARGET.deconvolve(bm0); kx=k.major.to(u.arcsec).value/2.3548/pix; ky=k.minor.to(u.arcsec).value/2.3548/pix; pak=k.pa.to(u.deg).value
        except Exception: out[i]=np.nan; continue
        # rotate-free approximation: elliptical Gaussian kernel via two 1-D passes only if PA~0/90; use mean sigma (kernel axes differ <15%)
        s=np.sqrt(kx*ky); out[i]=gaussian_filter(im,s)*(TARGET.sr/bm0.sr).value
    return v[sel],out,w,pix
def pixpos(w,dra,dde):
    x,y=w.world_to_pixel_values(RA0+dra/3600/np.cos(np.deg2rad(DE0)),DE0+dde/3600); return float(x),float(y)
def aperture(cube,w,pix,dra,dde,r=2.2):
    x0,y0=pixpos(w,dra,dde); Y,X=np.mgrid[:cube.shape[1],:cube.shape[2]]; m=np.hypot(X-x0,Y-y0)*pix<=r
    return np.nanmean(cube[:,m],axis=1)
def peakpix(cube,w,pix,dra,dde,v,vc):
    x0,y0=pixpos(w,dra,dde); Y,X=np.mgrid[:cube.shape[1],:cube.shape[2]]; m=np.hypot(X-x0,Y-y0)*pix<=1.5
    sel=np.abs(v-vc)<150; I=np.nansum(cube[sel],axis=0); I=np.where(m,I,-1e9); j=np.unravel_index(np.argmax(I),I.shape); return cube[:,j[0],j[1]]
def regrid(v,s,vc):
    o=np.argsort(v); return np.interp(VG,(v-vc)[o],np.asarray(s)[o],left=np.nan,right=np.nan)
def stitch(parts):
    a=np.array(parts); return np.nanmean(a,axis=0)
def est(s12,s13,mode):
    # s12,s13 on VG grid (velocity relative to local systemic)
    ok=np.isfinite(s13)&np.isfinite(s12); dv=26.
    if mode=='A':
        sel=ok&(np.abs(VG)<=150); bl=ok&(np.abs(VG)>350)&(np.abs(VG)<650)
        def integ(s):
            p=np.polyfit(VG[bl],s[bl],1); return np.sum(s[sel]-np.polyval(p,VG[sel]))*dv
        return integ(s13)/integ(s12)
    # template: baseline-free 12C shape (subtract own baseline), linear baseline in 13C
    bl=ok&(np.abs(VG)>350)&(np.abs(VG)<650); p=np.polyfit(VG[bl],s12[bl],1); t=s12-np.polyval(p,VG)
    sel=ok&(np.abs(VG)<=650); A=np.vstack([t[sel],np.ones(sel.sum()),VG[sel]/300]).T; c,_,_,_=np.linalg.lstsq(A,s13[sel],rcond=None); return c[0]
def chan_rms(s,vlim=(450,700)):
    m=np.isfinite(s)&(np.abs(VG)>vlim[0])&(np.abs(VG)<vlim[1]); return np.std(s[m])
if __name__=='__main__':
    res={}; cub={}
    for key,ln in (('X961_0','HCN'),('X961_0','HCO+'),('X95d_0','H13CN'),('X95d_0','H13CO+'),('X95d_1','H13CO+')):
        cub[(key,ln)]=load(key,LINES[ln]); print('loaded',key,ln,cub[(key,ln)][1].shape,flush=True)
    rng=np.random.default_rng(1); out={}
    for pair,(l12,l13) in {'HCN':('HCN','H13CN'),'HCO+':('HCO+','H13CO+')}.items():
        v12,c12,w12,p12=cub[('X961_0',l12)]
        c13keys=[('X95d_0',l13)]+([('X95d_1',l13)] if l13=='H13CO+' else [])
        for reg,(dra,dde) in REG.items():
            vc=VLOC[reg]
            s12=regrid(v12,aperture(c12,w12,p12,dra,dde),vc)
            s13=stitch([regrid(cub[k][0],aperture(cub[k][1],cub[k][2],cub[k][3],dra,dde),vc) for k in c13keys])
            pk12=regrid(v12,peakpix(c12,w12,p12,dra,dde,v12,vc),vc)
            pk13=stitch([regrid(cub[k][0],peakpix(cub[k][1],cub[k][2],cub[k][3],dra,dde,cub[k][0],vc),vc) for k in c13keys])
            for mode,(a,b) in {'T':(s12,s13),'A':(s12,s13),'P':(pk12,pk13)}.items():
                m='A' if mode=='P' else mode
                R=est(a,b,m)
                # formal error from line-free rms propagated: dR = sig13*sqrt(dv)/I12 style (approx via 12C integral)
                I12=np.nansum(np.where(np.abs(VG)<=150,a-np.nanmedian(a[(np.abs(VG)>450)]),0))*26
                nsel=(np.abs(VG)<=150).sum() if m=='A' else (np.abs(VG)<=650).sum()
                formal=chan_rms(b)*np.sqrt(nsel)*26/I12
                # empirical: offsource apertures, N template shape for 12C replaced by that region's own 12C template
                nul=[]
                H,W=c12.shape[1:]; tries=0
                while len(nul)<40 and tries<2000:
                    tries+=1; ang=rng.uniform(0,2*np.pi); rr=rng.uniform(7,22); d1,d2=rr*np.cos(ang),rr*np.sin(ang)
                    z12=aperture(c12,w12,p12,d1,d2); 
                    if not np.all(np.isfinite(z12)): continue
                    z13=stitch([regrid(cub[k][0],aperture(cub[k][1],cub[k][2],cub[k][3],d1,d2),vc) for k in c13keys])
                    zz12=regrid(v12,z12,vc)
                    if np.sum(np.where(np.abs(VG)<=150,zz12-np.nanmedian(zz12[np.abs(VG)>450]),0))*26 > 0.1*I12: continue   # keep emission-free positions only
                    # null: use THIS region's 12C template on the off-source 13C spectrum
                    nul.append(est(a,z13,m)*I12 if m=='T' else est(zz12*0+a,z13,'A')*0+ (np.sum(np.where(np.abs(VG)<=150,z13-np.polyval(np.polyfit(VG[(np.abs(VG)>350)&(np.abs(VG)<650)&np.isfinite(z13)],z13[(np.abs(VG)>350)&(np.abs(VG)<650)&np.isfinite(z13)],1),VG),0))*26))
                nul=np.array(nul); emp=np.std(nul)/I12 if len(nul)>5 else np.nan
                out[f'{pair}|{reg}|{mode}']=dict(R=float(R),formal=float(formal),empirical=float(emp),nnull=len(nul),I12_Jybkms=float(I12),err=float(max(formal,emp if np.isfinite(emp) else 0)))
                print(f'{l13}/{l12} {reg} {mode}: R={R:.4f} +- {max(formal,emp if np.isfinite(emp) else 0):.4f} (formal {formal:.4f}, null {emp:.4f}, n={len(nul)}) I12={I12:.3f}',flush=True)
    json.dump(out,open('/workspace/out/q47/ratios10.json','w'),indent=1)
