# Re-identification of the Xb4f ~254.0 GHz (rest) feature: velocity-centroid test vs HCO+ 3-2 in the same pixel.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import curve_fit
from regions import REGIONS
C=299792.458
def load(ous,spw):
    h=fits.open(f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.spw{spw}.cube.I.pbcor.fits',memmap=True)
    H=h[0].header; d=np.squeeze(h[0].data)
    nu=(H['CRVAL3']+(np.arange(H['NAXIS3'])-H['CRPIX3']+1)*H['CDELT3'])/1e9
    return d,H,nu
def spec(d,H,ra,dec):
    x,y=WCS(H).celestial.all_world2pix(ra,dec,0); return d[:,int(round(float(y))),int(round(float(x)))]*1e3
def gl(x,a,x0,w,b0,b1): return a*np.exp(-0.5*((x-x0)/(w/2.3548))**2)+b0+b1*x
ous=sys.argv[1] if len(sys.argv)>1 else 'Xb4f'
SP={'Xb4f':('20','16'),'Xb4d':('29','25'),'Xb4b':('29','25')}[ous]
dA,HA,nuA=load(ous,SP[0]); dB,HB,nuB=load(ous,SP[1])
HCO=267.557626
CAND={'CH3OH 2(0)-1(-1)E 254.01538':254.015377,'NS 253.968 (hf-weighted f)':253.9692,'SO2 253.957':253.956567,
      'c-HCCCH/c-C3H 254.014':254.0143,'SiS 14-13 254.103':254.103}
out=[]
print('ous',ous)
for reg in ('S','Bridge','N','Sblue'):
    ra,dec=REGIONS[reg][:2]
    sB=spec(dB,HB,ra,dec); vB=C*(1-nuB/HCO)
    m=(vB>2200)&(vB<3400); o=((vB>1500)&(vB<2100))|((vB>3500)&(vB<4100))
    # reference: HCO+ 3-2 Gaussian + linear baseline; baseline fixed from line-free channels
    p,cv=curve_fit(gl,vB[m|o],sB[m|o],p0=[sB[m].max(),2750,150,0,0])
    vh,evh,fw=p[1],np.sqrt(cv[1,1]),p[2]
    # 254 window, in velocity relative to the CH3OH rest freq; fixed baseline from line-free obs freq
    rest=254.015377; vA=C*(1-nuA/rest)
    sA=spec(dA,HA,ra,dec)
    sel=(vA>1200)&(vA<4300)
    line=(vA>2350)&(vA<3250)
    pf=np.polyfit(vA[sel&~line],sA[sel&~line],1); res=sA[sel&~line]-np.polyval(pf,vA[sel&~line])
    rms=res.std(); r1=np.corrcoef(res[:-1],res[1:])[0,1]
    s2=sA-np.polyval(pf,vA)
    # centroid in window of fixed width +-450 km/s about the HCO+ velocity (CH3OH frame), error by MC from scatter (lag-1 corr via AR(1))
    def cen(vc,hw,s):
        w=np.abs(vA-vc)<hw; ss=s[w]; return (ss*vA[w]).sum()/ss.sum(), ss.sum()*3.91e-3*0+ss.sum()
    # template-shift fit: observed HCO+ 3-2 profile (same pixel, normalised) shifted by dv and scaled by A; baseline fixed
    from scipy.interpolate import interp1d
    tmpl=interp1d(vB[::-1]-vh,(sB-np.polyval(np.polyfit(vB[o],sB[o],1),vB))[::-1]/ (sB[m].max()),bounds_error=False,fill_value=0.)
    fw_=(vA>1800)&(vA<3700)
    dvs=np.arange(-250,250.1,1.0); chi=[];amps=[]
    for dv in dvs:
        t=tmpl(vA[fw_]-vh-dv); A=(t*s2[fw_]).sum()/(t*t).sum(); amps.append(A); chi.append((((s2[fw_]-A*t)/rms)**2).sum())
    chi=np.array(chi); k=chi.argmin(); infl=(1+r1)/(1-r1) if r1>0 else 1.0
    ok=(chi-chi[k])*1.0/infl<=1.0; lo,hi=dvs[ok].min(),dvs[ok].max()
    dv0=dvs[k]; edv=max((hi-lo)/2,1.0)
    chi0=(((s2[fw_])/rms)**2).sum()
    print(f'{reg}: HCO+3-2 (same pixel) v0={vh:.0f} FWHM={fw:.0f}; 254 window rms={rms:.2f} mJy/ch lag1={r1:.2f}; template-fit peak {amps[k]:.2f} mJy/beam (HCO+ 3-2 peak {sB[m].max():.0f}); dchi2 vs null={(chi0-chi[k])/infl:.1f}; shift of feature in CH3OH frame vs HCO+ profile: {dv0:+.0f}+-{edv:.0f} km/s')
    q=[amps[k],vh+dv0,fw]; eq=[0,0,0]; eq[1]=edv; vh_use=vh; evh=1.0
    for kk,f in CAND.items():
        v=C*(1-rest*(1-q[1]/C)/f); d=v-vh; sig=np.hypot(eq[1],evh)
        print(f'    {kk:32s} implied v={v:7.0f}  minus HCO+ v0 = {d:+6.0f} +- {sig:.0f} km/s  ({d/sig:+.1f} sigma)')
        out.append(dict(region=reg,cand=kk,v=v,dv=d,err=sig,vhco=vh,A_mJy=q[0],dchi2=(chi0-chi[k])/infl,rms=rms,lag1=r1))
    # integrated intensity, window +-HCO+ fwhm-ish, CH3OH assumption
    w=np.abs(vA-q[1])<fw/2*1.5; 
    I=s2[w].sum()*abs(C*(HA['CDELT3']/1e9)/rest)*1e-3
    Ie=rms*np.sqrt(w.sum()*(1+r1)/(1-r1) if r1>0 else w.sum())*abs(C*(HA['CDELT3']/1e9)/rest)*1e-3
    print(f'    I(+-FWHM window)={I:.3f}+-{Ie:.3f} Jy/beam km/s')
    out[-1]['I']=I;out[-1]['eI']=Ie
json.dump(out,open(f'/workspace/out/s20_254_centroids_{ous}.json','w'),indent=1,default=float)
