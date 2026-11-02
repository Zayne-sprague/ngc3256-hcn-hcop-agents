# Is the Xb4f feature near 251.83 GHz rest SO 6(5)-5(4) (251.8257) or the CH3OH A-type 3(3)-2(2)-ladder blend (251.74-251.92)? Template-shift test vs HCO+ 3-2.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import interp1d
from regions import REGIONS
C=299792.458; HCO=267.557626; REF=251.825770
def load(spw):
    h=fits.open(f'/data/member.uid___A001_X1273_Xb4f.ngc_3256_sci.spw{spw}.cube.I.pbcor.fits',memmap=True)
    H=h[0].header; d=np.squeeze(h[0].data); nu=(H['CRVAL3']+(np.arange(H['NAXIS3'])-H['CRPIX3']+1)*H['CDELT3'])/1e9; return d,H,nu
def sp(d,H,r):
    x,y=WCS(H).celestial.all_world2pix(*REGIONS[r][:2],0); return d[:,int(round(float(y))),int(round(float(x)))]*1e3
dH,HH,nuH=load('16'); dA,HA,nuA=load('22')
CH=[251.738437,251.811956,251.866524,251.890886,251.895728,251.900452,251.905729,251.917065]
w_=np.array([0.02,0.05,0.13,0.05,0.02,0.13,0.18,0.18]); w25=np.array([1.91,2.63,3.0,2.63,1.91,3.0,2.41,2.42])
cen_ch={'T=10K':(np.array(CH)*w_).sum()/w_.sum(),'T=25K':(np.array(CH)*w25).sum()/w25.sum()}
print('CH3OH A-ladder weighted centre (GHz):',cen_ch)
for r in ('S','Bridge','N'):
    sH=sp(dH,HH,r); vH=C*(1-nuH/HCO); o=((vH>1500)&(vH<2100))|((vH>3500)&(vH<4100)); m=(vH>2200)&(vH<3400)
    prof=sH-np.polyval(np.polyfit(vH[o],sH[o],1),vH); prof/=prof[m].max()
    v0=(prof[m].clip(0)*vH[m]).sum()/prof[m].clip(0).sum()
    tm=interp1d(vH-v0,prof,bounds_error=False,fill_value=0.)
    s=sp(dA,HA,r); v=C*(1-nuA/REF); sel=(v>1000)&(v<4600); line=(v>2250)&(v<3350)
    pf=np.polyfit(v[sel&~line],s[sel&~line],1); res=s[sel&~line]-np.polyval(pf,v[sel&~line]); rms=res.std(); r1=np.corrcoef(res[:-1],res[1:])[0,1]; infl=(1+r1)/(1-r1)
    s2=s-np.polyval(pf,v); fw_=(v>1800)&(v<3700); dvs=np.arange(-300,300.1,2.); chi=[];am=[]
    for dv in dvs:
        t=tm(v[fw_]-(v0+dv)); A=(t*s2[fw_]).sum()/(t*t).sum(); am.append(A); chi.append((((s2[fw_]-A*t)/rms)**2).sum())
    chi=np.array(chi); k=chi.argmin(); ok=(chi-chi[k])/infl<=1; chi0=((s2[fw_]/rms)**2).sum()
    dv=dvs[k]; e=max((dvs[ok].max()-dvs[ok].min())/2,2.)
    print(f'{r}: rms {rms:.2f} mJy/ch lag1 {r1:.2f}; template amp {am[k]:.2f} mJy/beam; dchi2 {(chi0-chi[k])/infl:.1f}; centre vs HCO+ profile when assigned to SO 251.8257: {dv:+.0f}+-{e:.0f} km/s')
    for nm,f in (('SO 6(5)-5(4) 251.8257',REF),('CH3OH A-ladder (T=10K weights) %.4f'%cen_ch['T=10K'],cen_ch['T=10K']),('CH3OH A-ladder (T=25K) %.4f'%cen_ch['T=25K'],cen_ch['T=25K']),('CH3OH 5(3,2)-5(2,3) A 251.8120',251.811956)):
        vv=C*(1-REF*(1-(v0+dv)/C)/f)-v0
        print(f'    {nm:42s} offset from HCO+ = {vv:+.0f}+-{e:.0f} km/s ({vv/e:+.1f} sigma)')
