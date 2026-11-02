import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
c=299792.458
def spec(spw,ra,dec,ous='Xb4f'):
    h=fits.open(f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.spw{spw}.cube.I.pbcor.fits',memmap=True)
    H=h[0].header;d=h[0].data[0] if h[0].data.ndim==4 else h[0].data
    nu=(H['CRVAL3']+(np.arange(H['NAXIS3'])-H['CRPIX3']+1)*H['CDELT3'])/1e9
    x,y=WCS(H).celestial.all_world2pix(ra,dec,0);return nu,d[:,int(round(float(y))),int(round(float(x)))]*1e3
S=(156.96340,-43.90533);N=(156.96346,-43.90390)
for nm,p in (('S',S),('N',N)):
  for spw,tag in (('20','254.0 region'),):
    nu,s=spec(spw,*p)
    # HCO+/HCN 3-2 centroid in spw16/18 etc not needed; use velocity relative to 254.0153
    for rest,lab in ((254.01534,'CH3OH 2_0-1_-1E'),(253.96800,'NS 11/2-9/2 f')):
        v=(1-nu*1e0/ (rest))*c  # radio velocity w.r.t. observed freq nu
        # radio: v=c(1-nu/nu_rest*(1+z))... compute v_radio = c*(rest-nu)/rest  -> that's rest-frame freq w/o z; we use nu_obs: v=c*(1-nu/rest) is wrong sign conv; do properly
        v=c*(rest-nu)/rest  # frequency of obs vs rest
        # for z=0.0093, nu_obs ~ rest/(1+z): v_opt-ish ; use radio: v=c*(rest-nu)/rest
        m=(v>2400)&(v<3150)
        w=np.clip(s[m],0,None)
        # baseline from outer
        off=((v>1500)&(v<2300))|((v>3300)&(v<4200))
        b=np.median(s[off]);sd=s[off].std()
        ss=s[m]-b
        cen=(ss*v[m]).sum()/ss.sum()
        print(nm,lab,'centroid v=%.0f'%cen,'peak v=%.0f'%v[m][np.argmax(ss)],'peak %.1f rms %.1f'%(ss.max(),sd),'I=%.2f'%(ss.sum()*3.91e-3*c/ (nu.mean())*1e0/1e3*1e3/1e3))
