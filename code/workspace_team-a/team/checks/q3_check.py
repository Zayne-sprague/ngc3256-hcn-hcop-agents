import numpy as np,sys
from astropy.io import fits
from astropy.wcs import WCS
C=299792.458
RA0,DE0=156.96346,-43.90461
def run(spw,f0,tag):
    fn=f'/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.{spw}.cube.I.pbcor.fits'
    h=fits.open(fn,memmap=True); hd=h[0].header; d=h[0].data[0]
    w=WCS(hd).celestial
    fr=(hd['CRVAL3']+(np.arange(d.shape[0])-(hd['CRPIX3']-1))*hd['CDELT3'])/1e9
    v=(1-fr/f0)*C; dv=abs(np.median(np.diff(v)))
    out={}
    for name,(x,y) in {'Soutflow':(-0.21,-4.48),'ctrl_E':(4.0,-4.48),'ctrl_W':(-4.4,-4.48)}.items():
        px,py=w.world_to_pixel_values(RA0+x/3600/np.cos(np.deg2rad(DE0)),DE0+y/3600)
        px,py=int(round(float(px))),int(round(float(py)))
        r=int(round(0.5/0.15)); s=np.nan_to_num(d[:,py-r:py+r+1,px-r:px+r+1]).mean((1,2))*1e3  # mJy/beam, ~1" box mean
        out[name]=s
    s=out['Soutflow']
    free=(v>1400)&(v<2100); sd=1.4826*np.median(np.abs(s[free]-np.median(s[free])))
    print(tag,'chan rms (mJy)',round(sd,2))
    for lo in range(2175,3000,100):
        m=(v>=lo)&(v<lo+100); print(tag,lo,round(s[m].sum()*dv/1e3,3),'Jy/b km/s', 'ctrlE',round(out['ctrl_E'][m].sum()*dv/1e3,3),'ctrlW',round(out['ctrl_W'][m].sum()*dv/1e3,3))
run('spw27',265.88618,'HCN'); run('spw25',267.55762,'HCO+')
