import numpy as np, glob
from astropy.io import fits
from astropy.wcs import WCS
C=299792.458
def cube(spw):
    fn=glob.glob(f'/data/*Xb4d.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]; h=fits.open(fn); return h[0].data[0],h[0].header
RA0,DE0=156.96346,-43.90461
res={}
for k,spw,f0 in (('HCN',"spw27",265.886434),('HCO',"spw25",267.557626)):
    d,hd=cube(spw); n=hd['NAXIS3']; fr=(hd['CRVAL3']+(np.arange(n)-hd['CRPIX3']+1)*hd['CDELT3'])/1e9
    v=(1-fr/f0)*C; dv=abs(v[1]-v[0]); sel=(v>2430)&(v<3130)
    m=np.nansum(d[sel],0)*dv; res[k]=(m,hd)
    # noise: window 3600-4300 if within band
    print(k,v.min(),v.max(),sel.sum())
hd=res['HCN'][1]; w=WCS(hd).celestial
for name,(dra,dde) in {'N':None,'S':None}.items() if False else []: pass
m1=res['HCN'][0]; m2=res['HCO'][0]
for nm,box in (('N',(-1.5,0.5,0.5,2.5)),('S',(-1.5,0.5,-2,0))):
    pass
# find peaks of HCO+ in two halves split by dec offset
y,x=np.mgrid[:m2.shape[0],:m2.shape[1]]
ra,de=w.pixel_to_world_values(x,y)
dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600
r=np.hypot(dx,dy)
for nm,sel in (('north half',(dy>-0.5)&(r<5)),('south half',(dy<=-0.5)&(r<5))):
    i=np.nanargmax(np.where(sel,m2,-1e9)); print(nm,dx.flat[i],dy.flat[i],m1.flat[i],m2.flat[i],1.0126*m1.flat[i]/m2.flat[i])
