import numpy as np
from astropy.io import fits
d='/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw%d.cube.I.pbcor.fits'
c=299792.458
out={}
for sp,rest in((25,267.557626),(27,265.886434)):
    h=fits.open(d%sp,memmap=True); hd=h[0].header
    D=h[0].data  # (1,nchan,y,x) or (nchan,y,x)
    D=D.reshape(D.shape[-3:])
    nu=hd['CRVAL3']+(np.arange(D.shape[0])-(hd['CRPIX3']-1))*hd['CDELT3']
    v=c*(rest*1e9/nu-1) if False else c*(rest*1e9-nu)/(rest*1e9)  # radio
    v=c*(rest*1e9-nu)/(rest*1e9)
    print(sp,hd['CRVAL3'],hd['CDELT3'],v[0],v[-1],D.shape,hd['CDELT2']*3600)
    ch=np.where(abs(v-2780)<350)[0]
    lf=np.where((abs(v-2780)>700)&(abs(v-2780)<1100))[0]
    # restrict to central 300x300 pix
    sl=(slice(None),slice(270,570),slice(270,570))
    X=np.array(D[:,270:570,270:570],dtype=np.float32)
    X=np.nan_to_num(X)
    base=np.median(X[lf],axis=0)
    m0=(X[ch]-base).sum(0)*abs(np.mean(np.diff(v)))
    out[sp]=m0
    h.close()
np.savez('/workspace/team/checks/chk_q26_m0.npz',a=out[25],b=out[27])
