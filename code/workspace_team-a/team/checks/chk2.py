import numpy as np
from astropy.io import fits
c=299792.458
D='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.%s.cube.I.pbcor.fits'
pos={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15)}
out={}
for spw,f0,nm in [('spw25',267.557626e9,'HCO+'),('spw27',265.886434e9,'HCN')]:
    h=fits.open(D%spw,memmap=True);d=h[0].data;d=d[0] if d.ndim==4 else d;H=h[0].header
    f=H['CRVAL3']+(np.arange(d.shape[0])-H['CRPIX3']+1)*H['CDELT3']
    v=c*(1-f/f0); dv=abs(np.diff(v).mean())
    for r,(dx,dy) in pos.items():
        x=int(round(H['CRPIX1']-1+dx/H['CDELT1']/3600*-1*-1)) if False else int(round(135+dx/(-0.15)))
        y=int(round(135+dy/0.15))
        s=d[:,y-1:y+2,x-1:x+2].max(axis=(1,2))  # peak in 3x3
        s0=d[:,y,x]
        # find peak pixel
        for cen,tag in [(2780,'fixed2780'),(None,'centred')]:
            if cen is None:
                cen=(s0*v).sum()/s0.sum() if False else v[np.argmax(s0)]
            m=abs(v-cen)<350
            I=s0[m].sum()*dv/1e3 if False else s0[m].sum()*dv
            out[(nm,r,tag)]=I
        rms=np.nanstd(s0[(v<2780-900)|(v>2780+900)])
        print(nm,r,'dv',round(dv,2),'peakJy',s0.max(),'vpk',v[np.argmax(s0)],'I2780',round(out[(nm,r,'fixed2780')],2),'Icen',round(out[(nm,r,'centred')],2),'rms mJy',rms*1e3)
for r in pos: print(r,'ratio fixed',out[('HCN',r,'fixed2780')]/out[('HCO+',r,'fixed2780')],'centred',out[('HCN',r,'centred')]/out[('HCO+',r,'centred')])
