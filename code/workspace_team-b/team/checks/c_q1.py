import numpy as np, glob
from astropy.io import fits
from astropy.wcs import WCS
C=299792.458
def load(ous,spw):
    f=glob.glob(f'/data/*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]
    h=fits.open(f,memmap=True); d=h[0].data; w=WCS(h[0].header)
    d=d.reshape(d.shape[-3:]); hd=h[0].header
    nu=hd['CRVAL3']+(np.arange(d.shape[0])-(hd['CRPIX3']-1))*hd['CDELT3']
    return d,w.celestial,nu,abs(hd['CDELT1'])*3600
regs={'N':(156.9645,-43.9040),'S':(156.9625,-43.9055)}
import sys
exec(open('/workspace/scripts/regions.py').read().split('\n',0)[0]) if False else None
sys.path.insert(0,'/workspace/scripts'); from regions import REGIONS
for name in ('N','S','Bridge'):
    ra,de=REGIONS[name][:2]; out={}
    for mol,spw,f0 in (('HCN','spw27',265.88618e9),('HCO+','spw25',267.55762e9)):
        d,w,nu,ps=load('Xb4d',spw); x,y=w.world_to_pixel(__import__('astropy.coordinates',fromlist=['x']).SkyCoord(ra,de,unit='deg'))
        x,y=int(round(float(x))),int(round(float(y)))
        sub=d[:,y-3:y+4,x-3:x+4]
        v=(1-nu/(f0/(1+0.0))*1)*C  # placeholder
        v=C*(1-nu/f0)  # radio velocity
        base=np.nanmedian(d[:, y-3:y+4,x-3:x+4][(np.abs(v-2780)>800)],axis=0)
        sel=np.abs(v-2780)<350
        # peak pixel in 7x7 box of integrated map
        m=np.nansum(sub[sel]-base,axis=0)*abs(nu[1]-nu[0])/nu.mean()*C
        out[mol]=m.max()
        out[mol+'c']=m[3,3]
    print(name,out,'R_peak',out['HCN']/out['HCO+'],'R_pos',out['HCNc']/out['HCO+c'])
