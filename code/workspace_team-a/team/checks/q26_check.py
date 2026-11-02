import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from astropy.io import fits
import q26_measure as q
D=q.D
for u in (6,8):
  for o,J in q.SET if hasattr(q,'SET') else [('Xb55','2-1'),('Xb4f','3-2'),('Xb4d','3-2')]:
    nm=f'{o}_HCN_{J}_u{u}'
    im=fits.getdata(D+nm+'.image.fits'); h=fits.getheader(D+nm+'.image.fits')
    ps=fits.getdata(D+nm+'.psf.fits'); pb=fits.getdata(D+nm+'.pb.fits'); md=fits.getdata(D+nm+'.model.fits') if False else None
    ok=np.isfinite(im).all((1,2)); 
    print(nm, 'bmaj',h.get('BMAJ',0)*3600, 'chans with data',(np.nanmax(np.abs(im),(1,2))>0).sum(),'peak',np.nanmax(im),'min',np.nanmin(im),'psf max',np.nanmax(ps),'pbmin/max',np.nanmin(pb),np.nanmax(pb))
    s=np.nansum(im,(1,2))
    print('  sum spectrum mid chans', np.round(s[40:60],2))
