# extract spectra at nuclei positions for all cubes (pbcor), save npz + find peaks
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
import glob, json
NUC={'N':(156.96345,-43.90389),'S':(156.96340,-43.90533)}
out={}
for p in sorted(glob.glob(D+'*cube.I.pbcor.fits')+glob.glob(D+'*cube.I.manual.image.pbcor.fits')+glob.glob(D+'*tapered.image.pbcor.fits')):
    if 'Xb4b' in p and 'manual' not in p: continue  # big; handled separately
    if 'cont' in p: continue
    key=os.path.basename(p).replace('member.uid___A001_X1273_','').replace('.ngc_3256_sci','').replace('_ngc_3256_sci','_tap').replace('.fits','')
    h=fits.open(p,memmap=True); hd=h[0].header; w=WCS(hd).celestial
    cube=h[0].data[0] if h[0].data.ndim==4 else h[0].data
    fr=freqaxis(hd)
    spec={}
    for n,(ra,de) in NUC.items():
        x,y=w.world_to_pixel_values(ra,de); x=int(round(float(x))); y=int(round(float(y)))
        spec[n]=np.array(cube[:,y,x],float)
    # also peak-pixel spectrum within 3" of centre of the line-integrated map
    np.savez(f'/workspace/out/spec_{key}.npz',freq=fr,**spec)
    print(key,cube.shape, {n:f'{np.nanstd(s)*1e3:.2f}' for n,s in spec.items()})
