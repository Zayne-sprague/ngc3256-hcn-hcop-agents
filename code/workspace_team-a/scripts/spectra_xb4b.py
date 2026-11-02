# Xb4b full cubes: spectra at nuclear positions (peak pixel) and mean in 1" radius aperture (memmap, chunked)
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
NUC={'N':(156.96345,-43.90389),'S':(156.96340,-43.90533)}
for spw in ['spw25','spw27','spw29','spw31']:
    p=f('Xb4b',spw); h=fits.open(p,memmap=True); hd=h[0].header; w=WCS(hd).celestial
    cube=h[0].data[0]; fr=freqaxis(hd); out={'freq':fr}
    pix=hd['CDELT2']*3600
    for n,(ra,de) in NUC.items():
        x,y=[float(v) for v in w.world_to_pixel_values(ra,de)]; xi,yi=int(round(x)),int(round(y))
        r=int(np.ceil(1.0/abs(pix)))
        yy,xx=np.mgrid[-r:r+1,-r:r+1]; m=(xx**2+yy**2)*pix**2<=1.0
        sub=np.array(cube[:,yi-r:yi+r+1,xi-r:xi+r+1],float)
        out[n]=sub[:,r,r]; out[n+'_ap1']=np.nanmean(sub[:,m],axis=1)
    np.savez(f'/workspace/out/spec_Xb4b.{spw}.cube.I.pbcor.npz',**out); print(spw,'done',flush=True)
