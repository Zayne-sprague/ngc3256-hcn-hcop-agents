# Q13e: is the blue-lobe emission inside the Xb4d clean mask?  (uncleaned channels: flux scale is dirty-beam, not clean-beam)
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
C=299792.458
for sp,spw in (('HCN','spw27'),('HCO+','spw25')):
    h=fits.open(f('Xb4d',spw,'cube.I.mask','.'),memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9; v=vel(fr,F0[sp])
    d=h[0].data; d=d[0] if d.ndim==4 else d
    cx,cy=[float(a) for a in WCS(hd).celestial.world_to_pixel_values(RA0,DE0)]
    x,y=POS['Soutflow']; pix=abs(hd['CDELT2'])*3600
    px=int(round(cx+x/pix*-1)) ; py=int(round(cy+y/pix))
    yy,xx=np.mgrid[:d.shape[1],:d.shape[2]]
    off=np.hypot((xx-cx)*pix*-1-x,(yy-cy)*pix-y)<=1.0
    fr_=np.array([(d[i][off]>0.5).mean() for i in range(len(v))])
    for a,b in ((2375,2575),(2575,2875),(2300,2375),(2875,3100)):
        k=(v>=a)&(v<b); print(sp,f'{a}-{b}: mask fraction of r=1" lobe aperture, mean over channels {fr_[k].mean():.3f}, max {fr_[k].max():.3f}')
    nz=np.where(d.reshape(len(v),-1).max(1)>0.5)[0]; print(sp,'mask channels span v =',v[nz].min(),v[nz].max())
