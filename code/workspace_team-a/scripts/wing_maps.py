import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from plot_maps import offs
VSYS=2780.
fig,axs=plt.subplots(2,4,figsize=(24,12))
for r,(spw,nm,f0) in enumerate([('spw25','HCO+ 3-2',267.557626),('spw27','HCN 3-2',265.886434)]):
    h=fits.open(f('Xb4d',spw),memmap=True); hd=h[0].header; v=(1-freqaxis(hd)/1e9/f0)*C; dv=abs(np.median(np.diff(v))); w=WCS(hd).celestial
    cube=np.array(h[0].data[0],float); sig=np.nanstd(cube[:15],axis=0)
    dx,dy=offs(w,cube.shape[1:])
    for c,(lo,hi,t) in enumerate([(-450,-150,'blue -450..-150'),(-150,150,'core'),(150,450,'red 150..450'),(-450,450,'all')]):
        m=(v>=VSYS+lo)&(v<=VSYS+hi); m0=np.nansum(cube[m],axis=0)*dv; e=sig*dv*np.sqrt(m.sum())
        a=axs[r,c]; a.pcolormesh(dx,dy,m0/np.nanmedian(e),cmap='magma',vmin=-3,vmax=np.nanpercentile(m0/np.nanmedian(e),99.5)); a.set_title(f'Xb4d {nm} {t} (S/N units)')
        a.contour(dx,dy,m0/e,levels=[3,5],colors='c',linewidths=.5)
        a.set_xlim(12,-12); a.set_ylim(-12,12); a.set_aspect('equal')
plt.tight_layout(); plt.savefig('/workspace/out/wing_maps_xb4d.png',dpi=50)
