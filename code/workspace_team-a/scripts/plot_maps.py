import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
RA0,DE0=156.96346,-43.90461
def offs(w,shape,x0=0,y0=0):
    ny,nx=shape; X,Y=np.meshgrid(np.arange(nx)+x0,np.arange(ny)+y0)
    ra,de=w.pixel_to_world_values(X,Y); return (ra-RA0)*np.cos(np.deg2rad(DE0))*3600,(de-DE0)*3600
fig,axs=plt.subplots(2,3,figsize=(20,13))
hdb=fits.getheader(f('Xb4b','spw25')); wb=WCS(hdb).celestial
cont=np.squeeze(fits.getdata(f('Xb4b','spw25_27_29_31',kind='cont.I.pbcor')))
for a,(spw,nm) in zip(axs[0],[('spw25','HCO+ 3-2'),('spw27','HCN 3-2')]):
    d=np.load(f'/workspace/out/m0_Xb4b_{spw}.npz'); m=d['m0_core']
    ny,nx=m.shape; y0,x0=int(d['y0']),int(d['x0'])
    dx,dy=offs(wb,m.shape,x0,y0)
    a.pcolormesh(dx,dy,m,cmap='magma',vmax=np.nanpercentile(m,99.9)); a.set_title('Xb4b '+nm+' m0 (core)')
    c=cont[y0:y0+ny,x0:x0+nx]; a.contour(dx,dy,c,levels=[0.1e-3,0.3e-3,0.6e-3,1.2e-3,1.8e-3],colors='c',linewidths=.6)
    a.set_xlim(8,-8); a.set_ylim(-9,6); a.set_aspect('equal')
# Xb4d
hd=fits.getheader(f('Xb4d','spw25')); w=WCS(hd).celestial
cont4=np.squeeze(fits.getdata(f('Xb4d','spw25_27_29_31',kind='cont.I.pbcor')))
for a,(spw,nm) in zip(axs[1][:2],[('spw25','HCO+ 3-2'),('spw27','HCN 3-2')]):
    d=np.load(f'/workspace/out/m0_Xb4d_{spw}.npz'); m=d['m0']; dx,dy=offs(w,m.shape)
    a.pcolormesh(dx,dy,m,cmap='magma',vmax=np.nanpercentile(m,99.9)); a.set_title('Xb4d '+nm+' m0')
    a.contour(dx,dy,cont4,levels=[0.3e-3,0.6e-3,1.2e-3,2.4e-3,3.6e-3],colors='c',linewidths=.6)
    a.set_xlim(14,-14); a.set_ylim(-14,14); a.set_aspect('equal')
plt.tight_layout(); plt.savefig('/workspace/out/zoom_maps.png',dpi=60)
