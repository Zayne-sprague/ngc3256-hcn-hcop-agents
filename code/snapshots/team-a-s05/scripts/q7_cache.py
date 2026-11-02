# Q7: cache central-crop (+-HALF arcsec about the N-S midpoint) pbcor cubes of HCN/HCO+ 3-2 (Xb4b, Xb4d) as float32 npz (Jy/beam). python3 -I q7_cache.py
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
RA0,DE0=156.96346,-43.90461; HALF=7.0
K={'HCN':('spw27',265.886434),'HCO+':('spw25',267.557626)}
for ous in ('Xb4b','Xb4d'):
    for ln,(spw,f0) in K.items():
        h=fits.open(f(ous,spw),memmap=True); hd=h[0].header; w=WCS(hd).celestial
        cx,cy=[float(a) for a in w.world_to_pixel_values(RA0,DE0)]; pix=abs(hd['CDELT2'])*3600; hp=int(HALF/pix)
        x0,y0=int(round(cx))-hp,int(round(cy))-hp
        a=np.array(h[0].data[0,:,y0:y0+2*hp+1,x0:x0+2*hp+1],np.float32)
        fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C
        # per-channel beam (BEAMS table) if present
        nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
        np.savez(f'/workspace/out/q7/cube_{ous}_{ln}.npz',cube=a,v=v,x0=x0,y0=y0,pix=pix,beam=nb,cd1=hd['CDELT1'],cd2=hd['CDELT2'],
                 ra=np.array(w.pixel_to_world_values(np.arange(a.shape[2])+x0,np.zeros(a.shape[2])+y0)[0]),
                 de=np.array(w.pixel_to_world_values(np.zeros(a.shape[1])+x0,np.arange(a.shape[1])+y0)[1]))
        print(ous,ln,a.shape,nb,'nan frac',np.isnan(a).mean(),flush=True)
