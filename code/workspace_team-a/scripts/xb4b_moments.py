# Xb4b (0.24") HCO+ 3-2 / HCN 3-2 (and other windows) moment-0 maps, cropped to central +-12"
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
VSYS=2780.
LINES={'spw25':('HCO+ 3-2',267.557626),'spw27':('HCN 3-2',265.886434)}
R=255
for spw,(nm,f0) in LINES.items():
    h=fits.open(f('Xb4b',spw),memmap=True); hd=h[0].header; fr=freqaxis(hd); v=(1-fr/1e9/f0)*C; dv=abs(np.median(np.diff(v)))
    cx,cy=int(hd['CRPIX1']-1),int(hd['CRPIX2']-1)
    sl=(slice(cy-R,cy+R),slice(cx-R,cx+R))
    out={}
    for tag,(lo,hi) in {'core':(-350,350),'blue':(-350,-150),'red':(150,350)}.items():
        sel=np.where((v>=VSYS+lo)&(v<=VSYS+hi))[0]
        c=np.array(h[0].data[0,sel.min():sel.max()+1][(slice(None),)+sl],float)
        out['m0_'+tag]=np.nansum(c,axis=0)*dv; out['n_'+tag]=len(sel)
    off=np.array(h[0].data[0,0:20][(slice(None),)+sl],float); out['sig']=np.nanstd(off,axis=0)
    out['dv']=dv; out['x0']=cx-R; out['y0']=cy-R
    np.savez(f'/workspace/out/m0_Xb4b_{spw}.npz',**out); print(spw,nm,'done',dv,flush=True)
