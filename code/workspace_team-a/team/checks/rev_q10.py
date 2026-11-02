import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
from q7_lib import *
from astropy.io import fits
for ous in ('Xb4b',):
  z0=np.load(f'/workspace/out/q7/cube_{ous}_HCN.npz'); x0,y0=int(z0['x0']),int(z0['y0']); ny,nx=z0['cube'].shape[1:]
  d=fits.getdata(f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')[0,0][y0:y0+ny,x0:x0+nx]
  for ln in ('HCN','HCO+'):
    z=np.load(f'/workspace/out/q7/cube_{ous}_{ln}.npz'); cube=z['cube'].astype(float); v=z['v']-VSYS
    X=Cube(ous,ln)
    for reg,(rx,ry) in {'N':(0,2.6),'S':(-0.15,-2.6)}.items():
      m=np.hypot(X.X-rx,X.Y-ry)<0.2
      sp=cube[:,m].mean(1)*1e3
      cont=d[m].mean()*1e3
      base=np.abs(v)>520
      print(ln,reg,'pipeline cont',round(cont,3),'mJy; v range',v.min(),v.max())
      # print spectrum binned 100 km/s
      for a in range(-1000,1000,100):
        k=(v>=a)&(v<a+100)
        if k.any(): print('  ',a,round(sp[k].mean(),2),end=';')
      print()
      # baseline fit using only |v|>800 and compare
      for lim in (520,700,900):
        k=np.abs(v)>lim
        if k.sum()>4: 
          A=np.vstack([np.ones(k.sum()),v[k]/500]).T; c=np.linalg.lstsq(A,sp[k],rcond=None)[0]; print('   fit |v|>',lim,'c0=',round(c[0],3),'slope',round(c[1],3))
