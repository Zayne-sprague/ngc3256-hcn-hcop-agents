"""r11: N/S HCN/HCO+ 3-2 contrast vs estimator (peak pixel, apertures r=0.25,0.5,0.75,beam) on the s192 uv-tapered cubes
(/workspace/out/q34/*_{1p5,2p2}_{comb,c12}.pbcor.fits, made by scripts/s191_q34_casa_image.py).  Same conventions as scripts/s192_q34_measure.py:
centre = HCO+ 3-2 +-350 km/s moment-0 peak within 0.4" of nominal, +-350 km/s window, baseline from line-free channels (v<-420,>440), RJ T_b.
Cell is 0.3": apertures are evaluated on a 0.05" grid by cubic-spline interpolation of each channel (beam >= 5 pixels, so interpolation is safe).
Usage: python3 -I r11_apertures.py"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from regions import REGIONS
O='/workspace/out/q34/'; C=299792.458; V0=2780.
F0={'HCN32':265.88618,'HCOp32':267.55762}
def load(tag):
    h=fits.open(O+tag+'.pbcor.fits'); d=np.nan_to_num(h[0].data.astype(float)).squeeze(); hd=h[0].header; w=WCS(hd).celestial
    nu=(hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3'])/1e9
    v=C*(1-nu/F0[tag.split('_')[0]])-V0 if hd['CUNIT3'].strip()=='Hz' else hd['CRVAL3']/1e3+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']/1e3-V0
    bm=hd['BMAJ']*3600; pix=abs(hd['CDELT2'])*3600
    Tb=d*1.222e6/(F0[tag.split('_')[0]]**2*bm**2); lf=(v<-420)|(v>440)
    return Tb-np.median(Tb[lf],axis=0),v,pix,w,lf
def kw(v,dv): return np.clip(np.minimum(v+dv/2,350)-np.maximum(v-dv/2,-350),0,dv)/dv
SUB=0.05
res={}
for bn,th in (('1p5',1.5),('2p2',2.2)):
  for vn in ('comb','c12'):
    D={l:load(f'{l}_{bn}_{vn}') for l in F0}
    Tb,v,pix,w,lf=D['HCOp32']; ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]; dv=abs(v[1]-v[0]); k=kw(v,dv)
    m0=np.tensordot(k,Tb,axes=(0,0))*dv
    for rn in ('N','S'):
        px,py=[float(t) for t in w.world_to_pixel_values(*REGIONS[rn][:2])]
        mm=np.where(np.hypot(X-px,Y-py)*pix<0.4,m0,-np.inf); j=np.argmax(mm); cx,cy=float(X.flat[j]),float(Y.flat[j])
        out={}
        g=np.arange(-th*1.2,th*1.2+1e-9,SUB)/pix; gy,gx=np.meshgrid(g,g,indexing='ij'); rr=np.hypot(gy,gx)*pix
        for ln in F0:
            T=D[ln][0]
            sp_pk=T[:,int(cy),int(cx)]
            cube=np.stack([map_coordinates(T[i],[gy+cy,gx+cx],order=3) for i in range(T.shape[0])])
            for lab,sp in [('peak',sp_pk)]+[(f'r{r}',cube[:,rr<=r].mean(axis=1)) for r in (0.25,0.5,0.75,th)]:
                sig=float(np.std(sp[lf])); kk=k
                out[f'{ln}|{lab}']=[float((kk*sp).sum()*dv),float(sig*np.sqrt((kk**2).sum())*dv)]
        res[f'{bn}|{vn}|{rn}']=out
import itertools
for k_,v_ in res.items():
    for a_,b_ in v_.items():
        if b_[0]==0: print("ZERO",k_,a_)
labs=['peak','r0.25','r0.5','r0.75','rbeam']
tab={}; lines=[]
for bn,th in (('1p5',1.5),('2p2',2.2)):
  for vn in ('comb','c12'):
    for lab in labs:
      key=lambda ln,rn: res[f'{bn}|{vn}|{rn}'][f'{ln}|'+(f'r{th}' if lab=='rbeam' else lab)]
      row={}
      for rn in ('N','S'):
        h,o=key('HCN32',rn),key('HCOp32',rn); R=h[0]/o[0]; row[rn]=(R,R*np.hypot(h[1]/h[0],o[1]/o[0]))
      c=row['N'][0]/row['S'][0]; ce=c*np.hypot(row['N'][1]/row['N'][0],row['S'][1]/row['S'][0])
      tab[f'{bn}|{vn}|{lab}']=dict(N=row['N'],S=row['S'],NS=(c,ce))
      lines.append(f'{th}\" {vn:4s} {lab:6s} N {row["N"][0]:.3f}+-{row["N"][1]:.3f}  S {row["S"][0]:.3f}+-{row["S"][1]:.3f}  N/S {c:.3f}+-{ce:.3f}')
open('r11_apertures.txt','w').write('\n'.join(lines)); print('\n'.join(lines))
json.dump(dict(res=res,tab=tab),open('r11_apertures.json','w'))
