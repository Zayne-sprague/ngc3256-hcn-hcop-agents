# Q1: HCN/HCO+ 3-2 at the native 0.25" (Xb4b) and 0.94" (Xb4d) beams from the s5 moment-0 maps (|v-VSYS|<350, baseline-free window sums), same pixel for both lines
# (pixel of the nearest HCO+ peak within 0.3" (Xb4b) / 1" (Xb4d) of the region centre). Output out/s34_native_ratios.txt
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *; from regions import REGIONS
L=[]
for ous,sp,rad in (('Xb4b',0.047,0.3),('Xb4d',0.15,1.0)):
    z={m:np.load(f'/workspace/out/mom_{ous}.{s}_{m}32.npz') for m,s in (('HCO+','spw25'),('HCN','spw27'))}
    fn=files(f'*{ous}.ngc_3256_sci.spw25.cube.I.pbcor.fits')[0]; hd=fits.getheader(fn); wc=WCS(hd).celestial
    for reg in ('N','S','Bridge','E1','E2','W','Sblue'):
        px,py=wc.world_to_pixel_values(*REGIONS[reg][:2]); px-=int(z['HCO+']['x0']); py-=int(z['HCO+']['y0'])
        n=z['HCO+']['core'].shape[0]; Y,X=np.mgrid[:n,:n]; mk=np.hypot(X-px,Y-py)*sp*3600*(1 if ous=='Xb4b' else 1)<=rad
        # pixel scale in arcsec = CDELT
        ps=abs(hd['CDELT1'])*3600; mk=np.hypot(X-px,Y-py)*ps<=rad
        idx=np.unravel_index(np.argmax(np.where(mk,z['HCO+']['core'],-1e9)),(n,n))
        o,eo=z['HCO+']['core'][idx],z['HCO+']['core_err'][idx]; h,eh=z['HCN']['core'][idx],z['HCN']['core_err'][idx]
        r=h/o; e=r*np.hypot(eh/h,eo/o) if h>0 else np.nan
        L.append(f'{ous} {reg:7s} HCN {h:7.3f}+-{eh:.3f} HCO+ {o:7.3f}+-{eo:.3f} Jy/b km/s  R_Jy={r:.2f}+-{e:.2f} (stat; +5% cal)  S/N HCN {h/eh:.0f}')
open('/workspace/out/s34_native_ratios.txt','w').write('\n'.join(L)); print('\n'.join(L))
