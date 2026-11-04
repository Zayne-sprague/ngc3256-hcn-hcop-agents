"""r12: J-independence D = C(3-2)/C(1-0), C = N/S HCN/HCO+ T_b-ratio contrast, same estimator and +-350 km/s window.
3-2 from r11_apertures.json (s192 uv-tapered cubes, 12m+7m 'comb' and 12m-only 'c12'; beams 1.5, 2.2); 1-0 from r12_apertures_10.json (X145 1.5, 2.2; X961 2.2)."""
import json, numpy as np
a=json.load(open('../r11/r11_apertures.json'))['tab']; b=json.load(open('r12_apertures_10.json'))
out=[]
for ds,th,bn in (('X145',1.5,'1p5'),('X145',2.2,'2p2'),('X961',2.2,'2p2')):
  for vn in ('comb','c12'):
    for lab in ('peak','r0.25','r0.5','r0.75','rbeam'):
      c32,e32=a[f'{bn}|{vn}|{lab}']['NS']; c10,e10=b[f'NS|{ds}|{th}|w350|{lab}'][:2]
      D=c32/c10; out.append(f'{ds} {th}" 3-2 {vn:4s} {lab:6s} C32 {c32:.3f}+-{e32:.3f} C10 {c10:.3f}+-{e10:.3f}  D {D:.3f}+-{D*np.hypot(e32/c32,e10/c10):.3f}')
open('r12_D.txt','w').write('\n'.join(out)); print('\n'.join(out))
