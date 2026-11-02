# Q42 part B (model-free): 13C/12C ratios at ONE common 2.2" circular taper (s185 = out/q29/s185.json/.txt; X961 12C + X95d 13C, 2015.1.00993 short-spacing complete; the native beam is
# 2.1x1.7" so 1.7" itself cannot be reached) for N, S, C(Bridge), three extractions, and the N/S/C contrasts of the thin-limit 12C/13C = 1/R (MC with gaussian errors; negative/low R of S kept as is).
# Harada+2018 Table 13 (1.7", peak pixel, +-150 km/s) H13CO+/HCO+ N = 9.40/205.25 and C = 3.44/139.46 as comparison.
# usage: python3 -I s250e_q42_common_taper_ratios.py -> out/q42/common_taper_ratios.txt
import numpy as np
rng=np.random.default_rng(1); NS=200000
D={ # (R, err)  s185.txt / s189 / report entries 200-201
 'H13CO+/HCO+  template aperture (r<2.2", |v|<500)':{'N':(0.0231,0.0031),'S':(0.0184,0.0040),'C':(0.0185,0.0033)},
 'H13CO+/HCO+  aperture +-150 km/s':{'N':(0.0271,0.0049),'S':(0.0150,0.0047),'C':(0.0182,0.0045)},
 'H13CO+/HCO+  peak pixel +-150 km/s':{'N':(0.0372,0.0066),'S':(0.0141,0.0088),'C':(0.0222,0.0087)},
 'H13CN/HCN    template aperture':{'N':(0.0207,0.0034),'S':(0.0057,0.0059),'C':(0.0237,0.0044)},
 'H13CN/HCN    aperture +-150 km/s':{'N':(0.0220,0.0046),'S':(0.0048,0.0069),'C':(0.0234,0.0058)},
 'H13CN/HCN    peak pixel +-150 km/s':{'N':(0.0302,0.0058),'S':(-0.0006,0.0130),'C':(0.0265,0.0117)},
}
L=[]
def P(x=''): print(x); L.append(x)
P('Harada+2018 Table 13 (1.7", peak pixel): H13CO+/HCO+ N 0.0458+-0.0049 (HCO+/H13CO+ = 22+-3), C 0.0247+-0.0069 (40 +-11), S not tabulated.')
for k,v in D.items():
    P(f'\n{k}')
    s={p:rng.normal(*v[p],NS) for p in v}
    for p in v:
        R,e=v[p]; thin=1/R if R>0 else np.inf
        lo=1/(R+e) if R+e>0 else np.inf; hi=1/(R-e) if R-e>0 else np.inf
        lim3=1/(R+3*e) if R+3*e>0 else np.inf
        P(f'  {p}: R = {R:.4f}+-{e:.4f} ({R/e:.1f} sigma) -> thin 12C/13C = {thin:.0f} [{lo:.0f}, {hi if np.isfinite(hi) else float("inf"):.0f}]; 3-sigma lower limit on 12C/13C > {lim3:.0f}')
    for a,b in (('N','S'),('N','C'),('C','S')):
        d=s[a]-s[b]; rr=np.where(s[b]>0,s[a]/np.where(s[b]>0,s[b],1),np.inf)
        P(f'  R_{a}-R_{b} = {d.mean():+.4f}+-{d.std():.4f} ({d.mean()/d.std():+.1f} sigma); P(R_{a}>R_{b}) = {np.mean(d>0):.2f}; ratio R_{a}/R_{b} median {np.median(rr):.2f} [16-84%: {np.percentile(rr,16):.2f}, {np.percentile(rr,84):.2f}]  => 12C/13C(thin) S/N-type contrast X_{b}/X_{a} = {np.median(rr):.2f}')
# Harada N (1.7") vs our 2.2" at the same N peak pixel and aperture
h=(0.0458,0.0049); o=(0.0372,0.0066)
P(f'\nBeam/aperture dependence at N: Harada 1.7" peak {h[0]:.4f}+-{h[1]:.4f}; ours 2.2" peak {o[0]:.4f}+-{o[1]:.4f} ({(h[0]-o[0])/np.hypot(h[1],o[1]):.1f} sigma diff.); 2.2" aperture +-150 0.0271+-0.0049 ({(h[0]-0.0271)/np.hypot(h[1],0.0049):.1f} sigma below Harada); template aperture 0.0231+-0.0031 ({(h[0]-0.0231)/np.hypot(h[1],0.0031):.1f} sigma). The N value falls by 2x from the 1.7" peak pixel to the 2.2" aperture, i.e. 12C/13C(thin) rises 22 -> 43 within N alone.')
open('/workspace/out/q42/common_taper_ratios.txt','w').write('\n'.join(L))
