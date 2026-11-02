"""Q44(a) summary of s260a: radial H13CO+/HCO+ and H13CN/HCN (2.2" beam, +-150 km/s) at N, S, Bridge: aperture ratios, annulus ratios (from area-weighted differences of nested aperture means; the
beam-smoothed core light is NOT removed, so annulus ratios are biased toward the core value = gradient underestimated), 8% cross-project scale added, and the N 0.3" vs 3.3" ratio difference.
usage: python3 -I s260f_q44_radial_summary.py -> out/q44/s260f.txt"""
import json,numpy as np
R=json.load(open('/workspace/out/q44/s260a.json')); L=[]
def P(s=''): print(s); L.append(s)
rads=[0.3,1.1,2.2,3.3]
for l,den in (('H13COp','HCOp'),('H13CN','HCN')):
    for rn in ('N','S','Bridge'):
        row=[]
        for r in rads:
            d=R[f'{r}|{rn}|{l}']; rat=d['amp']/d['ref']; e=np.hypot(d['err']/d['ref'],0.08*abs(rat)); row.append((rat,e,d))
        P(f'{l}/{den} {rn:6s}: '+'  '.join(f'r<{r}": {x[0]:.4f}+-{x[1]:.4f}' for r,x in zip(rads,row)))
        # annulus between 1.1 and 3.3
        a=[R[f'{r}|{rn}|{l}'] for r in (1.1,3.3)]; A=(3.3**2-1.1**2)
        num=(a[1]['amp']*3.3**2-a[0]['amp']*1.1**2)/A; ref=(a[1]['ref']*3.3**2-a[0]['ref']*1.1**2)/A
        en=np.hypot(a[1]['err']*3.3**2,a[0]['err']*1.1**2)/A
        c=[R[f'{r}|{rn}|{l}'] for r in (0.3,1.1)]  # core r<1.1 (nested)
        P(f'     ring 1.1-3.3": {num/ref:.4f}+-{en/ref:.4f} (H13C {num:.2f}+-{en:.2f}, 12C {ref:.1f}); core r<1.1": {c[1]["amp"]/c[1]["ref"]:.4f}')
for l in ('H13COp','H13CN'):
    a=R[f'0.3|N|{l}']; b=R[f'3.3|N|{l}']; ra=a['amp']/a['ref']; rb=b['amp']/b['ref']
    ea=a['err']/a['ref']; eb=b['err']/b['ref']
    P(f'N {l}: ratio(0.3") - ratio(3.3") = {ra-rb:+.4f} +- {np.hypot(ea,eb):.4f} ({(ra-rb)/np.hypot(ea,eb):.1f} sigma, independent-noise treatment; nested apertures overlap so this is conservative)  [8% cross-project scale cancels in a within-line comparison]; ratio factor {ra/rb:.2f}')
P('Harada+2018 Table 13 N (1.7", peak pixel): 0.0458+-0.0049; ours r<0.3" 2.2" beam 0.0372+-0.0059 (formal) -> consistent at 1.1 sigma; the +-150 km/s extraction is SiO-free.')
open('/workspace/out/q44/s260f.txt','w').write('\n'.join(L))
