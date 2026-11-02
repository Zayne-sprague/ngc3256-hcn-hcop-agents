# Q36 step 2: free-free + dust decomposition of the nuclear continuum (matched 10.5" and 6.5" brightness) and FIR check.
# ff: S~nu^-0.1 ; dust: S~nu^(2+beta); fits per nucleus; compared with Michiyama+2020 (ApJ 895,85) Fig.7 / Table note: ff fraction at 99 GHz 76% (N), 90% (S) (2.5"/3.5" apertures; dust index fixed 4).
import json,numpy as np
R=json.load(open('/workspace/out/q36/q36_cont.json')); out={}; L=[]
nus=np.array([170.7,257.7,346.5])
for tgt in ('10.5','6.5'):
  for r in 'NS':
    ks=[k for k in ('B5','B6f','B7') if f'{k}_{tgt}_{r}' in R]; nu=np.array([R[f'{k}_{tgt}_{r}']['nu'] for k in ks]); S=np.array([R[f'{k}_{tgt}_{r}']['mJy'] for k in ks])
    # 10% flux-scale error per band (ALMA) + 3% noise
    sig=0.10*S
    for beta in (1.5,2.0):
        best=None
        for ff in np.linspace(0,1.5*S[0] if tgt=='10.5' else 1.5*S[0],4000):
            a=ff*(nu/257.7)**-0.1; X=(nu/257.7)**(2+beta); 
            d=np.sum((S-a)*X/sig**2)/np.sum(X*X/sig**2); d=max(d,0); m=a+d*X; c=np.sum(((S-m)/sig)**2)
            if best is None or c<best[0]: best=(c,ff,d)
        c,ff,d=best; fd=lambda v: d*(v/257.7)**(2+beta)/(ff*(v/257.7)**-0.1+d*(v/257.7)**(2+beta))
        key=f'{r}_{tgt}_b{beta}'
        out[key]=dict(chi2=c,ff258_mJy=ff,dust258_mJy=d,fdust_258=fd(257.7),fdust_346=fd(346.5),fdust_170=fd(170.7),fdust_99=fd(99.0),ff_frac_99=1-fd(99.0))
        L.append(f'{r} {tgt}" beta={beta}: chi2={c:.2f}(1 dof)  S_ff(258)={ff:.1f} mJy, S_dust(258)={d:.1f} mJy; dust fraction 170/258/346 GHz = {fd(170.7):.2f}/{fd(257.7):.2f}/{fd(346.5):.2f}; implied ff fraction at 99 GHz = {1-fd(99.0):.2f}  (Michiyama+20: N 0.76, S 0.90)')
# FIR: global modified BB (S_nu ~ nu^beta B_nu(T)) from the Michiyama+2020 flux densities (S70,S100,S160,S250 = 120.3,145.4,93.48,93.48 Jy; S250 is printed equal to S160, probably a typo -> 250um excluded)
from scipy.optimize import brentq
def ratio(T,beta,l1,l2): n1,n2=3e14/l1,3e14/l2; return (n1/n2)**(3+beta)*np.expm1(4.799e-11*n2/T)/np.expm1(4.799e-11*n1/T)
for beta in (1.5,2.0):
    for (l1,l2,s1,s2) in ((100,160,145.4,93.48),(70,160,120.3,93.48)):
        T=brentq(lambda T: ratio(T,beta,l1,l2)-s1/s2,10,200); L.append(f'global FIR S{l1}/S{l2} -> Td={T:.1f} K (beta={beta}); (70um also has a hotter warm-dust component)'); out[f'FIR_{l1}_{l2}_b{beta}']=T
print('\n'.join(L)); open('/workspace/out/q36/q36_sed.log','w').write('\n'.join(L)); json.dump(out,open('/workspace/out/q36/q36_sed.json','w'),indent=1)
