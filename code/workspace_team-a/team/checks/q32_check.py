import json,numpy as np
from scipy.optimize import curve_fit
d=json.load(open('out/q23/co_hires.json'));k=['native','0.26','0.38','0.57','0.75','1.0','1.7']
th=np.array([0.178,.26,.38,.57,.75,1.,1.7])
for c in 'SN':
    T=np.array([d['cores'][c]['beams'][x]['Tb_RJ'] for x in k]);print(c,T.round(1))
    print(' flux ratio (1.0 vs 0.26 aperture-equivalent Tb*th^2):',T[5]*1/(T[1]*.26**2))
    # pedestal-free power-law-in-size alt: two-Gaussian core+extended with free extended size
    m=lambda t,a,s1,b,s2:a*s1**2/(s1**2+t**2)+b*s2**2/(s2**2+t**2)
    try:
        p,_=curve_fit(m,th,T,p0=[30,.4,5,3],bounds=([0,.05,0,1.5],[300,3,50,30]),maxfev=20000);print(' 2G',p.round(2),'ff26',p[1]**2/(p[1]**2+.26**2))
    except Exception as e:print(e)
    for s in (0.75,1.37): print(' nopedestal ff(0.26) s=',s,s**2/(s**2+.26**2))
