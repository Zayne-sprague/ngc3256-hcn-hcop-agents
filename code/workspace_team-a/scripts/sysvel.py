import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
C=299792.458
LINES={'Xb4f.spw16':('HCO+ 3-2',267.557626),'Xb4f.spw18':('HCN 3-2',265.886434),'Xb47.spw16':('HCN 4-3',354.505477),'Xb47.spw18':('HCO+ 4-3',356.734223),
 'Xb47.spw22':('CS 7-6',342.882850),'Xb55.spw16':('HCN 2-1',177.261111),'Xb55.spw18':('HCO+ 2-1',178.375065),'Xb4d.spw25':('HCO+ 3-2',267.557626),'Xb4d.spw27':('HCN 3-2',265.886434)}
for k,(nm,f0) in LINES.items():
    d=np.load(f'/workspace/out/spec_{k}.cube.I.pbcor.npz'); fr=d['freq']/1e9; v=(1-fr/f0)*C
    for n in 'NS':
        s=d[n]; m=(v>2780-450)&(v<2780+450)
        # robust window: iterate
        for _ in range(3):
            thr=0.05*s[m].max(); mm=m&(s>thr); lo,hi=v[mm].min(),v[mm].max()
        w=(v>=lo-30)&(v<=hi+30); 
        m1=np.sum(v[w]*s[w])/np.sum(s[w]); I=np.sum(s[w])*abs(np.median(np.diff(v)))
        print(f'{k} {nm:9s} {n}: v1={m1:7.1f} km/s  range {lo:.0f}..{hi:.0f}  peak={s[m].max()*1e3:.1f} mJy/b  I={I:.2f} Jy/b km/s  vmid={(lo+hi)/2:.0f}')
