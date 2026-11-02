# Q13d: window/aperture/array dependence of the lobe HCN/HCO+ 3-2 (Tb ratio = flux ratio, same beam) -- reconciles the 1.0-1.4 (fixed blue window), 1.7-1.8 (E window, r<=1") and 0.7-0.8 (whole broad component).
# Reads out/q13/q13a_maps_*.npz.  Windows (radio, own frame): E 2375-2575; blueB 2300-2575; total = blueB + core(2575-2875) = 2300-2875; core 2575-2875.
import numpy as np, pandas as pd, sys
sys.path.insert(0,'/workspace/scripts')
from q3_pos import POS
O='/workspace/out/q13/'
M={o:dict(np.load(O+f'q13a_maps_{o}.npz')) for o in ('Xb4d','Xb4b')}
rng=np.random.default_rng(11); ex=[POS['N'],POS['S'],POS['Bridge'],POS['Soutflow'],(-1.8,-2.4)]
rob=lambda a: 1.4826*np.median(np.abs(a-np.median(a)))
def am(Mp,key,x,y,rad):
    m=np.hypot(Mp['dx']-x,Mp['dy']-y)<=rad; return Mp[key][m].mean()
def nul(rad,n=300):
    out=[]
    while len(out)<n:
        x,y=rng.uniform(-7.5,7.5,2)
        if np.hypot(x,y)>7.5 or min(np.hypot(x-a,y-b) for a,b in ex)<3.0: continue
        if am(M['Xb4d'],'HCO+_core',x,y,rad)>0.25: continue
        out.append((x,y))
    return out
W={'E-window 2375-2575':lambda Mp,l,x,y,r: am(Mp,f'{l}_blue',x,y,r),
   'blueB 2300-2575':lambda Mp,l,x,y,r: am(Mp,f'{l}_blueB',x,y,r),
   'core 2575-2875':lambda Mp,l,x,y,r: am(Mp,f'{l}_core',x,y,r),
   'whole broad 2300-2875':lambda Mp,l,x,y,r: am(Mp,f'{l}_blueB',x,y,r)+am(Mp,f'{l}_core',x,y,r)}
rows=[]
for rad in (0.5,1.0,1.5):
    P=nul(rad); x,y=POS['Soutflow']
    for o in ('Xb4d','Xb4b'):
        for wn,fn in W.items():
            h=fn(M[o],'HCN',x,y,rad); c=fn(M[o],'HCO+',x,y,rad)
            sh=rob(np.array([fn(M[o],'HCN',*p,rad) for p in P])); sc=rob(np.array([fn(M[o],'HCO+',*p,rad) for p in P]))
            r=h/c; er=abs(r)*np.hypot(sh/h,sc/c) if (h>2*sh and c>2*sc) else np.nan
            rows.append(dict(rad=rad,ous=o,window=wn,HCN=h,eHCN=sh,HCO=c,eHCO=sc,ratio=r,e_ratio=er,HCO_snr=c/sc))
T=pd.DataFrame(rows); T.round(3).to_csv(O+'q13d_ratio.csv',index=False)
pd.set_option('display.width',200); print(T.round(3).to_string(index=False))
