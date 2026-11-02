# Q35: HCN/CS and HCO+/CS at matched 1.7" (Harada+2018 Table 13, peak pixel, +-150 km/s; parsed in lit/harada_table13_parsed.csv) and N/S, C/S, OS/S, SW/S double ratios; 
# plus H13CN/H13CO+ (thin-isotopologue) from our Q29 numbers. Errors: cube rms of each entry + 10% flux scale per line (uncorrelated within a transition pair taken from the same band: 5% ratio) -> 7% added per ratio.
# usage: python3 -I s201_q35_harada_cs.py   -> out/q35/s201.txt
import pandas as pd, numpy as np, re
d=pd.read_csv('/workspace/lit/harada_table13_parsed.csv'); d['sp']=d.species_transition.str.replace(' ','')
def g(p,s): r=d[(d.pos==p)&(d.sp.str.startswith(s))]; return (r.I_Kkms.iloc[0],r.eI_Kkms.iloc[0]) if len(r) else None
T=[]
def P(s): print(s); T.append(s)
pairs={'1-0/2-1':('HCN1','HCO+1','CS2'),'3-2/5-4':('HCN3','HCO+3','CS5')}
res={}
for lab,(a,b,c) in pairs.items():
    P(f'== HCN/CS, HCO+/CS ({lab}), K km/s ratios at 1.7"')
    for p in ('N','S','C','OS','SW'):
        A,B,Cc=g(p,a),g(p,b),g(p,c)
        if None in (A,B,Cc): continue
        out={}
        for nm,X in (('HCN/CS',A),('HCO+/CS',B)):
            r=X[0]/Cc[0]; e=r*np.sqrt((X[1]/X[0])**2+(Cc[1]/Cc[0])**2+0.07**2); out[nm]=(r,e)
        r=A[0]/B[0]; out['HCN/HCO+']=(r,r*np.sqrt((A[1]/A[0])**2+(B[1]/B[0])**2+0.05**2)); res[(lab,p)]=out
        P(f'{p}: '+'; '.join(f'{k} {v[0]:.2f}+-{v[1]:.2f}' for k,v in out.items()))
    for p in ('N','C','OS','SW'):
        if (lab,p) in res:
            s=[]
            for k in ('HCN/CS','HCO+/CS','HCN/HCO+'):
                x=res[(lab,p)][k][0]/res[(lab,'S')][k][0]; e=x*np.hypot(res[(lab,p)][k][1]/res[(lab,p)][k][0],res[(lab,'S')][k][1]/res[(lab,'S')][k][0]); s.append(f'{k} {p}/S = {x:.2f}+-{e:.2f}')
            P('  '+'; '.join(s))
P('\n== Thin isotopologues (Q29, 2.2" aperture, K km/s; out/q29/s185.txt): H13CN/H13CO+ intensity ratio')
for p,(a,ea,b,eb) in {'N':(1.88,.31,2.50,.34),'S':(0.31,.32,1.66,.36),'Bridge':(1.80,.34,1.92,.35)}.items():
    r=a/b; e=np.hypot(ea/a,eb/b)*r if a>0 else float('nan'); P(f'{p}: {r:.2f}+-{e:.2f}  (3 sigma upper at S: {(a+3*ea)/b:.2f})')
P('N/S (3 sigma lower limit using S upper): >%.1f'%((1.88/2.50)/((0.31+3*0.32)/1.66)))
open('/workspace/out/q35/s201.txt','w').write('\n'.join(T))
