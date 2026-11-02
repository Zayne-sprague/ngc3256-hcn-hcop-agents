# Q1: HCN/HCO+ 3-2 (Tb) in blue/core/red velocity windows (S-local centre 2722 km/s, +-150 / 300-km/s windows as in p0001 wing_null.py) at a matched 1.05" beam, Xb4d.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import pandas as pd
W={'blue(2272-2572)':(2272,2572),'core(2572-2872)':(2572,2872),'red(2872-3172)':(2872,3172)}
POS={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15),'Wclump':(-1.8,-2.4),'Soutflow':(-0.15,-4.5)}
rows=[]
for wn,win in W.items():
    M={}
    for sp in ('HCN','HCO+'):
        k=f'{sp} 3-2@Xb4d'; d=load(k,window=win); sm=smooth(d,1.05); pr=noise_profile(d,1.05); M[sp]=(d,sm,pr)
    cT=Tb_factor(LINES['HCN 3-2@Xb4d'][2]/(1+VSYS/C),1.05)/Tb_factor(LINES['HCO+ 3-2@Xb4d'][2]/(1+VSYS/C),1.05)
    for rn,(x,y) in POS.items():
        r=np.hypot(x,y); v={}
        for sp in M: d,sm,pr=M[sp]; v[sp]=(sample(sm,d,x,y),sigma_at(pr,r))
        (a,ea),(b,eb)=v['HCN'],v['HCO+']
        row=dict(window=wn,region=rn,I_HCN=a,e_HCN=ea,I_HCO=b,e_HCO=eb,snr_HCN=a/ea,snr_HCO=b/eb)
        if a/ea>=3 and b/eb>=3: R=a/b*cT; row.update(R=R,eR=R*np.hypot(ea/a,eb/b),kind='value')
        elif a/ea>=3: row.update(R=a/(3*eb)*cT,kind='lower limit')
        elif b/eb>=3: row.update(R=3*ea/b*cT,kind='upper limit')
        rows.append(row)
T=pd.DataFrame(rows); T.to_csv('/workspace/out/q1/q1_wing_ratios.csv',index=False)
pd.set_option('display.width',250); print(T.round(3).to_string())
