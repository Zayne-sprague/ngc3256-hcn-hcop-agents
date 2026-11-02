# Q13a2: Xb4b(1.05"-smoothed)/Xb4d recovery per line, window and aperture, and the likelihood of the recovery fraction R that makes the
# Xb4b excess E consistent with the Xb4d excess.  Reads out/q13/q13a_maps_*.npz (run q13a_xb4_compare.py first).
# Recovery R = S(Xb4b)/S(Xb4d) of the same window sum in the same aperture, both at 1.05"; errors from random-aperture scatter in each array.
import numpy as np, pandas as pd, sys
sys.path.insert(0,'/workspace/scripts')
from q3_pos import POS
O='/workspace/out/q13/'
M={o:dict(np.load(O+f'q13a_maps_{o}.npz')) for o in ('Xb4d','Xb4b')}
rng=np.random.default_rng(7)
ex=[POS['N'],POS['S'],POS['Bridge'],POS['Soutflow'],(-1.8,-2.4)]
def am(Mp,key,x,y,rad):
    m=np.hypot(Mp['dx']-x,Mp['dy']-y)<=rad; return Mp[key][m].mean()
def nullpos(rad,n=400):
    out=[]
    while len(out)<n:
        x,y=rng.uniform(-7.5,7.5,2)
        if np.hypot(x,y)>7.5 or min(np.hypot(x-a,y-b) for a,b in ex)<3.0: continue
        if am(M['Xb4d'],'HCO+_core',x,y,rad)>0.25: continue
        out.append((x,y))
    return out
rob=lambda a: 1.4826*np.median(np.abs(a-np.median(a)))
rows=[]
for rad in (1.0,1.5,2.0):
    P=nullpos(rad)
    for pn in ('Soutflow','S','N','Bridge'):
        x,y=POS[pn]
        for line in ('HCN','HCO+'):
            for w in ('blue','blueB','core','red'):
                key=f'{line}_{w}'
                d=am(M['Xb4d'],key,x,y,rad); b=am(M['Xb4b'],key,x,y,rad)
                sd=rob(np.array([am(M['Xb4d'],key,*p,rad) for p in P])); sb=rob(np.array([am(M['Xb4b'],key,*p,rad) for p in P]))
                R=b/d if abs(d)>3*sd else np.nan
                rows.append(dict(rad=rad,pos=pn,line=line,win=w,Xb4d=d,e_d=sd,Xb4b=b,e_b=sb,R=R,e_R=abs(R)*np.hypot(sd/d,sb/b) if R==R else np.nan))
T=pd.DataFrame(rows); T.round(4).to_csv(O+'q13a2_recovery.csv',index=False)
pd.set_option('display.width',200)
print(T[(T.pos=='Soutflow')].round(3).to_string(index=False))
print(T[(T.win=='core')&(T.rad==1.0)].round(3).to_string(index=False))
# excess E in each array at several apertures with k fixed 0.574; recovery likelihood of R
k=0.574; out=[]
for rad in (0.5,1.0,1.5,2.0):
    P=nullpos(rad,300); x,y=POS['Soutflow']
    for wn in ('blue','blueB'):
        def E(o,p): return am(M[o],f'HCN_{wn}',*p,rad)-k*am(M[o],f'HCO+_{wn}',*p,rad)
        Ed=E('Xb4d',(x,y)); Eb=E('Xb4b',(x,y))
        sd=rob(np.array([E('Xb4d',p) for p in P])); sb=rob(np.array([E('Xb4b',p) for p in P]))
        # R with Xb4b E = R*Xb4d E : weighted LS using Xb4b only, with Xb4d E as template (its error propagates)
        R=Eb/Ed; eR=np.hypot(sb/Ed,Eb*sd/Ed**2)
        out.append(dict(rad=rad,win=wn,E_Xb4d=Ed,sig_d=sd,snr_d=Ed/sd,E_Xb4b=Eb,sig_b=sb,snr_b=Eb/sb,expect_snr_b_if_R1=Ed/sb,R=R,e_R=eR,tension_R1_sigma=(1-R)/eR,diff_sigma=(Ed-Eb)/np.hypot(sd,sb)))
U=pd.DataFrame(out); U.round(4).to_csv(O+'q13a2_E_apertures.csv',index=False); print(U.round(3).to_string(index=False))
