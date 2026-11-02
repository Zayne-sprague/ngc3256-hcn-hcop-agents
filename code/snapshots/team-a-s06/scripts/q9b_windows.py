# Q9b fix 1+5: HCN/HCO+ 3-2 ratio of velocity-window-integrated fluxes (no S/N clipping), blue/red separately, N/S/Bridge,
# Xb4d@1.05" and Xb4b@0.30"; errors = null (random positions x line-free sliding windows) + linear-baseline-fit variance; censored limits.
# Also aperture-centre shift test (HCN vs HCO+ peak offset, +-0.1"). python3 -I q9b_windows.py
import sys; sys.path.insert(0,'/workspace/scripts')
from q9b_lib import *
import pandas as pd
rng=np.random.default_rng(99)
BINS=[(0,60),(60,100),(100,150),(150,250),(100,260)]
WINS=[(a,b) for a,b in BINS]+[(-b,-a) for a,b in BINS]
CASES=[('Xb4d',(1.05,1.05,0.),{'N':[(0,.5),(.5,1),(1,2)],'S':[(0,.5),(.5,1),(1,2)],'Bridge':[(0,1)]}),
       ('Xb4b',(0.30,0.30,0.),{'N':[(0,.25),(.25,.5),(.5,1)],'S':[(0,.25),(.25,.5),(.5,1)]}),
       ('Xb4b',(1.05,1.05,0.),{'N':[(0,.5),(.5,1),(1,2)],'S':[(0,.5),(.5,1),(1,2)]})]
rows=[]
for ous,tb,regs in CASES:
    for reg,ann in regs.items():
        for r1,r2 in ann:
            (x0,y0),res=window_table(ous,tb,reg,r1,r2,WINS+[(-60,60)],rng)
            # core |v|<60 as one window, plus stacked blue+red
            def get(lo,hi): a=res[('HCN',lo,hi)]; b=res[('HCO+',lo,hi)]; return a,b
            for (lo,hi) in WINS+[(-60,60)]:
                a,b=get(lo,hi); R,eR,k,L=ratio(a[0],a[1],b[0],b[1])
                rows.append(dict(ous=ous,beam=tb[0],reg=reg,r1=r1,r2=r2,x0=round(x0,2),y0=round(y0,2),vlo=lo,vhi=hi,HCN=a[0],eHCN=a[1],HCO=b[0],eHCO=b[1],R=R,eR=eR,kind=k,limit=L,snrH=a[0]/a[1],snrO=b[0]/b[1]))
            for lo,hi in BINS:   # stacked red+blue wing (errors added in quadrature, null windows independent)
                a1,b1=get(lo,hi); a2,b2=get(-hi,-lo)
                A,eA=a1[0]+a2[0],np.hypot(a1[1],a2[1]); B,eB=b1[0]+b2[0],np.hypot(b1[1],b2[1]); R,eR,k,L=ratio(A,eA,B,eB)
                rows.append(dict(ous=ous,beam=tb[0],reg=reg,r1=r1,r2=r2,x0=round(x0,2),y0=round(y0,2),vlo=lo,vhi=hi,HCN=A,eHCN=eA,HCO=B,eHCO=eB,R=R,eR=eR,kind=k,limit=L,snrH=A/eA,snrO=B/eB,stack=True))
    print(ous,tb[0],'done',flush=True)
df=pd.DataFrame(rows); df['stack']=df['stack'].fillna(False); df.to_csv('/workspace/out/q9b/windows.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',900)
d=df.copy(); d['Tb']=d.R*TBF
for (ous,beam,reg,r1,r2),g in d.groupby(['ous','beam','reg','r1','r2'],sort=False):
    print(f'\n{ous}@{beam}" {reg} r={r1}-{r2}"  centre ({g.x0.iloc[0]},{g.y0.iloc[0]})')
    print(g.sort_values(['stack','vlo'])[['vlo','vhi','stack','HCN','eHCN','HCO','eHCO','R','eR','kind','limit']].round(4).to_string(index=False))
# aperture-centre shift test: Xb4b@0.30 N and S r<0.25, shift +-0.1" in x and y; window ratio core(+-60), blue 60-150, red 60-150
print('\n== centre-shift test (Xb4b@0.30, r<0.25") ==')
sh=[]
for reg in ('N','S'):
    x0,y0=centre('Xb4b',(0.3,0.3,0.),reg); hx,hy=centre('Xb4b',(0.3,0.3,0.),reg,'HCN')
    print(reg,'HCO+ peak',round(x0,3),round(y0,3),' HCN peak',round(hx,3),round(hy,3),'offset (arcsec)',round(hx-x0,3),round(hy-y0,3))
    for dx,dy in ((0,0),(.1,0),(-.1,0),(0,.1),(0,-.1),(hx-x0,hy-y0)):
        v={}
        for ln in ('HCN','HCO+'):
            c=cube('Xb4b',ln); s=smcube('Xb4b',ln,(.3,.3,0.)); m=apmask(c,x0+dx,y0+dy,0,.25); sp=apspec(c,s,m,(.3,.3,0.))
            for nm,(lo,hi) in {'core':(-60,60),'blue':(-150,-60),'red':(60,150)}.items(): v[(ln,nm)]=winflux(c,sp,lo,hi)
        sh.append(dict(reg=reg,dx=round(dx,3),dy=round(dy,3),**{nm:v[('HCN',nm)]/v[('HCO+',nm)] for nm in ('core','blue','red')}))
sh=pd.DataFrame(sh); print(sh.round(3).to_string(index=False)); sh.to_csv('/workspace/out/q9b/centre_shift.csv',index=False)
