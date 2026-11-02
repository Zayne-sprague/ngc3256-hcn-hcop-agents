"""Q54: CN (3-2, 1-0 hyperfine complexes, frequency windows with a linear baseline) vs HCN/HNC at matched circular beams, Tb (K km/s).
CN integrals are over the full complex(es) (both J groups); HCN/HNC use the fixed velocity window of q54_chem.py. Output out/q54/ratios.csv"""
import sys,os,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q4_lib, q54_chem as c
from q4_lib import spectrum,kb,vel
C=299792.458
W={'CN3-2_X1df1':dict(win=(336.68,337.30),bl=[(336.35,336.55),(337.45,337.75)],rest=340.2477),
   'CN1-0_X4d':dict(win=(111.98,112.58),bl=[(111.60,111.90),(112.70,113.30)],rest=113.490982),
   'CN1-0_main':dict(win=(112.34,112.58),bl=[(111.60,111.90),(112.70,113.30)],rest=113.490982,key='CN1-0_X4d')}
SH={'N':0.,'S':+1.0}   # S line at lower velocity: shift windows by (2785-2750)/c*f
def cn(key,bm,reg):
    w=W[key]; path,_=c.L[w.get('key',key)]; q4_lib.f=lambda *a,**k: path
    fr,s,b,nb=spectrum('x','x',reg,tb=(bm,bm,0.)); y=s*kb(fr.mean(),bm,bm)
    sh=SH[reg]*35./C*fr.mean(); df=abs(np.median(np.diff(fr))); dv=df/fr.mean()*C
    m=np.zeros(len(fr),bool)
    for lo,hi in w['bl']: m|=(fr>lo)&(fr<hi)
    cc=np.polyfit(fr[m]-fr.mean(),y[m],1); r=y[m]-np.polyval(cc,fr[m]-fr.mean()); sg=1.4826*np.median(np.abs(r-np.median(r)))
    wm=(fr>w['win'][0]+sh)&(fr<w['win'][1]+sh); n=wm.sum()
    I=np.sum(y[wm]-np.polyval(cc,fr[wm]-fr.mean()))*dv
    e=np.hypot(sg*dv*np.sqrt(n),sg*dv*n/np.sqrt(m.sum()))
    return I,e,sg
rows=[]
HCN={('1-0',2.5):['HCN1-0_00993','HCN1-0_00412'],('3-2',1.0):None}
def get(key,bm,reg):
    if key in W: return cn(key,bm,reg)
    I,e,nb,sg,_,_=c.meas(key,bm,reg); return I,e,sg
def ratio(a,b):
    r=a[0]/b[0]; return r,abs(r)*np.hypot(a[1]/a[0],b[1]/b[0])
PAIRS=[('CN1-0_X4d','HCN1-0_00993',[2.5]),('CN1-0_X4d','HCN1-0_00412',[2.5]),('CN1-0_main','HCN1-0_00993',[2.5]),
       ('HNC1-0_00993','HCN1-0_00993',[2.5]),('HNC1-0_00412','HCN1-0_00412',[1.4,1.9,2.5]),
       ('CN3-2_X1df1','HCN4-3_X1df1',[1.0,1.4,1.9,2.5]),('HNC3-2_00412','HCN3-2_Xb4d',[1.1,1.4,1.9,2.5]),
       ('CN1-0_X4d','HNC1-0_00993',[2.5]),
       ('CN3-2_X1df1','HCN3-2_Xb4d',[1.1,1.4,1.9,2.5]),('CN3-2_X1df1','CN1-0_X4d',[2.5]),('HCN3-2_Xb4d','HCN1-0_00993',[2.5]),
       ('HCN3-2_Xb4d','HCN1-0_00412',[2.5]),('HCN4-3_X1df1','HCN3-2_Xb4d',[1.4,2.5]),('HNC3-2_00412','HNC1-0_00412',[2.5]),
       ('HNC3-2_00412','HNC1-0_00993',[2.5])]
for a,b,bms in PAIRS:
    for bm in bms:
        for reg in 'NS':
            try: A=get(a,bm,reg); B=get(b,bm,reg)
            except (ValueError,KeyError,FileNotFoundError,OSError) as ex: print('skip',a,b,bm,reg,type(ex).__name__,ex); continue
            r,er=ratio(A,B); rows.append(dict(num=a,den=b,beam=bm,reg=reg,Inum=A[0],enum=A[1],Iden=B[0],eden=B[1],ratio=r,err=er)); print(a,b,bm,reg,round(r,3),round(er,3),flush=True)
pd.DataFrame(rows).round(4).to_csv('/workspace/out/q54/ratios.csv',index=False)
