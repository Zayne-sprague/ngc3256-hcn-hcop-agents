"""Q54 r3: CN 1-0 (2018.1.00223 X4b 0.49", X4d, X4f compact) vs HCN 1-0 (00965 0.76x0.63", 00993) at N,S in 0.9-2.5" circular beams; missing flux CN X4b vs X4d vs X4f.
Output out/q54/cn_subarcsec.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
src=open('/workspace/scripts/q54b_cn_ratios.py').read().split('PAIRS=')[0]
exec(src)
A2='/archive/files/2018.1.00223.S/member.uid___A001_'
c.L['CN1-0_X4b']=(A2+'X133f_X4b.NGC3256_sci.spw27.cube.I.pbcor.fits',113.490982)
c.L['CN1-0_X4f']=(A2+'X133f_X4f.NGC3256_sci.spw18.cube.I.pbcor.fits',113.490982)
for k in ['X4b','X4f']:
    W['CN1-0_'+k]=dict(W['CN1-0_X4d'])
P=[('CN1-0_X4b','CN1-0_X4d',[0.9,1.4,2.5]),('CN1-0_X4b','HCN1-0_00965',[0.9,1.0,1.4,2.5]),('CN1-0_X4d','HCN1-0_00965',[0.9,1.0,1.4,2.5]),
   ('CN1-0_X4b','HCN1-0_00993',[1.4,2.5]),('CN1-0_X4d','HCN1-0_00993',[1.4,2.5]),('HCN1-0_00965','HCN1-0_00993',[0.9,1.4,2.5]),
   ('CN1-0_X4f','CN1-0_X4d',[17.]),('CN1-0_X4f','HCN1-0_00993',[17.]),('HCN1-0_00993','HCN1-0_00965',[17.])]
rows=[]
for a,b,bms in P:
    for bm in bms:
        for reg in 'NS':
            try: Aa=get(a,bm,reg); B=get(b,bm,reg)
            except Exception as ex: print('skip',a,b,bm,reg,type(ex).__name__,ex); continue
            r,er=ratio(Aa,B); rows.append(dict(num=a,den=b,beam=bm,reg=reg,Inum=Aa[0],enum=Aa[1],Iden=B[0],eden=B[1],ratio=r,err=er)); print(a,b,bm,reg,round(r,3),round(er,3),flush=True)
pd.DataFrame(rows).round(4).to_csv('/workspace/out/q54/cn_subarcsec.csv',index=False)
