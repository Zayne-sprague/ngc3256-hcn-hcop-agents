"""Q54 r2 fix4: ratios at the finest common beams (0.8-1.0") for the 3-2 pairs; HCN 3-2 from Xb4b (0.25") tests missing flux vs Xb4d.
Output out/q54/subarcsec.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q4_lib, q54_chem as c, q54b_cn_ratios as q
src=open('/workspace/scripts/q54b_cn_ratios.py').read().split('PAIRS=')[0]
exec(src)   # defines cn(), get(), ratio(), W, c
c.L['HCN3-2_Xb4b']=('/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw27.cube.I.pbcor.fits',265.886434)
P=[('CN3-2_X1df1','HCN3-2_Xb4d',[0.95,1.1]),('CN3-2_X1df1','HCN3-2_Xb4b',[0.95,1.1]),
   ('HNC3-2_00412','HCN3-2_Xb4b',[0.8,0.95,1.1]),('HNC3-2_00412','HCN3-2_Xb4d',[0.95,1.1]),
   ('HCN3-2_Xb4b','HCN3-2_Xb4d',[0.95,1.1])]
rows=[]
for a,b,bms in P:
    for bm in bms:
        for reg in 'NS':
            try: A=get(a,bm,reg); B=get(b,bm,reg)
            except Exception as ex: print('skip',a,b,bm,reg,ex); continue
            r,er=ratio(A,B); rows.append(dict(num=a,den=b,beam=bm,reg=reg,Inum=A[0],enum=A[1],Iden=B[0],eden=B[1],ratio=r,err=er)); print(a,b,bm,reg,round(r,3),round(er,3),flush=True)
pd.DataFrame(rows).round(4).to_csv('/workspace/out/q54/subarcsec.csv',index=False)
