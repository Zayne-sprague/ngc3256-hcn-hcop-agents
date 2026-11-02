"""Q54 r2 fix1b: HCN and HNC 1-0 integrals vs circular beam for 2015.1.00993 (X961) and 2015.1.00412 (X145/X149): does the 00412 deficit (HNC) grow with aperture (missing extended flux)?
Output out/q54/beamtrend.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q54_chem as c
rows=[]
for key in ('HCN1-0_00993','HCN1-0_00412','HNC1-0_00993','HNC1-0_00412'):
    for bm in (2.2,2.5,3.5,5.0,7.0):
        for reg in 'NS':
            try: I,e,nb,sg,_,_=c.meas(key,bm,reg)
            except Exception as ex: print(key,bm,reg,ex); continue
            rows.append(dict(line=key,beam=bm,reg=reg,I=I,e=e)); print(key,bm,reg,round(I,1),round(e,1),flush=True)
d=pd.DataFrame(rows).round(2); d.to_csv('/workspace/out/q54/beamtrend.csv',index=False)
p=d.pivot_table(index=['reg','beam'],columns='line',values='I')
p['HCN 412/993']=p['HCN1-0_00412']/p['HCN1-0_00993']; p['HNC 412/993']=p['HNC1-0_00412']/p['HNC1-0_00993']
p['HNC/HCN 993']=p['HNC1-0_00993']/p['HCN1-0_00993']; p['HNC/HCN 412']=p['HNC1-0_00412']/p['HCN1-0_00412']
print(p.round(3).to_string())
