# Re-assign every blind detection (out/blind_detections_all.csv) with the bug-free matcher (ident.py): local line velocity of the region (HCO+ Gaussian centre of the same OUS),
# tolerance +-max(80, 0.5*FWHM_region) km/s, all species with Eu<200 K, no Aij cut, ranked by LTE strength at 30 K; plus the same at the systemic velocity 2780.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from ident import ident,fmt
from ch3oh_lines import ref
b=pd.read_csv('/workspace/out/blind_detections_all.csv').drop(columns=['candidates'])
rows=[]
for r in b.itertuples():
    v0,fw=ref(r.ous,r.region); tol=max(80.,0.5*fw)
    s=ident(r.fsky,v0,tol,emax=200,top=4); s2=ident(r.fsky,2780.,tol,emax=200,top=4)
    rows.append(dict(v_local=v0,tol=tol,id_local=fmt(s[s.cur]) if s.cur.any() else 'none among curated; '+fmt(s.head(1)),id_vsys=fmt(s2[s2.cur]) if s2.cur.any() else 'none among curated'))
o=pd.concat([b,pd.DataFrame(rows)],axis=1)
# short assignment: top-ranked curated species name+freq at local velocity
def short(t):
    return t.split(';')[0][:70]
o['assign']=o.id_local.map(short)
o.to_csv('/workspace/out/blind_detections_reassigned.csv',index=False)
pd.set_option('display.width',300,'display.max_rows',500,'display.max_colwidth',75)
print(o[['ous','region','fsky','width','snr','v_local','assign']].round(3).sort_values(['ous','fsky']).to_string())
