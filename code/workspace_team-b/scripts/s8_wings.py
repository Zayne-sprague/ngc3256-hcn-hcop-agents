# High-velocity wing fluxes (blue: -450..-150, red: +150..+450 km/s rel. VSYS) in region beam apertures; empirical box-sum noise from line-free channels
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
import pandas as pd
SETS={('Xb4b','spw25'):('HCO+ 3-2',267.55763),('Xb4b','spw27'):('HCN 3-2',265.88618),('Xb4d','spw25'):('HCO+ 3-2',267.55763),('Xb4d','spw27'):('HCN 3-2',265.88618),
      ('Xb4f','spw16'):('HCO+ 3-2',267.55763),('Xb4f','spw18'):('HCN 3-2',265.88618),('Xb47','spw18'):('HCO+ 4-3',356.73424),('Xb47','spw16'):('HCN 4-3',354.50547),
      ('Xb55','spw18'):('HCO+ 2-1',178.37506),('Xb55','spw16'):('HCN 2-1',177.26111)}
rows=[]
for (ous,spw),(lab,f0) in SETS.items():
    z=np.load(f'/workspace/out/regspec/{ous}.{spw}.cube.pbcor.npz'); v=vrad(z['nu'],f0*1e9)-VSYS; dv=abs(np.median(np.diff(v)))
    for reg in REGIONS:
        if ous in ('Xb4f','Xb47','Xb55') and reg not in ('N','S'): continue
        y=z[reg+'_ap']; ok=np.isfinite(y)
        fl=ok&(abs(v)>500)&(abs(v)<1000)
        if fl.sum()<30: continue
        yb=y[fl]
        for _ in range(4):
            med=np.median(yb); s=1.4826*np.median(np.abs(yb-med)); yb=yb[abs(yb-med)<3.5*s]
        base=np.median(yb); yy=y-base; out=dict(ous=ous,line=lab,region=reg)
        for nm,(a,b) in {'blue':(-450,-150),'red':(150,450),'core':(-150,150)}.items():
            m=ok&(v>=a)&(v<=b); n=m.sum()
            if n<0.8*(b-a)/dv: continue
            sums=[yy[i:i+n].sum() for i in range(0,len(yy)-n) if fl[i:i+n].all()]
            sig=(np.std(sums) if len(sums)>4 else np.std(yb)*np.sqrt(n))*dv
            out[f'W_{nm}']=yy[m].sum()*dv; out[f'sig_{nm}']=sig
        if 'W_blue' in out and 'W_core' in out:
            out['blue/core']=out['W_blue']/out['W_core']; out['red/core']=out['W_red']/out['W_core'] if 'W_red' in out else np.nan
        rows.append(out)
df=pd.DataFrame(rows); df.to_csv('/workspace/out/s8_wings.csv',index=False)
pd.set_option('display.width',250); print(df.round(3).to_string())
