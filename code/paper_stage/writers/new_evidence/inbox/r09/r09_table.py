# Reads out/q66s/uvmatched_32_43.csv (scripts/s438_q66_32_vs_43.py) ; python3 -I r09_table.py
import pandas as pd,numpy as np
T=pd.read_csv('/workspace/out/q66s/uvmatched_32_43.csv').set_index(['arr','var'])
pairs=[('b5','b5','5.0" all-uv(4-3) / >8 kl(3-2)'),('b6p5','b6p5','6.5" native/native'),('b5u15','b5u11','5.0" >15 kl / >11.3 kl'),('b5u20','b5u15','5.0" >20 kl / >15 kl')]
sg=lambda r:(r.NS_hi-r.NS_lo)/2
for a,b,l in pairs:
    x=T.loc[('X4f32',b)]; y=T.loc[('X47',a)]
    R=y.NS/x.NS; e=R*np.hypot(sg(x)/x.NS,sg(y)/y.NS)
    print(l,'| 3-2 NS %.2f +-%.2f | 4-3 NS %.2f +-%.2f | HCN N/S 3-2 %.2f 4-3 %.2f | HCO+ N/S 3-2 %.2f 4-3 %.2f | ratio %.2f +-%.2f | deb A/B/C ratio %s'%(
      x.NS,sg(x),y.NS,sg(y),x.I_HCN_N/x.I_HCN_S,y.I_HCN_N/y.I_HCN_S,x.I_HCOp_N/x.I_HCOp_S,y.I_HCOp_N/y.I_HCOp_S,R,e,
      '/'.join('%.2f(%.2f->%.2f)'%(y['deb_'+m]/x['deb_'+m],x.NS,x['deb_'+m]) for m in 'ABC')))
