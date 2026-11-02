# Q16 (4): tabulate the profile chi2(T) deltas from s73 outputs (dchi2 relative to each variant's own minimum), ranges at dchi2<4 (grid-limited), and P(T>=80) exclusion level. -> out/q16/s74.txt
import glob,pandas as pd,numpy as np
D=pd.concat([pd.read_csv(f) for f in glob.glob('/workspace/out/q16/s73_*.csv')])
L=[]
for (g,reg,var),x in D.groupby(['grid','region','variant']):
    q=x.groupby('T').chi2.min(); m=q.min(); d=q-m; ok=[int(t) for t in q.index[d<4]]
    L.append(f'{g:6s} {reg} {var:12s} chi2min {m:5.1f} at T={int(q.idxmin()):3d} K | dchi2<4: {min(ok)}-{max(ok)} K ({len(ok)} grid pts) | dchi2(80)={d[80]:.1f} dchi2(100)={d[100]:.1f} dchi2(150)={d[150]:.1f} dchi2(20)={d[20]:.1f} dchi2(30)={d[30]:.1f}')
open('/workspace/out/q16/s74.txt','w').write('\n'.join(sorted(set(L)))); print('\n'.join(sorted(set(L))))
