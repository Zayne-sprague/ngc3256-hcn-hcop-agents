# Q16r2 (9): one table of S (and N) profile dchi2(T) from s81 (relative to each variant's own minimum, class-pooled min as in p0012), allowed T range (dchi2<4), dchi2 at 80/100 K.  -> out/q16b/s84.txt
import glob,pandas as pd,numpy as np
D=pd.concat([pd.read_csv(f) for f in sorted(glob.glob('/workspace/out/q16b/s81_*.csv'))]); L=[]
order=['base','aper32','aper43','aper','rec','decB10','decB30','decN10','decN20','decN30','decN30rec','brk']
for g in ('base','tbg5'):
    for reg in ('S','N'):
        for v in order:
            x=D[(D.grid==g)&(D.region==reg)&(D.variant==v)]
            if x.empty: continue
            q=x.groupby('T').chi2.min(); m=q.min(); d=q-m; ok=[int(t) for t in q.index[d<4]]
            L.append(f'Tbg {"2.73" if g=="base" else "5   "} K  {reg} {v:9s} chi2min {m:5.1f} at {int(q.idxmin()):3d} K | dchi2<4: {min(ok):3d}-{max(ok):3d} K | dchi2: 20K {d[20]:5.1f} 30K {d[30]:5.1f} 50K {d[50]:5.1f} 65K {d[65]:5.1f} 80K {d[80]:5.1f} 100K {d[100]:5.1f} 150K {d[150]:5.1f}')
    L.append('')
open('/workspace/out/q16b/s84.txt','w').write('\n'.join(L)); print('\n'.join(L))
