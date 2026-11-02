# Q11: effective primary-beam response vs radius from the pbcor null-window noise (rms ∝ 1/PB) for each map; gives the PB>0.5 radius.
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
for k in ('HCN 2-1','HCO+ 2-1','HCN 3-2@Xb4f','HCO+ 3-2@Xb4f','HCN 3-2@Xb4d','HCO+ 3-2@Xb4d'):
    d=load(k,nullspec=True); r=np.hypot(d['dx'],d['dy']); N=np.array([n for _,n in d['nulls']])
    rms=[]; edges=[0,4,8,12,16,20,26,32,40,48,60]
    for a,b in zip(edges[:-1],edges[1:]):
        m=(r>=a)&(r<b)
        rms.append(np.sqrt(np.mean(N[:,m]**2)) if m.sum()>3 else np.nan)
    rms=np.array(rms); pb=rms[0]/rms
    print(k,'beam',np.round(d['beam'],1),'nnull',len(N),'rmax',round(r.max()))
    print('   r  ',edges[1:]); print('   PB~',np.round(pb,2))
