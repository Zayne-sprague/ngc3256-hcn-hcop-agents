# Q15g: HCN/HCO+ 4-3 (Xb47, native 4.1x2.7") vs 3-2 (Xb4d) red-wing check at the ridge, both smoothed to a 4.5" circular beam (Q13c machinery).
# The ridge is only 1.7" from N and 3.3" from S: at 4.5" it is fully blended with the N/S wings, so the test is the Tb ratio r43/32 of the red wing (+230..+470) in an r=1" aperture at the ridge vs
# the same ratio of the blue wing at the lobe and at N, and the null scatter of random apertures. Baseline: HCN43 free 1790-2100 only (order 0 and 1); HCO43 free 2147-2300 and 3450-3660 (order 1).
import sys; sys.path.insert(0,'/workspace/scripts')
from q13c_hcn43 import *
O15='/workspace/out/q15/'; log=open(O15+'q15g.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
Sx,Sy=POS['S']; RIDGE=(0.3,Sy+3.3); LOBE=POS['Soutflow']; N=POS['N']
FRE={'HCN32':L['HCN32'][3],'HCO32':L['HCO32'][3],'HCN43':[(1790,2100)],'HCO43':[(2147,2300),(3450,3660)]}
WR=(2960,3200); WB=(2260,2500)
Dd={n:load(n) for n in L}
nu32=265.886434/1.0093; nu43=354.505477/1.0093
rng=np.random.default_rng(7); ex=[N,POS['S'],RIDGE,LOBE,POS['Bridge']]
P=[]
while len(P)<200:
    x,y=rng.uniform(-9,9,2)
    if np.hypot(x,y)<=9 and min(np.hypot(x-a,y-b) for a,b in ex)>=4.5: P.append((x,y))
for od43 in (0,1):
  pr(f'### baseline order for 4-3 HCN: {od43} (3-2 order 1)')
  W={}
  for n,c in Dd.items():
    v=c['v']; free=FRE[n]; od=od43 if n=='HCN43' else 1
    m=np.zeros(len(v),bool)
    for a,b in free: m|=(v>=a)&(v<=b)
    X=np.vander(v-v[m].mean(),od+1); A=c['a'][m].reshape(m.sum(),-1); coef=np.linalg.lstsq(X[m],A,rcond=None)[0]
    r=c['a']-(X@coef).reshape(len(v),*c['a'].shape[1:])
    for w,(a,b) in (('red',WR),('blue',WB)):
        mm=(v>=a)&(v<b); W[(n,w)]=r[mm].sum(0)*c['dv']
  def ap(n,w,x,y): c=Dd[n]; return W[(n,w)][np.hypot(c['dx']-x,c['dy']-y)<=1.0].mean()
  for line,(h32,h43) in (('HCN',('HCN32','HCN43')),('HCO+',('HCO32','HCO43'))):
    for w,pn,(x,y) in (('red','ridge',RIDGE),('blue','lobe',LOBE),('red','N',N),('red','lobe(ctrl)',LOBE)):
        s32=ap(h32,w,x,y); s43=ap(h43,w,x,y)
        n32=np.array([ap(h32,w,*q) for q in P]); n43=np.array([ap(h43,w,*q) for q in P])
        e32=n32.std(); e43=n43.std()
        r=s43/s32*(nu32/nu43)**2 if s32>0 else np.nan
        pr(f'{line:5s} {pn:12s} {w:4s}: S32={s32:+.3f}+-{e32:.3f} (null mean {n32.mean():+.3f})  S43={s43:+.3f}+-{e43:.3f} (null mean {n43.mean():+.3f}) Jy/b km/s; z32={s32/e32:.1f} z43={s43/e43:.1f}; Tb r43/32 = {r:.2f} +- {abs(r)*np.hypot(e32/max(s32,1e-9),e43/max(abs(s43),1e-9)):.2f}')
