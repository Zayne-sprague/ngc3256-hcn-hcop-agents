# Q25: uv-recovery test. Apply a cosine-tapered high-pass (u1->u2 kl, angular) to the best-recovered maps (Michiyama 1-0, Xb55 2-1) and see how annulus HCN/HCO+ and the 8-12/4-8 decline respond;
# check whether the Harada 1-0 result (cumulative 12" flux 0.75-0.78 of Michiyama) is reproduced by a high-pass of Michiyama.
import sys; sys.path.insert(0,'/workspace/scripts')
from q25_ratio10 import *
from q11_likeforlike_lib import bp
M={}
for nm in ('Harada','Michiyama'):
    M[nm]=[ (lambda d:(d,sm(d)))(load10(P[nm],f0)) for f0 in REST.values()]
X=[(lambda d:(d,smooth(d,10.)))(load(k,nullspec=False)) for k in ('HCN 2-1','HCO+ 2-1')]
def run(name,pair,u1,u2):
    f=[]; 
    for d,s in pair:
        sb=bp(s,d['pix'],u1,u2) if u2>0 else s; f.append(ann(d,sb)[0])
    h,o=f; R=h/o; return R,h,o
for name,pair in (('Michiyama 1-0',M['Michiyama']),('Harada 1-0',M['Harada']),('Xb55 2-1',X)):
    print('==',name)
    for u1,u2 in ((0,0),(0.5,1.5),(1,2),(1.5,2.5),(2,3),(3,4)):
        R,h,o=run(name,pair,u1,u2)
        k=ann(pair[0][0],pair[0][1])[0]
        print(f' hp {u1}-{u2} kl: HCN {h.round(2)} HCO+ {o.round(2)}  R {R.round(3)}  8-12/4-8 {R[3]/R[1]:.3f}  flux kept(0-12 HCN) {h[[0,1,3]].sum()/k[[0,1,3]].sum():.2f}')
