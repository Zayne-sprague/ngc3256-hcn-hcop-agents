# Q15d: baseline systematics of the red wing. FREE channels are all blue of the line (HCN 1380-2100; HCO+ 1280-1680,1880-2100) so the red side is a linear EXTRAPOLATION.
# Re-derive ridge/lobe wing fluxes with baseline order 0/1/2, with alternative free ranges, and with red-edge channels added to the baseline (assumes wing ends before +520), Xb4d, r=1" mean.
import sys; sys.path.insert(0,'/workspace/scripts')
from q13a_xb4_compare import *
from q3_pos import POS
O='/workspace/out/q15/'; log=open(O+'q15d.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
Sx,Sy=POS['S']; RIDGE=(Sx+0.0,Sy+3.3); LOBE=POS['Soutflow']; N=POS['N']
V0=2730.
for sp in ('HCN','HCO+'):
    c=vcube(sp,'Xb4d',1.05,9.); v=c['v']; dx,dy,dv=c['dx'],c['dy'],c['dv']; a=c['a']
    pr(f'## {sp}: cube v {v.min():.0f}..{v.max():.0f} km/s (radio, own rest freq); channel {dv:.2f} km/s')
    fr0=FREE[sp]
    variants={'free-blue only, order1 (baseline used in Q13/Q15b)':(fr0,1),'free-blue, order0':(fr0,0),'free-blue, order2':(fr0,2),
              'blue + red edge >3290 (order1)':(fr0+[(3290,v.max())],1),'blue + red edge >3290 (order2)':(fr0+[(3290,v.max())],2),
              'blue + 3200-edge (order1)':(fr0+[(3200,v.max())],1)}
    for nm,(free,od) in variants.items():
        m=np.zeros(len(v),bool)
        for p,q in free: m|=(v>=p)&(v<=q)
        X=np.vander(v-v[m].mean(),od+1); A=a[m].reshape(m.sum(),-1); coef=np.linalg.lstsq(X[m],A,rcond=None)[0]
        r=a-(X@coef).reshape(len(v),*a.shape[1:])
        out=[]
        for pn,(x,y) in (('ridge',RIDGE),('lobe',LOBE),('N',N),('ctrl(+3.5,-1.5)',(3.5,-1.5)),('ctrl(-3,-6)',(-3,-6))):
            ap=np.hypot(dx-x,dy-y)<=1.0; s=r[:,ap].mean(1)
            for wn,(p,q) in (('red230-470',(2960,3200)),('blue-mirror',(2260,2500))):
                out.append(f'{pn}:{wn}={s[(v>=p)&(v<q)].sum()*dv:+.3f}')
        pr(f'{nm:55s} '+' '.join(out))
