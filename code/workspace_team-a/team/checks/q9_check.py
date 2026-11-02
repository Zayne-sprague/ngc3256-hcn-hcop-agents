import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from q7_lib import *
H=Cube('Xb4d','HCN'); O=Cube('Xb4d','HCO+')
for nm,(x0,y0) in {'N':(-0.15,2.55)}.items():
    pass
# HCO+ peak near N, S
m0=O.m0(-100,100)
def pk(x,y,r=1.5):
    return peak_near(O,m0,x,y,r)
N=pk(-0.15,2.55); S=pk(0.0,-2.0,2.5) if False else None
print('N peak',N)
Rr=lambda c,x0,y0,r1,r2,lo,hi:(lambda m:(c.cube[(c.v>=lo)&(c.v<hi)].sum(0)*c.dv)[m].sum())((np.hypot(c.X-x0,c.Y-y0)>=r1)&(np.hypot(c.X-x0,c.Y-y0)<r2))
for lo,hi in [(-140,-100),(-100,-60),(-60,-20),(-20,20),(20,60),(60,100),(100,140),(140,260),(-260,-140)]:
    print(lo,hi,[round(Rr(H,*N,0,r,lo,hi)/Rr(O,*N,0,r,lo,hi),2) for r in (0.5,1.0,2.0)])
