import numpy as np
C=299792.458
def load(ln):
    z=np.load(f'/workspace/out/q7/cube_Xb4d_{ln}.npz'); c=z['cube']; v=z['v']-2780; 
    bm=(abs(v)>600)&(abs(v)<1500); c=c-np.median(c[bm],axis=0)
    dv=abs(np.median(np.diff(v))); m0=c[abs(v)<350].sum(0)*dv
    X,Y=np.meshgrid((z['ra']-156.96346)*np.cos(np.deg2rad(-43.90461))*3600,(z['de']+43.90461)*3600); return m0,X,Y
a,X,Y=load('HCN'); b,_,_=load('HCO+')
for nm,(x0,y0) in {'N':(0,2.63),'S':(-0.15,-2.64)}.items():
    r=np.hypot(X-x0,Y-y0)
    out=[]
    for r1,r2 in [(0,.5),(.5,1),(1,1.5),(1.5,2),(2,3)]:
        m=(r>=r1)&(r<r2); out.append(round(a[m].sum()/b[m].sum(),3))
    print(nm,out)
