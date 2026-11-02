import sys,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s320_q52_chem_ratios as m
from s180_q29_load import to_tb_conv
def prof(J,ln,th,rn):
    path,f0=m.LINES[J][ln]; pos,v0=m.POS[rn]
    cut,nu,pix,bm,pxy=m.getcube(path,th); v=m.vel(nu,f0,v0)
    sel=(v>-1500)&(v<1500)
    Tb=to_tb_conv(cut[sel],nu[sel],pix,bm,th,float(np.mean(nu))); vs=v[sel]
    x,y=pxy(pos); sp=Tb[:,int(round(y)),int(round(x))]
    free=(abs(vs)>600)&(abs(vs)<1200)
    p=np.polyfit(vs[free],sp[free],1); return vs,sp-np.polyval(p,vs)
edges=[-350,-200,-100,-50,0,50,100,200,350]
for rn in 'NS':
    P={l:prof('10',l,2.5,rn) for l in ['HCN','HCO+','HNC']}
    vs=np.linspace(-350,350,1401); P={k:(vs,np.interp(vs,*(a[0],a[1]))) if np.all(np.diff(a[0])>0) else (vs,np.interp(vs,a[0][::-1],a[1][::-1])) for k,a in P.items()}; print(rn,'v_sys',m.POS[rn][1])
    for a,b in zip(edges[:-1],edges[1:]):
        w=(vs>=a)&(vs<b); dv=0.5
        h=P['HCO+'][1][w].sum()*dv
        print(f'{a:5d}..{b:4d} HCO+ {h:7.1f} HCN/HCO+ {P["HCN"][1][w].sum()*dv/h:6.2f} HNC/HCO+ {P["HNC"][1][w].sum()*dv/h:6.2f}')
