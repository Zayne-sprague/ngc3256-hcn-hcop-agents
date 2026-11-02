"""[helper copy of s340 functions] Q56: R32/10 (K km/s, RJ Tb) for HCN and HCO+ at N and S vs circular beam 1.4-3.5", Xb4d image-plane 3-2 (uv-incomplete) vs X961 1-0 (12m+ACA, 2015.1.00993).
Same pixel positions, velocity windows (+-100, +-150, +-350 km/s), linear baseline, as s332 (Q53). Output out/q56/s340.txt/json.
usage: python3 -I s340_q56_sweep.py"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
import s320_q52_chem_ratios as m
from s180_q29_load import to_tb_conv
TB={}
def spec(J,ln,th,rn):
    k=(J,ln,th)
    if k not in TB:
        path,f0=m.LINES[J][ln]; cut,nu,pix,bm,pxy=m.getcube(path,th)
        if th<bm[0]*0.999: TB[k]=None
        else: TB[k]=(to_tb_conv(cut,nu,pix,bm,th,float(np.mean(nu))),nu,pxy)
    if TB[k] is None: return None
    Tb,nu,pxy=TB[k]; path,f0=m.LINES[J][ln]; pos,v0=m.POS[rn]; x,y=pxy(pos)
    sp=Tb[:,int(round(y)),int(round(x))].astype(float); v=m.vel(nu,f0,v0); o=np.argsort(v); v,sp=v[o],sp[o]; nuo=nu[o]
    oth={'HCN':'HCO+','HCO+':'HCN'}[ln]; g0=m.LINES[J][oth][1]; vg=m.vel(nuo,g0,v0)
    free=(abs(v)>600)&(abs(v)<1200)&~((vg>-650)&(vg<650))
    p=np.polyfit(v[free],sp[free],1); s=sp-np.polyval(p,v); sig=np.std(s[free]); return v,s,sig,np.gradient(v)
def integ(v,s,dv,a,b):
    e=np.concatenate([[v[0]-dv[0]/2],(v[1:]+v[:-1])/2,[v[-1]+dv[-1]/2]]); F=np.concatenate([[0],np.cumsum(s*dv)]); return np.interp(b,e,F)-np.interp(a,e,F)
def I(J,ln,th,rn,w):
    r=spec(J,ln,th,rn)
    if r is None: return None
    v,s,sig,dv=r; return integ(v,s,dv,-w,w),sig*np.sqrt(2*w*np.mean(dv[abs(v)<300]))
