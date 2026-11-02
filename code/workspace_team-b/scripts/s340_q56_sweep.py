"""Q56: R32/10 (K km/s, RJ Tb) for HCN and HCO+ at N and S vs circular beam 1.4-3.5", Xb4d image-plane 3-2 (uv-incomplete) vs X961 1-0 (12m+ACA, 2015.1.00993).
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
out=[]; res={}
def P(t): print(t,flush=True); out.append(t)
for w in (100,150,350):
    P(f'\n== window |v|<{w} km/s: R32/10 = I(3-2 Xb4d image)/I(1-0), Tb K km/s')
    for th in (1.4,1.7,2.0,2.5,3.0,3.5):
        row=f'{th:.1f}" '
        for ln in ('HCN','HCO+'):
            R={}
            for rn in 'NS':
                a=I('32',ln,th,rn,w); b=I('10',ln,th,rn,w)
                if a is None or b is None: R[rn]=None; continue
                R[rn]=(a[0]/b[0],abs(a[0]/b[0])*np.hypot(a[1]/a[0],b[1]/b[0]),a[0],b[0])
            if None in R.values(): row+=f'{ln}: n/a (beam < native) | '; continue
            row+=f'{ln}: N {R["N"][0]:.3f} S {R["S"][0]:.3f} S/N {R["S"][0]/R["N"][0]:.2f} (I10 N {R["N"][3]:.1f} S {R["S"][3]:.1f}) | '
            res[f'{w}|{th}|{ln}']={rn:R[rn][:2] for rn in 'NS'}
        P(row)
json.dump(res,open('/workspace/out/q56/s340.json','w'),indent=1,default=float); open('/workspace/out/q56/s340.txt','w').write('\n'.join(out))
