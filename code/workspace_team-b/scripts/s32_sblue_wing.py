# Q1/Q10: HCN/HCO+ 3-2 in the Sblue blue-wing blob (S nucleus outflow side) vs core, Xb4d native cubes (beams 0.94x0.74 / 0.93x0.74, matched to 1%).
# Fixed aperture centre = S nucleus + (dE,dN)=(-0.2,-1.85)" (mean of the HCO+/HCN wing centroids in s11b_geometry.txt), r = 0.5, 0.75, 1.0, 1.5".
# Wing window -400..-150 km/s, core |v|<150 about VSYS (Sblue and S nucleus). Noise: scatter of the same-size aperture sums in the same window at 400 random
# line-free sky positions (r 6-18" of the midpoint, 2.5 beams from known emission peaks). Also the Xb4b (0.25") cube convolved to the Xb4d beam is not used here (see s11).
import sys; sys.path.insert(0,'/workspace/scripts')
from s30_matchedJ_ratios import *
import json
rng=np.random.default_rng(1)
o=run('Xb4d',3,None,windows=(350,))
S0=REGIONS['S'][:2]; dec=S0[1]-1.85/3600; ra=S0[0]-(-0.2)/3600/np.cos(np.deg2rad(S0[1]))*(-1)   # offset east=-0.2" -> ra decreases... East is +RA: -0.2" E => ra - 0.2"/cos
ra=S0[0]-0.2/3600/np.cos(np.deg2rad(S0[1]))
res={}; txt=['## Sblue / S-nucleus wing and core, Xb4d native, HCN/HCO+ (Jy km/s per beam; R in Jy; K-ratio = x1.0126)']
def win(R,a,b): m=(R['v']>=a)&(R['v']<=b); return R['c'][m].sum(0)*R['dv']
for name,(rr_,dd_,rra,rdec) in {'Sblue':(0,0,ra,dec),'S':(0,0,S0[0],S0[1]),'N':(0,0,*REGIONS['N'][:2])}.items():
    for wn,(a,b) in {'blue(-400,-150)':(-400,-150),'core(-150,150)':(-150,150),'red(150,400)':(150,400)}.items():
        for rad in (0.5,0.75,1.0,1.5):
            vals={}
            for mol in ('HCN','HCO+'):
                R=o[mol]; m=win(R,a,b); n=m.shape[0]; Y,X=np.mgrid[:n,:n]; px,py=R['pxy'](rra,rdec); ps=R['ps']
                ap=lambda x,y:np.hypot(X-x,Y-y)*ps<=rad
                s=m[ap(px,py)].sum()/ (1.1331*0.94*0.74/ps**2)
                # random positions
                cx,cy=R['pxy'](*CEN); rr=R['rr']; bright=R['m0_350']>5*R['err_350']
                from scipy.ndimage import binary_dilation
                bd=binary_dilation(bright,iterations=int(2.5*0.94/ps))
                ok=[];
                while len(ok)<400:
                    x,y=rng.uniform(0,n,2)
                    if 6<np.hypot(x-cx,y-cy)*ps<18 and not bd[int(y),int(x)]: ok.append((x,y))
                bg=np.array([m[ap(x,y)].sum()/(1.1331*0.94*0.74/ps**2) for x,y in ok])
                vals[mol]=(s,bg.std())
            (h,eh),(c,ec)=vals['HCN'],vals['HCO+']
            r=h/c if c>0 else np.nan; e=abs(r)*np.hypot(eh/h,ec/c) if h>0 and c>0 else np.nan
            res[f'{name}_{wn}_r{rad}']=dict(HCN=h,eHCN=eh,HCOp=c,eHCOp=ec,R=r,eR=e)
            txt.append(f'  {name:6s} {wn:16s} r={rad:4.2f}"  HCN {h:7.3f}+-{eh:.3f}  HCO+ {c:7.3f}+-{ec:.3f}  R={r:.2f}+-{e:.2f}')
json.dump(res,open('/workspace/out/s32_sblue_wing.json','w'),indent=1,default=float)
open('/workspace/out/s32_sblue_wing.txt','w').write('\n'.join(txt)); print('\n'.join(txt))
