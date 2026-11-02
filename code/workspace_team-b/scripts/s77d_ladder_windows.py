# Q2 round 2, fix 1 (inputs to the non-LTE fit): wing fluxes of HCN and HCO+ at the lobe in the 3-2 line (Xb4d, matched to the Michiyama+2018 South-OF-Blue window/aperture) and the
# ladder constraints from the low-resolution cubes (Xb55 J=2-1, Xb4f J=3-2, Xb47 J=4-3) convolved to a common 10" circular beam, window -405..-150 and -360..-180, point value (Jy per 10" beam = total flux for a source compact vs 10") at the S+lobe midpoint.
# Per-pixel baseline median(-900..-450), flux rescaled by beam-area ratio, noise from random field apertures. The aggregate low-res wing contains S-nucleus wings + lobe + Bridge: UPPER limits for the lobe.
# Michiyama+2018 (arXiv:1810.04821) Table 3: HCN(1-0) 0.07+-0.02, HCO+(1-0) 0.10+-0.02 Jy km/s (South-OF-Blue, [-360,-180], inner 2.1"); text sect 3.1: 0.12 / 0.20 Jy km/s for the whole blue wing map.
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.ndimage import map_coordinates
OUT='/workspace/out/q2r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); rng=np.random.default_rng(3)
def wflux(v,dv,sp,lo,hi):
    ov=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return float((sp*ov).sum()*dv)
res={}
M18=(10+27/60+51.17/3600)*15, -(43+54/60+21.43/3600)          # South-OF-Blue (Table 1)
N18=((10+27/60+51.22/3600)*15,-(43+54/60+19.20/3600))          # M18 South nucleus
# --- (1) 3-2 in the M18 position/window
off=((M18[0]-N18[0])*np.cos(np.deg2rad(N18[1]))*3600,(M18[1]-N18[1])*3600); P(f'M18 South-OF-Blue is at dE={off[0]:+.2f}", dN={off[1]:+.2f}" from the M18 South nucleus')
for fw in (1.05,2.0):
  for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,fw,12.0); v=R['v']; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0); fm=field_pixels(R); ys,xs=np.where(fm); Y,X=np.mgrid[:R['n'],:R['n']]
    for tag,pos in (('M18pos',offset_pos(S0,off[0],off[1])),('lobe',LOBE)):
        for rad in (1.05,1.5,2.5):
            F=wflux(v,R['dv'],aspec(c,R['ba'],apmask(R['n'],R['ps'],R['pxy'](*pos),rad)),-360,-180)
            nn=[]
            while len(nn)<40:
                i=rng.integers(len(ys)); m=np.hypot(X-xs[i],Y-ys[i])*R['ps']<=rad
                if fm[m].all(): nn.append(wflux(v,R['dv'],aspec(c,R['ba'],m),-360,-180))
            res[f'32|{fw}|{mol}|{tag}|r{rad}']=(F,float(np.std(nn)))
            P(f'3-2 Xb4d@{fw}" {mol:5s} {tag:6s} r={rad}" window -360..-180: {F:.3f} +- {np.std(nn):.3f} Jy km/s')
# --- (2) ladder at 10"
BEAM=10.0
LAD={('Xb55',2):('spw16','spw18'),('Xb4f',3):('spw18','spw16'),('Xb4d',3):('spw27','spw25'),('Xb47',4):('spw16','spw18')}
mid=offset_pos(S0,0.1,-0.95)
for (ous,J),(sh,so) in LAD.items():
  for mol,spw in (('HCN',sh),('HCO+',so)):
    c,v,ps,wc,x0,y0,b0=load_cube(ous,spw,F0[mol][J],-1000,700,half=24.0 if ous!='Xb55' else 24.0)
    c=conv(c,ps,b0,(BEAM,BEAM,0.0)); c=c-np.median(c[(v>-900)&(v<-450)],axis=0); n=c.shape[1]; ba=1.1331*BEAM**2/ps**2; dv=abs(np.median(np.diff(v))); pxy=pxy_fn(wc,x0,y0)
    m0=c[np.abs(v)<350].sum(0)*dv; Y,X=np.mgrid[:n,:n]; cx,cy=pxy(*CEN)
    # field apertures r=3" away from the centre (>=12" from the mid-point), no brightness cut possible at 10" beam -> use ring 12-20" around centre
    fld=[(x,y) for x,y in zip(X.ravel(),Y.ravel()) if 12<np.hypot(x-cx,y-cy)*ps<20 and 3/ps<x<n-3/ps and 3/ps<y<n-3/ps]
    for lo,hi in ((-405,-150),(-360,-180),(-205,95)):
        px_,py_=pxy(*mid); pt=lambda x,y: np.array([map_coordinates(c[i],[[y],[x]],order=1)[0] for i in range(len(v))])   # Jy/(10" beam) at a point = total flux of a source compact vs 10"
        F=wflux(v,dv,pt(px_,py_),lo,hi)
        idx=rng.choice(len(fld),size=min(60,len(fld)),replace=False); nn=[wflux(v,dv,pt(*fld[i]),lo,hi) for i in idx]
        res[f'lad|{ous}|J{J}|{mol}|{lo}_{hi}']=(F,float(np.std(nn)),float(np.mean(nn)))
        P(f'LADDER {ous} J={J}-{J-1} {mol:5s} @10" point value at S+lobe midpoint (= total flux of a source compact vs 10"), window {lo}..{hi}: {F:.3f} +- {np.std(nn):.3f} Jy km/s (field mean {np.mean(nn):+.3f}); channel width {dv:.1f} km/s; native beam {b0[0]:.1f}x{b0[1]:.1f}')
json.dump(res,open(OUT+'/s77d.json','w'),indent=1); open(OUT+'/s77d.txt','w').write('\n'.join(txt))
