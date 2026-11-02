# Q15h: kinematic consequences of the three origins. Core centroid velocities at S, N, mid-Bridge, ridge; wing velocities; profile of wing flux along the N-S axis (red +230..+470 vs blue -470..-230); dynamical times.
# Uses the Q15e baseline C (blue free + red edge >3290, order 1), Xb4d 1.05". D=40 Mpc is NOT measured here: 1"=194 pc (assumed, as elsewhere in the report; check q-lib).
import sys; sys.path.insert(0,'/workspace/scripts')
from q13a_xb4_compare import *
from q3_pos import POS
O='/workspace/out/q15/'; log=open(O+'q15h.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
Sx,Sy=POS['S']; N=POS['N']; RIDGE=(0.3,Sy+3.3); LOBE=POS['Soutflow']; PCPA=194.
def red(sp):
    c=vcube(sp,'Xb4d',1.05,9.); v=c['v']; a=c['a']
    m=((v>=FREE[sp][0][0])&(v<=FREE[sp][-1][1]))|(v>=3290)
    X=np.vander(v-v[m].mean(),2); A=a[m].reshape(m.sum(),-1); coef=np.linalg.lstsq(X[m],A,rcond=None)[0]
    return c,a-(X@coef).reshape(len(v),*a.shape[1:]),v
for sp in ('HCN','HCO+'):
    c,r,v=red(sp); dx,dy,dv=c['dx'],c['dy'],c['dv']
    pr(f'## {sp}: core (2575-2925) flux-weighted centroid velocity and wing mean velocity (r=0.5" apertures)')
    for nm,(x,y) in (('S',POS['S']),('N',N),('mid S-N (0.15,0.0)',(0.15,0.0)),('Bridge nominal',POS['Bridge']),('ridge',RIDGE),('lobe',LOBE)):
        ap=np.hypot(dx-x,dy-y)<=0.5; s=r[:,ap].mean(1); w=(v>=2575)&(v<2925)
        vc=(s[w]*v[w]).sum()/s[w].sum(); wr=(v>=2960)&(v<3300); wb=(v>=2160)&(v<2500)
        pr(f'{nm:20s} core centroid {vc:.0f} km/s; core peak {1e3*s[w].max():.1f} mJy/b; red-wing(2960-3300) mean v {((s[wr]*v[wr]).sum()/s[wr].sum()-vc) if s[wr].sum()>0 else np.nan:+.0f} rel core (flux {s[wr].sum()*dv:+.3f}); blue-wing(2160-2500) flux {s[wb].sum()*dv:+.3f}')
    # profile along N-S axis, strip |dx-0.15|<0.6, rows 0.5" wide
    pr(f'## {sp}: wing flux (Jy/b km/s, mean of strip) vs dN from S: red +230..+470 | blue -470..-230')
    for yo in np.arange(-5.0,8.1,0.75):
        y=Sy+yo; m=(np.abs(dx-0.0)<=0.6)&(np.abs(dy-y)<=0.375)
        R=r[(v>=2960)&(v<3200)][:,m].mean(1).sum()*dv; B=r[(v>=2260)&(v<2500)][:,m].mean(1).sum()*dv
        pr(f'  dN_S={yo:+5.2f}": red {R:+.3f}  blue {B:+.3f}')
pr(f'Geometry: lobe at {np.hypot(LOBE[0]-Sx,LOBE[1]-Sy):.2f}" = {np.hypot(LOBE[0]-Sx,LOBE[1]-Sy)*PCPA:.0f} pc S of S; ridge at {np.hypot(RIDGE[0]-Sx,RIDGE[1]-Sy):.2f}" = {np.hypot(RIDGE[0]-Sx,RIDGE[1]-Sy)*PCPA:.0f} pc N of S, {np.hypot(RIDGE[0]-N[0],RIDGE[1]-N[1]):.2f}" = {np.hypot(RIDGE[0]-N[0],RIDGE[1]-N[1])*PCPA:.0f} pc S of N; S-N separation {np.hypot(*(np.array(N)-np.array(POS["S"]))):.2f}" = {np.hypot(*(np.array(N)-np.array(POS["S"])))*PCPA:.0f} pc')
for nm,d,vv in (('blue lobe',1.93,262),('red ridge',3.4,290)):
    pr(f'{nm}: projected dynamical time d/v = {d*PCPA*3.0857e13/vv/3.156e13:.2f} Myr (no inclination correction)')
