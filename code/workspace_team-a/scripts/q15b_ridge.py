# Q15b: is the red-wing dense-gas ridge (~3.3" N of S, team-b p0017) real? Xb4d and Xb4b (both 1.05" circular beam), HCN and HCO+ 3-2.
# Windows (radio km/s, own rest frame; v0=2730 = line-core centre): red wing +155..+355 (2885-3085), blue mirror 2375-2575 (-355..-155; = Q13 blue lobe window),
#   team-b window +180..+320 (2910-3050) reported too. Aperture = disk r=1.0" (mean, Jy/beam km/s) as Q13g.
# (1) flux at the ridge (peak pixel in a 1.05" beam and r=1" mean), both lines, both arrays
# (2) blind search: local maxima of red map over field r<=7.5"; global p from sliding line-free 200-km/s windows (same noise treatment) -> trials-corrected sigma
# (3) core-matched / mirrored controls: D = red - blue(mirror) removes symmetric core leakage; scatter of D vs core brightness at grid points excl. ridge/lobe/N/S
# (4) 2D Gaussian centroid/size of the ridge; offsets from M18 OF-Red
# Outputs: out/q15/q15b.log, q15b_*.csv
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from scipy.signal import fftconvolve
from scipy.ndimage import maximum_filter
from scipy.optimize import curve_fit
from scipy.stats import norm
from q3_pos import POS
from q3b_lib import FREE
O='/workspace/out/q15/'; log=open(O+'q15b.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
rob=lambda a: 1.4826*np.median(np.abs(a-np.median(a)))
V0=2730.; RMAX=7.5; RAD=1.0
WR=(2885,3085); WB=(2375,2575); WT=(2910,3050); WCORE=(2575,2875)
def load(ous):
    z=np.load(O+f'cube_{ous}.npz'); return {k:z[k] for k in z.files}
def disk_mean(M,dx,dy,pix):
    n=int(np.ceil(RAD/pix)); yy,xx=np.mgrid[-n:n+1,-n:n+1]; k=(np.hypot(xx,yy)*pix<=RAD).astype(float); k/=k.sum()
    return fftconvolve(M,k,mode='same')
def win(C,sp,a,b):
    s,v=(C['HCN'],C['vH']) if sp=='HCN' else (C['HCO'],C['vO']); m=(v>=a)&(v<b); return s[m].sum(0)*C['dv']
def ap_at(M,C,x,y,rad=RAD):
    m=np.hypot(C['dx']-x,C['dy']-y)<=rad; return M[m].mean()
Sx,Sy=POS['S']; Nx,Ny=POS['N']; RIDGE0=(Sx+0.0,Sy+3.3); LOBE0=POS['Soutflow']
res={}
for ous in ('Xb4d','Xb4b'):
    C=load(ous)
    if ous=='Xb4b':   # resample to 0.15" grid by block-mean of the 1.05"-smoothed 0.047" cube (strided sampling, beam >> pixel)
        st=3; C={**C,'dx':C['dx'][::st,::st],'dy':C['dy'][::st,::st],'HCN':C['HCN'][:,::st,::st],'HCO':C['HCO'][:,::st,::st],'pix':C['pix']*st}
    dx,dy,pix=C['dx'],C['dy'],C['pix']; field=np.hypot(dx,dy)<=RMAX
    n=int(np.ceil(RAD/pix))
    pr(f'=== {ous}: pix {pix:.3f}", dv {C["dv"]:.2f} km/s; grid {dx.shape}')
    for sp in ('HCN','HCO'):
        Mr=win(C,sp,*WR); Mb=win(C,sp,*WB); Mt=win(C,sp,*WT); Mc=win(C,sp,*WCORE)
        # per-channel rms from line-free channels (formal) and window noise from random-aperture scatter in line-free windows
        sr=disk_mean(Mr,dx,dy,pix); sb=disk_mean(Mb,dx,dy,pix); st_=disk_mean(Mt,dx,dy,pix)
        # find ridge peak in red map within 2" of nominal
        near=np.hypot(dx-RIDGE0[0],dy-RIDGE0[1])<=2.0
        i=np.unravel_index(np.argmax(np.where(near,sr,-9)),sr.shape); xp,yp=dx[i],dy[i]
        near2=np.hypot(dx-LOBE0[0],dy-LOBE0[1])<=2.0; j=np.unravel_index(np.argmax(np.where(near2,sb,-9)),sb.shape); xl,yl=dx[j],dy[j]
        # sliding null windows (width 200, line-free)
        v=C['vH'] if sp=='HCN' else C['vO']; fr=FREE['HCN' if sp=='HCN' else 'HCO+']
        starts=[a for a in np.arange(1300,2100,25) if any(a>=r0 and a+200<=r1 for r0,r1 in fr)]
        nul=[]; nulmax=[]
        for a in starts:
            W=win(C,sp,a,a+200); s=disk_mean(W,dx,dy,pix)
            mx=(s==maximum_filter(s,size=2*n+1))&field; nul.append(s[field]); nulmax.append(s[mx].max())
        nul=np.concatenate(nul); sg=nul.std()
        pk=(sr==maximum_filter(sr,size=2*n+1))&field; pkv=sr[pk]; rank=(pkv>sr[i]+1e-9).sum()+1
        pkb=(sb==maximum_filter(sb,size=2*n+1))&field; rankb=(sb[pkb]>sb[j]+1e-9).sum()+1
        Nind=np.pi*RMAX**2/(1.1331*1.05**2)
        z=sr[i]/sg; p1=norm.sf(z); pg=1-(1-p1)**Nind
        zb=sb[j]/sg; pgb=1-(1-norm.sf(zb))**Nind
        pr(f'{sp}: line-free sliding-window null (r=1" mean, n_win={len(starts)}): sd {sg:.4f} Jy/b km/s, max over field per window {np.round(nulmax,3).tolist()}')
        pr(f'  red ridge peak at ({xp:.2f},{yp:.2f}) = {xp-Sx:+.2f},{yp-Sy:+.2f} from S: r=1" mean {sr[i]:.3f}; peak pix (1.05" beam, native) {Mr[i]:.3f}; team-b window {st_[i]:.3f}; pointwise z={z:.1f}; rank {rank}/{pk.sum()} maxima; global p={pg:.2e} ({norm.isf(pg):.1f} sigma one-sided, N_ind={Nind:.0f})')
        pr(f'  blue mirror at the same position: {sb[i]:.3f}; red at lobe position: {sr[j]:.3f}')
        pr(f'  blue lobe peak at ({xl:.2f},{yl:.2f}) = {xl-Sx:+.2f},{yl-Sy:+.2f} from S: r=1" mean {sb[j]:.3f}; z={zb:.1f}; rank {rankb}/{pkb.sum()}; global {norm.isf(pgb):.1f} sigma; ridge/lobe flux ratio {sr[i]/sb[j]:.2f}')
        top=np.argsort(-np.where(pk,sr,-9).ravel())[:6]
        pr('  top red maxima (dE,dN abs; sigma units):',[(round(dx.ravel()[t],1),round(dy.ravel()[t],1),round(sr.ravel()[t]/sg,1)) for t in top])
        res[(ous,sp)]=dict(sr=sr,sb=sb,Mr=Mr,Mb=Mb,Mc=Mc,sg=sg,i=i,j=j,C=C)
    # combined HCN+HCO+ red-wing map (sum) with sum noise
    sr=res[(ous,'HCN')]['sr']+res[(ous,'HCO')]['sr']; sb=res[(ous,'HCN')]['sb']+res[(ous,'HCO')]['sb']
    nul=[]
    for a in [a for a in np.arange(1300,2100,25) if any(a>=r0 and a+200<=r1 for r0,r1 in FREE['HCN']) and any(a>=r0 and a+200<=r1 for r0,r1 in FREE['HCO+'])]:
        nul.append((disk_mean(win(C,'HCN',a,a+200),dx,dy,pix)+disk_mean(win(C,'HCO',a,a+200),dx,dy,pix))[field])
    sg=np.concatenate(nul).std(); Nind=np.pi*RMAX**2/(1.1331*1.05**2)
    pk=(sr==maximum_filter(sr,size=2*n+1))&field
    near=np.hypot(dx-RIDGE0[0],dy-RIDGE0[1])<=2.0; i=np.unravel_index(np.argmax(np.where(near,sr,-9)),sr.shape); rk=(sr[pk]>sr[i]+1e-9).sum()+1
    z=sr[i]/sg; pg=1-(1-norm.sf(z))**Nind
    pr(f'HCN+HCO+ summed red map ridge: {sr[i]:.3f} Jy/b km/s, null sd {sg:.4f}, z={z:.1f}, rank {rk}/{pk.sum()}, global {norm.isf(pg):.1f} sigma; blue mirror at ridge {sb[i]:.3f}')
    res[(ous,'sum')]=dict(sr=sr,sb=sb,sg=sg,i=i,C=C)
import pickle; pickle.dump({k:{kk:vv for kk,vv in v.items() if kk!='C'} for k,v in res.items()},open(O+'q15b_maps.pkl','wb'))
# (3) mirrored / core-matched controls, Xb4d
for sp in ('HCN','HCO'):
    r=res[('Xb4d',sp)]; C=r['C']; dx,dy=C['dx'],C['dy']
    sr,sb=r['sr'],r['sb']; core=disk_mean(r['Mc'],dx,dy,C['pix']); D=sr-sb
    field=np.hypot(dx,dy)<=RMAX
    ex=np.zeros_like(field)
    for p in (RIDGE0,LOBE0,POS['N'],POS['S']): ex|=np.hypot(dx-p[0],dy-p[1])<1.5
    gx,gy=np.meshgrid(np.arange(-7,7.1,1.0),np.arange(-7,7.1,1.0)); pts=[(x,y) for x,y in zip(gx.ravel(),gy.ravel()) if np.hypot(x,y)<=RMAX and not (min(np.hypot(x-a,y-b) for a,b in (RIDGE0,LOBE0,POS['N'],POS['S']))<1.5)]
    cc=np.array([ap_at(core,C,x,y) for x,y in pts]); dd=np.array([ap_at(D,C,x,y) for x,y in pts]); rr=np.array([ap_at(sr,C,x,y) for x,y in pts]); bb=np.array([ap_at(sb,C,x,y) for x,y in pts])
    ci=ap_at(core,C,*RIDGE0); Dri=ap_at(D,C,*RIDGE0)
    pr(f'--- {sp} Xb4d controls: ridge core(2575-2875) {ci:.2f}, D=red-blue at ridge {Dri:.3f}; {len(pts)} grid points (1" spacing, >1.5" from ridge/lobe/N/S)')
    sel=(cc>0.5*ci)&(cc<2*ci)
    pr(f'    core-matched (0.5-2x ridge core) control points N={sel.sum()}: D mean {dd[sel].mean() if sel.any() else np.nan:.3f} sd {dd[sel].std() if sel.sum()>1 else np.nan:.3f}; red mean {rr[sel].mean() if sel.any() else np.nan:.3f}; blue mean {bb[sel].mean() if sel.any() else np.nan:.3f}')
    sel2=cc>0.25*ci
    pr(f'    brighter-core control points (>0.25x ridge core) N={sel2.sum()}: D mean {dd[sel2].mean():.3f} sd {dd[sel2].std():.3f} max {dd[sel2].max():.3f}; red max {rr[sel2].max():.3f}')
    A=np.vstack([cc,np.ones_like(cc)]).T
    for nm,y in (('red',rr),('blue',bb),('D',dd)):
        co=np.linalg.lstsq(A,y,rcond=None)[0]; res_=y-A@co
        pr(f'    {nm} wing vs core slope {co[0]:.4f} (+/-) offset {co[1]:.3f}; resid sd {res_.std():.3f}; predicted at ridge core {co[0]*ci+co[1]:.3f} vs obs {ap_at({"red":sr,"blue":sb,"D":D}[nm],C,*RIDGE0):.3f}')
    co=np.linalg.lstsq(A,rr,rcond=None)[0]; sd=(rr-A@co).std()
    pr(f'    ridge excess over core-scaled prediction: {(ap_at(sr,C,*RIDGE0)-(co[0]*ci+co[1]))/sd:.1f} sigma_scatter (N={len(pts)})')
    # symmetric leakage test: blue wing at ridge position / red wing at lobe
    pr(f'    leakage mirror: blue wing at ridge {ap_at(sb,C,*RIDGE0):.3f}; red wing at lobe {ap_at(sr,C,*LOBE0):.3f}; blue at lobe {ap_at(sb,C,*LOBE0):.3f}; red at ridge {ap_at(sr,C,*RIDGE0):.3f}')
# (4) 2D Gaussian (beam-convolved, free amp, centre, intrinsic sigma, fixed 1.05" FWHM beam), red wing HCN+HCO+ summed, Xb4d, fit within 4" of ridge
def g2(xy,A,x0,y0,s,b):
    x,y=xy; st=np.sqrt(s**2+(1.05/2.3548)**2); return A*np.exp(-((x-x0)**2+(y-y0)**2)/(2*st**2))+b
for sp in ('HCN','HCO','sum'):
    r=res[('Xb4d',sp)]; C=r['C']; dx,dy=C['dx'],C['dy']
    M=r['sr'] if sp=='sum' else r['Mr']
    if sp=='sum': M=res[('Xb4d','HCN')]['Mr']+res[('Xb4d','HCO')]['Mr']
    m=np.hypot(dx-RIDGE0[0],dy-RIDGE0[1])<=3.0
    # point-sampled (non-smoothed) map; pixels correlated by beam -> use sigma from line-free pixel scatter
    sig=res[('Xb4d','HCN' if sp!='HCO' else 'HCO')]['sg']*0+0.03
    try:
        p,cv=curve_fit(g2,(dx[m],dy[m]),M[m],p0=[0.2,RIDGE0[0],RIDGE0[1],0.5,0],bounds=([0,-3,-6,0,-0.1],[2,3,2,2,0.1]))
        e=np.sqrt(np.diag(cv)); n_eff=m.sum()/(1.1331*1.05**2/0.15**2)
        e=e*np.sqrt(m.sum()/n_eff)/1.0
        pr(f'2D gaussian ridge red wing ({sp}, Xb4d): centre dE={p[1]-Sx:+.2f}+-{e[1]:.2f} dN={p[2]-Sy:+.2f}+-{e[2]:.2f} from S (=({p[1]:+.2f},{p[2]:+.2f}) from phase centre; dN from N {p[2]-Ny:+.2f}); intrinsic sigma {p[3]:.2f}+-{e[3]:.2f}" (FWHM {2.3548*p[3]:.2f}"); peak {p[0]:.3f}; offset from M18 OF-Red ({POS["M18 OF-Red"][0]:+.2f},{POS["M18 OF-Red"][1]:+.2f}): d=({p[1]-POS["M18 OF-Red"][0]:+.2f},{p[2]-POS["M18 OF-Red"][1]:+.2f})')
    except Exception as ex_: pr('fit failed',sp,ex_)
