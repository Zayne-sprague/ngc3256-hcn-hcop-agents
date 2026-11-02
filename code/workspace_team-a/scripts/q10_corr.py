# Q10: per-pixel/annulus HCN/HCO+ 3-2 vs line-free continuum (from the cube baselines) at matched beams, N and S separately.
# python3 -I q10_corr.py  (needs out/q7/cube_*.npz from q7_cache.py).  Output: out/q10/q10_corr.log, q10_corr.json
import sys,json,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
from scipy.stats import spearmanr,pearsonr
from scipy.signal import fftconvolve
from q7_lib import *
rng=np.random.default_rng(10)
out={}
def raw_cont(ous,ln,mode='lin'):
    """line-free continuum (Jy/beam, pbcor) at v=VSYS from the baseline channels, and the raw cube"""
    z=np.load(f'/workspace/out/q7/cube_{ous}_{ln}.npz'); cube=z['cube'].astype(np.float64); v=z['v']-VSYS
    m=np.zeros(len(v),bool)
    for a,b in BASE: m|=(v>=a)&(v<=b)
    if mode=='lin': X=np.vstack([np.ones(m.sum()),v[m]/500.]).T; c=np.tensordot(np.linalg.pinv(X),cube[m],axes=(1,0)); return c[0],cube,v,m,c
    c0=cube[m].mean(0); return c0,cube,v,m,None
C={}
for ous in ('Xb4b','Xb4d'):
    for ln in ('HCN','HCO+'): C[(ous,ln)]=Cube(ous,ln)
def cont_map(ous):
    # archive 4-spw continuum (257.8 GHz, line-free channels selected by the pipeline), pbcor, cropped to the cube cache grid
    from astropy.io import fits
    z=np.load(f'/workspace/out/q7/cube_{ous}_HCN.npz'); x0,y0=int(z['x0']),int(z['y0']); ny,nx=z['cube'].shape[1:]
    d=fits.getdata(f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')[0,0]
    return np.nan_to_num(d[y0:y0+ny,x0:x0+nx].astype(np.float64))
def cdiff(ous):
    a,b=raw_cont(ous,'HCN')[0],raw_cont(ous,'HCO+')[0]; return a,b
REG={'N':(0.0,2.63),'S':(-0.15,-2.64)}
def null_noise(c,w,tb,box=(-5.5,5.5)):
    n=smooth(null_m0(c,w),c,tb); return n
exec(open('/workspace/scripts/q7_profiles.py').read().split('def null_m0')[1].join(['def null_m0','']).split('def centres')[0]) if False else None
# null window helper (same as q7_profiles)
NULLW=(-760,-540)
def null_m0(c,w,wn=NULLW):
    m=np.zeros(len(c.v),bool)
    for a,b in BASE: m|=(c.v>=a)&(c.v<=b)
    m&=~((c.v>=wn[0]-10)&(c.v<=wn[1]+10)); raw=c.cube
    X=np.vstack([np.ones(m.sum()),c.v[m]/500.]).T; P=np.linalg.pinv(X); cc=np.tensordot(P,raw[m],axes=(1,0))
    k=(c.v>=wn[0])&(c.v<=wn[1]); Xk=np.vstack([np.ones(k.sum()),c.v[k]/500.]).T
    res=raw[k]-np.tensordot(Xk,cc,axes=(1,0)); return res.sum(0)*c.dv*np.sqrt((w[1]-w[0])/(wn[1]-wn[0]))
def mad(a): a=a[np.isfinite(a)]; return 1.4826*np.median(np.abs(a-np.median(a)))
CASES={'Xb4b@0.30':('Xb4b',(0.30,0.30,0.)),'Xb4b@1.05':('Xb4b',(1.05,1.05,0.)),'Xb4d@1.05':('Xb4d',(1.05,1.05,0.))}
W=(-350,350)
log=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); log.append(s)
def smooth_c(img,c,tb):  # continuum map (Jy/beam at native beam) -> target beam
    return smooth(img,c,tb)
res={}
for case,(ous,tb) in CASES.items():
    H,O=C[(ous,'HCN')],C[(ous,'HCO+')]
    c0=cont_map(ous); cH,cO=cdiff(ous)
    c0_=c0
    Sc=smooth_c(c0,H,tb)*1e3     # mJy/beam(tb)
    mH,mO=smooth(H.m0(*W),H,tb),smooth(O.m0(*W),O,tb)
    sH,sO=mad(smooth(null_m0(H,W),H,tb)),mad(smooth(null_m0(O,W),O,tb))
    # continuum noise: half-difference of the independent HCN-cube and HCO+-cube continuum estimates, smoothed
    sC=mad(Sc[np.hypot(H.X-0,H.Y-0)>4.5])
    P(f'\n=== {case}: beam {tb}; noise: HCN m0 {sH*1e3:.2f}, HCO+ m0 {sO*1e3:.2f} mJy/b km/s; cont {sC:.3f} mJy/b (MAD of the smoothed continuum map beyond 4.5 arcsec of the centre)')
    R=mH/np.where(mO>0,mO,np.nan)
    for reg,(x0,y0) in REG.items():
        pk=peak_near(O,mO,x0,y0); r=np.hypot(O.X-pk[0],O.Y-pk[1])
        cpk=peak_near(O,Sc,x0,y0,0.6)
        P(f'{reg}: HCO+ peak {pk[0]:+.2f},{pk[1]:+.2f}; continuum peak {cpk[0]:+.2f},{cpk[1]:+.2f} ({np.hypot(cpk[0]-pk[0],cpk[1]-pk[1]):.2f}" apart); Sc(peak)={Sc[np.unravel_index(np.argmax(np.where(np.hypot(O.X-x0,O.Y-y0)<0.6,Sc,-1e9)),Sc.shape)]:.3f} mJy/b')
        # annulus profile: mean cont, mean line, R
        edges=[0,0.15,0.25,0.5,1,1.5,2,3] if tb[0]<0.5 else [0,0.5,1,1.5,2,3]
        prof=[]
        for a,b in zip(edges[:-1],edges[1:]):
            m=(r>=a)&(r<b); prof.append((a,b,Sc[m].mean(),mH[m].sum()/mO[m].sum(),mO[m].mean()))
        P(' annulus  <Sc> mJy/b  R=HCN/HCO+  <HCO+> Jy/b km/s')
        for a,b,s,rr,io in prof: P(f'  {a:.2f}-{b:.2f}  {s:7.3f}  {rr:6.3f}  {io:7.3f}')
        pa=np.array([[p[2],p[3]] for p in prof])
        if len(pa)>=4:
            rho=spearmanr(pa[:,0],pa[:,1]); P(f' annulus-level Spearman(Sc,R) rho={rho[0]:.2f} p={rho[1]:.3g} n={len(pa)}  (annuli are not independent: descriptive)')
        # per-pixel on independent (beam-spaced) grid
        step=max(1,int(round(tb[0]/H.pix)))
        for rmax in (1.5,2.5):
          sel=(r<rmax)&(mH>3*sH)&(mO>3*sO)&(Sc>3*sC)
          # exclude the other nucleus: only pixels closer to this nucleus
          other=[v for k,v in REG.items() if k!=reg][0]; ob=peak_near(O,mO,*other)
          sel&=(r<np.hypot(O.X-ob[0],O.Y-ob[1]))
          best=None; rs=[]
          for ox in range(step):
            for oy in range(step):
              g=np.zeros_like(sel); g[oy::step,ox::step]=True; s=sel&g
              if s.sum()>=6:
                  rs.append((spearmanr(Sc[s],R[s])[0],spearmanr(mO[s],R[s])[0],s.sum()))
          if not rs: P(f' r<{rmax}: too few'); continue
          rs=np.array(rs)
          g=np.zeros_like(sel); g[0::step,0::step]=True; s=sel&g
          # bootstrap on a single grid
          x,y=Sc[s],R[s]; bs=[]
          for _ in range(2000):
              i=rng.integers(0,len(x),len(x)); bs.append(spearmanr(x[i],y[i])[0])
          # partial correlation of R with Sc controlling for HCO+ brightness and radius
          def resid(a,*cs):
              X=np.vstack([np.ones(len(a))]+list(cs)).T; return a-X@np.linalg.lstsq(X,a,rcond=None)[0]
          rr_=r[s]; pc_cont=pearsonr(resid(y,mO[s],rr_),resid(x,mO[s],rr_))[0]
          pc_hco=pearsonr(resid(y,x,rr_),resid(mO[s],x,rr_))[0]
          # R in log; also slope dR/dSc
          sl=np.polyfit(x,y,1)
          P(f' r<{rmax}": n_indep={s.sum()}; Spearman(Sc,R)={spearmanr(x,y)[0]:+.2f} [boot 16-84: {np.percentile(bs,16):+.2f},{np.percentile(bs,84):+.2f}]; offset-grid mean {rs[:,0].mean():+.2f}; Spearman(HCO+,R)={spearmanr(mO[s],y)[0]:+.2f}; partial(R,Sc|HCO+,r)={pc_cont:+.2f}; partial(R,HCO+|Sc,r)={pc_hco:+.2f}; slope dR/dSc={sl[0]:.3f} per mJy/b; Sc range {x.min():.2f}-{x.max():.2f}, R range {y.min():.2f}-{y.max():.2f}')
          res[(case,reg,rmax)]=dict(n=int(s.sum()),rho=float(spearmanr(x,y)[0]),boot=[float(np.percentile(bs,16)),float(np.percentile(bs,84))],pc_cont=float(pc_cont),pc_hco=float(pc_hco),slope=float(sl[0]))
    # --- continuum subtraction effects at the N peak (0.30 beam only)
    if case=='Xb4b@0.30':
        for reg,(x0,y0) in REG.items():
            pk=peak_near(O,mO,x0,y0); 
            for mode in ('lin','const','none'):
                Fl={}
                for ln in ('HCN','HCO+'):
                    c=C[(ous,ln)]; z=raw_cont(ous,ln,mode)
                    cube=z[1]; v=z[2]
                    if mode=='lin': base=0*cube; b=np.tensordot(np.vstack([np.ones(len(v)),v/500.]).T,z[4],axes=(1,0)); sub=cube-b
                    elif mode=='const': sub=cube-z[0][None]
                    else: sub=cube-z[0][None]+c0[None]   # continuum put back, no baseline removal
                    m=(v>=W[0])&(v<=W[1]); img=sub[m].sum(0)*c.dv; img=smooth(img,c,tb)
                    mm=np.hypot(c.X-pk[0],c.Y-pk[1])<0.15; Fl[ln]=img[mm].mean()
                P(f' {reg} core(r<0.15") continuum treatment {mode:5s}: HCN {Fl["HCN"]:.3f} HCO+ {Fl["HCO+"]:.3f} Jy/b km/s  R={Fl["HCN"]/Fl["HCO+"]:.3f}')
        # window continuum/line magnitude at the N peak
        pk=peak_near(O,mO,*REG['N']); mm=np.hypot(O.X-pk[0],O.Y-pk[1])<0.15
        P(f' N core continuum {Sc[mm].mean():.3f} mJy/b x 700 km/s = {Sc[mm].mean()*0.7:.2f} Jy/b km/s-equivalent vs HCO+ line {mO[mm].mean():.2f}')
        for reg2,(xx,yy) in REG.items():
            pk2=peak_near(O,mO,xx,yy); m2=np.hypot(O.X-pk2[0],O.Y-pk2[1])<0.15
            P(f' {reg2} core residual baseline level in the (continuum-subtracted) cubes: HCN {smooth(cH,H,tb)[m2].mean()*1e3:+.4f}, HCO+ {smooth(cO,H,tb)[m2].mean()*1e3:+.4f} mJy/b vs continuum {Sc[m2].mean():.3f} mJy/b -> leftover <= {100*abs(smooth(cO,H,tb)[m2].mean()*1e3)/max(Sc[m2].mean(),1e-9):.0f}% of continuum')
        # HCN v2=1f in the HCO+ cube at +401 km/s
        for reg,(x0,y0) in REG.items():
            pk=peak_near(O,mO,x0,y0); mm=np.hypot(O.X-pk[0],O.Y-pk[1])<0.15
            for lo,hi in ((351,451),(300,500)):
                w=(O.v>=lo)&(O.v<=hi); im=smooth(O.cube[w].sum(0)*O.dv,O,tb); n_=smooth(null_m0(O,(lo,hi)),O,tb)
                P(f' {reg} HCO+ cube {lo}..{hi} km/s (HCN v2=1f 3-2 at +401): I={im[mm].mean()*1e3:.1f} mJy/b km/s, null rms {mad(n_)*1e3:.1f}; /HCN(v=0) core = {im[mm].mean()/(mH[mm].mean()):.3f}')
open('/workspace/out/q10/q10_corr.log','w').write('\n'.join(log)); json.dump({'|'.join(map(str,k)):v for k,v in res.items()},open('/workspace/out/q10/q10_corr.json','w'),indent=1)
