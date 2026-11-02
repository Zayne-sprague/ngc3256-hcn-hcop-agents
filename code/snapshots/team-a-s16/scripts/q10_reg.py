# Q10 v2: pooled N+S regression of R=HCN/HCO+ on log continuum, log HCO+, radius; shifted-continuum null; equal-brightness N/S comparison;
# individual HCN & HCO+ radial profiles vs continuum; velocity-resolved R at the cores; uv-recovery check of the continuum.
# python3 -I q10_reg.py -> out/q10v2/reg.log
import sys,json,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
from astropy.io import fits
from q7_lib import *
rng=np.random.default_rng(11)
log=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); log.append(s)
def mad(a): a=a[np.isfinite(a)]; return 1.4826*np.median(np.abs(a-np.median(a)))
BASEW=(-350,350); REG={'N':(0.0,2.63),'S':(-0.15,-2.64)}
NULLW=(-760,-540)
def null_m0(c,w,wn=NULLW):
    m=np.zeros(len(c.v),bool)
    for a,b in BASE: m|=(c.v>=a)&(c.v<=b)
    m&=~((c.v>=wn[0]-10)&(c.v<=wn[1]+10)); raw=c.cube
    X=np.vstack([np.ones(m.sum()),c.v[m]/500.]).T; Pm=np.linalg.pinv(X); cc=np.tensordot(Pm,raw[m],axes=(1,0))
    k=(c.v>=wn[0])&(c.v<=wn[1]); Xk=np.vstack([np.ones(k.sum()),c.v[k]/500.]).T
    res=raw[k]-np.tensordot(Xk,cc,axes=(1,0)); return res.sum(0)*c.dv*np.sqrt((w[1]-w[0])/(wn[1]-wn[0]))
def cont_map(ous):
    z=np.load(f'/workspace/out/q7/cube_{ous}_HCN.npz'); x0,y0=int(z['x0']),int(z['y0']); ny,nx=z['cube'].shape[1:]
    d=fits.getdata(f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')[0,0]
    return np.nan_to_num(d[y0:y0+ny,x0:x0+nx].astype(np.float64))
def ols(X,y):
    X=np.column_stack([np.ones(len(y))]+X); b=np.linalg.lstsq(X,y,rcond=None)[0]; res=y-X@b; return b,(res**2).sum()
def zs(a): return (a-a.mean())/a.std()
CASES={'Xb4b@0.30':('Xb4b',(0.30,0.30,0.),2.0),'Xb4d@1.05':('Xb4d',(1.05,1.05,0.),3.0)}
for case,(ous,tb,rmax) in CASES.items():
    H,O=Cube(ous,'HCN'),Cube(ous,'HCO+'); c0=cont_map(ous)
    Sc=smooth(c0,H,tb)*1e3; mH=smooth(H.m0(*BASEW),H,tb); mO=smooth(O.m0(*BASEW),O,tb)
    sH,sO=mad(smooth(null_m0(H,BASEW),H,tb)),mad(smooth(null_m0(O,BASEW),O,tb)); sC=mad(Sc[np.hypot(H.X,H.Y)>4.5])
    pk={k:peak_near(O,mO,*c) for k,c in REG.items()}
    rN=np.hypot(O.X-pk['N'][0],O.Y-pk['N'][1]); rS=np.hypot(O.X-pk['S'][0],O.Y-pk['S'][1])
    near=np.where(rN<rS,0,1); r=np.minimum(rN,rS); step=max(1,int(round(tb[0]/H.pix)))
    P(f'\n=== {case}: beam {tb[0]}", noise HCN {sH*1e3:.1f} HCO+ {sO*1e3:.1f} mJy/b km/s, cont {sC:.3f} mJy/b; beam-spaced step {step} px; r<{rmax}"')
    def pool(Scm,off):
        g=np.zeros(Scm.shape,bool); g[off[1]::step,off[0]::step]=True
        s=g&(r<rmax)&(mH>3*sH)&(mO>3*sO)&(Scm>3*sC)
        return s
    def fit_all(Scm,s):
        y=(mH/mO)[s]; ls=np.log10(Scm[s]); lo=np.log10(mO[s]*1e3); rr=r[s]; ind=near[s].astype(float)
        out={}
        for name,Xs in (('logSc',[zs(ls)]),('logHCO+',[zs(lo)]),('r',[zs(rr)]),('logSc+logHCO++r',[zs(ls),zs(lo),zs(rr)]),('logSc+logHCO++r+Snuc',[zs(ls),zs(lo),zs(rr),ind-ind.mean()]),('Snuc only',[ind-ind.mean()])):
            b,rss=ols(Xs,y); out[name]=(b,1-rss/((y-y.mean())**2).sum())
        return out,len(y)
    # observed over all grid offsets
    allfits={}; ns=[]
    for ox in range(step):
        for oy in range(step):
            s=pool(Sc,(ox,oy))
            if s.sum()<12: continue
            f,n=fit_all(Sc,s); ns.append(n)
            for k,(b,r2) in f.items(): allfits.setdefault(k,[]).append((b,r2))
    P(f' pooled pixels per grid offset n={np.mean(ns):.0f} (N {int((pool(Sc,(0,0))&(near==0)).sum())}, S {int((pool(Sc,(0,0))&(near==1)).sum())} at offset 0), {len(ns)} offsets')
    s0=pool(Sc,(0,0)); y0=(mH/mO)[s0]
    for k,v in allfits.items():
        R2=np.mean([x[1] for x in v]); bs=np.mean([x[0] for x in v],axis=0)
        # block bootstrap on offset (0,0) grid
        idx=np.where(s0)[0] if False else None
        P(f'  model R ~ {k:24s}: <R^2>={R2:.2f}, standardised coeffs (const, ...) = {np.round(bs,3).tolist()}')
    # bootstrap of the 4-term model on offset 0
    ii=np.argwhere(s0); bb=[]
    for _ in range(1000):
        j=rng.integers(0,len(ii),len(ii)); sel=ii[j]
        y=(mH/mO)[sel[:,0],sel[:,1]]; ls=np.log10(Sc[sel[:,0],sel[:,1]]); lo=np.log10(mO[sel[:,0],sel[:,1]]*1e3); rr=r[sel[:,0],sel[:,1]]; ind=near[sel[:,0],sel[:,1]].astype(float)
        b,_=ols([zs(ls),zs(lo),zs(rr),ind-ind.mean()],y); bb.append(b)
    bb=np.array(bb); nm=['const','logSc','logHCO+','r','Snuc']
    P('  bootstrap (offset 0, 1000) of R ~ logSc+logHCO++r+Snuc [std. coeff 16/50/84]: '+'; '.join(f'{n} {np.percentile(bb[:,i],16):+.3f}/{np.percentile(bb[:,i],50):+.3f}/{np.percentile(bb[:,i],84):+.3f}' for i,n in enumerate(nm)))
    # null: shift the continuum map by 3-8 beams (random direction), refit logSc coeff in the 3-term and 4-term model
    nb=[];nb4=[]
    for _ in range(300):
        sh=int(rng.uniform(3,8)*tb[0]/H.pix); th=rng.uniform(0,2*np.pi); dy,dx=int(sh*np.sin(th)),int(sh*np.cos(th))
        Scs=np.roll(np.roll(Sc,dy,0),dx,1); s=pool(Scs,(0,0))
        if s.sum()<12: continue
        f,_=fit_all(Scs,s); nb.append(f['logSc+logHCO++r'][0][1]); nb4.append(f['logSc+logHCO++r+Snuc'][0][1])
    P(f'  NULL (continuum map shifted by 3-8 beams, n={len(nb)}): std. coeff of logSc in 3-term model {np.mean(nb):+.3f}+-{np.std(nb):.3f}; in 4-term model (with S-nucleus indicator) {np.mean(nb4):+.3f}+-{np.std(nb4):.3f}; observed 3-term {np.mean([x[0][1] for x in allfits["logSc+logHCO++r"]]):+.3f}, 4-term {np.mean([x[0][1] for x in allfits["logSc+logHCO++r+Snuc"]]):+.3f}')
    # per-nucleus regressions (standardised slopes, bootstrap on offset-0 grid, shifted-continuum null)
    for nm_,k in (('N',0),('S',1)):
        sk=s0&(near==k); ii=np.argwhere(sk)
        def f2(sel,Scm):
            y=(mH/mO)[sel[:,0],sel[:,1]]; ls=np.log10(Scm[sel[:,0],sel[:,1]]); lo=np.log10(mO[sel[:,0],sel[:,1]]*1e3); rr=r[sel[:,0],sel[:,1]]
            return ols([zs(ls)],y)[0][1],ols([zs(lo)],y)[0][1],ols([zs(rr)],y)[0][1],ols([zs(ls),zs(lo),zs(rr)],y)[0][1]
        obs=f2(ii,Sc); bs=np.array([f2(ii[rng.integers(0,len(ii),len(ii))],Sc) for _ in range(500)])
        nl=[]
        for _ in range(200):
            sh=int(rng.uniform(3,8)*tb[0]/H.pix); th=rng.uniform(0,2*np.pi); Scs=np.roll(np.roll(Sc,int(sh*np.sin(th)),0),int(sh*np.cos(th)),1)
            sel=np.argwhere(sk&(Scs>0)); 
            if len(sel)>8: nl.append(f2(sel,Scs)[0])
        P(f'  {nm_}-only (n={len(ii)}): std. slope of R on logSc {obs[0]:+.3f} [{np.percentile(bs[:,0],16):+.3f},{np.percentile(bs[:,0],84):+.3f}], on logHCO+ {obs[1]:+.3f}, on r {obs[2]:+.3f}; partial logSc coeff (with logHCO+, r) {obs[3]:+.3f} [{np.percentile(bs[:,3],16):+.3f},{np.percentile(bs[:,3],84):+.3f}]; shifted-continuum null for the logSc slope {np.mean(nl):+.3f}+-{np.std(nl):.3f}')
    # equal-brightness comparison: bins of HCO+ and of Sc
    P('  equal-brightness comparison (all pixels r<%.1f" with HCN,HCO+ >3sigma, no continuum cut; mean R +-sem; beam-spaced offset-0 grid):'%rmax)
    g=np.zeros(Sc.shape,bool); g[0::step,0::step]=True; sx=g&(r<rmax)&(mH>3*sH)&(mO>3*sO)
    for lab,arr,edges in (('HCO+ Jy/b km/s',mO,np.array([0.1,0.2,0.4,0.8,1.6]) if tb[0]<0.5 else np.array([0.3,0.6,1.2,2.4,5.])),('Sc mJy/b',Sc,np.array([0.05,0.15,0.3,0.6,1.2,2.5]) if tb[0]<0.5 else np.array([0.2,0.5,1.0,2.0,4.0,8.]))):
        for a,b in zip(edges[:-1],edges[1:]):
            row=f'   {lab} {a:.2f}-{b:.2f}:'
            for nm_,k in (('N',0),('S',1)):
                m=sx&(near==k)&(arr>=a)&(arr<b); R_=(mH/mO)[m]
                row+=f'  {nm_} n={m.sum():2d} R={R_.mean() if m.sum() else np.nan:.2f}+-{R_.std(ddof=1)/np.sqrt(m.sum()) if m.sum()>1 else np.nan:.2f};'
            P(row)
    if case=='Xb4b@0.30':
        # individual HCN and HCO+ profiles normalised to the core, vs continuum
        P('  radial profiles normalised to the core value (r<0.15"): HCN, HCO+, continuum, R')
        edges=[0,0.15,0.25,0.5,1,1.5,2]
        for nm_,k in (('N',0),('S',1)):
            prof=[]
            for a,b in zip(edges[:-1],edges[1:]):
                m=(near==k)&(r>=a)&(r<b); prof.append((mH[m].mean(),mO[m].mean(),Sc[m].mean()))
            prof=np.array(prof); P(f'   {nm_}: r-edges {edges}\n     HCN/core {np.round(prof[:,0]/prof[0,0],3).tolist()}\n     HCO+/core {np.round(prof[:,1]/prof[0,1],3).tolist()}\n     cont/core {np.round(prof[:,2]/prof[0,2],3).tolist()}\n     abs core: HCN {prof[0,0]:.3f}, HCO+ {prof[0,1]:.3f} Jy/b km/s, cont {prof[0,2]:.3f} mJy/b; HCN/cont {np.round(prof[:,0]/prof[:,2],2).tolist()} HCO+/cont {np.round(prof[:,1]/prof[:,2],2).tolist()}')
        # velocity-resolved R at the cores (aperture r<0.25" at 0.30" beam), 44-km/s bins, blue/red wings
        P('  velocity-resolved R=HCN/HCO+ in r<0.25" of each nucleus (Xb4b smoothed to 0.30"), bins of 8 channels (~35 km/s), v rel. VSYS:')
        vs=np.arange(-300,301,40)
        for nm_,k in (('N',0),('S',1)):
            m=(near==k)&(r<0.25); out=[]
            Hs=np.array([smooth(H.cube[i],H,tb)[m].mean() for i in range(len(H.v))]); Os=np.array([smooth(O.cube[i],O,tb)[m].mean() for i in range(len(O.v))])
            # noise per bin from line-free channels
            for lo in vs[:-1]:
                w=(H.v>=lo)&(H.v<lo+40); wo=(O.v>=lo)&(O.v<lo+40); h_=Hs[w].mean(); o_=Os[wo].mean()
                out.append((lo+20,h_*1e3,o_*1e3))
            bn=(np.abs(H.v)>520)&(np.abs(H.v)<900); bno=(np.abs(O.v)>520)&(np.abs(O.v)<900); eh=np.std(Hs[bn])/np.sqrt(9); eo=np.std(Os[bno])/np.sqrt(9)
            P(f'      peak velocity (rel. VSYS): HCN {H.v[np.argmax(Hs)]:+.0f}, HCO+ {O.v[np.argmax(Os)]:+.0f} km/s; channel width HCN {H.dv:.1f} HCO+ {O.dv:.1f} km/s')
            P(f'   {nm_}: bin noise HCN {eh*1e3:.2f}, HCO+ {eo*1e3:.2f} mJy/b; '+' '.join(f'{c:+.0f}:{h:.1f}/{o:.1f}={h/o:.2f}' for c,h,o in out if o>3*eo*1e3))
            for lab,(a,b) in (('blue 40-200',(-200,-40)),('core +-40',(-40,40)),('red 40-200',(40,200)),('blue 200-350',(-350,-200)),('red 200-350',(200,350))):
                w=(H.v>=a)&(H.v<b); wo=(O.v>=a)&(O.v<b); h_=Hs[w].mean(); o_=Os[wo].mean(); e_h=eh*np.sqrt(9/w.sum()); e_o=eo*np.sqrt(9/wo.sum())
                Rv=h_/o_; eR=Rv*np.hypot(e_h/h_,e_o/o_)
                P(f'      {lab:13s}: HCN {h_*1e3:.2f}+-{e_h*1e3:.2f} HCO+ {o_*1e3:.2f}+-{e_o*1e3:.2f} mJy/b  R={Rv:.2f}+-{eR:.2f}')
# ---- uv recovery check: continuum in Xb4d vs Xb4f (compact array, 5.8x3.5") smoothed to the common beam; and Band 3 MRS from header
from scipy.signal import fftconvolve
P('\n=== uv-coverage check of the 98/258 GHz continuum (image-plane): flux in r<3" of the field centre, all smoothed to 6x6"')
def gk(sig):
    n=int(4*sig)+1; y,x=np.mgrid[-n:n+1,-n:n+1]; k=np.exp(-0.5*(x**2+y**2)/sig**2); return k/k.sum()
def aperflux(path,rad,beam_in=None,tb_=6.0,center=(RA0,DE0)):
    from astropy.wcs import WCS
    h=fits.getheader(path); d=np.nan_to_num(fits.getdata(path).squeeze().astype(float)); w=WCS(h).celestial
    pixs=abs(h['CDELT1'])*3600; bmaj,bmin=h['BMAJ']*3600,h['BMIN']*3600
    sig=np.sqrt(max(tb_**2-bmaj*bmin,0))/2.3548/pixs
    sm=fftconvolve(d,gk(sig),mode='same')*(tb_**2)/(bmaj*bmin)   # Jy/beam(6")
    x,y=w.world_to_pixel_values(*center); yy,xx=np.mgrid[0:d.shape[0],0:d.shape[1]]
    ap=np.hypot(xx-x,yy-y)*pixs<=rad
    return d[ap].sum()/(1.1331*bmaj*bmin/pixs**2), sm[int(round(float(y))),int(round(float(x)))], bmaj,bmin
for lab,p in (('Xb4b 0.24" (258)','/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'),('Xb4d 0.9" (258)','/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'),('Xb4f 5.8x3.5" (258)','/data/member.uid___A001_X1273_Xb4f.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits'),('2016.1.00965 Band3 (98)','/archive/files/2016.1.00965.S/member.uid___A001_X879_X94a.NGC_3256_sci.spw29.mfs.I.pbcor.fits')):
    for rad in (3.,):
        F,c6,bj,bn_=aperflux(p,rad)
        P(f' {lab:26s} beam {bj:.2f}x{bn_:.2f}: total flux r<{rad}" = {F*1e3:.2f} mJy; centre value after smoothing to 6" = {c6*1e3:.2f} mJy/6"-beam')
open('/workspace/out/q10v2/reg.log','w').write('\n'.join(log))
