# Q10 part 2: H40alpha (2016.1.00965.S, 0.65x0.56" cube) + 3.1 mm continuum vs HCN/HCO+ 3-2 at matched 1.05" (Xb4d), N and S.
# python3 -I q10_h40a.py ; outputs out/q10/q10_h40a.log/json
import sys,json,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.stats import spearmanr
from q7_lib import Cube,smooth,peak_near,RA0,DE0
D='/archive/files/2016.1.00965.S/member.uid___A001_X879_X94a.NGC_3256_sci.spw29.'
hc=fits.getheader(D+'cube.I.pbcor.fits'); cub=np.nan_to_num(fits.getdata(D+'cube.I.pbcor.fits')[0].astype(np.float64))
hm=fits.getheader(D+'mfs.I.pbcor.fits'); mfs=np.nan_to_num(fits.getdata(D+'mfs.I.pbcor.fits')[0,0].astype(np.float64))
nu=hc['CRVAL3']+hc['CDELT3']*np.arange(cub.shape[0]); dnu=abs(hc['CDELT3'])
REST=99.02295e9; VSYS=2780.
vopt=lambda f: 299792.458*(REST/f-1)-0   # optical-def velocity of H40a
vel=vopt(nu)  # km/s (optical, ~LSRK)
dv=abs(np.median(np.diff(vel)))
bm=(hc['BMAJ']*3600,hc['BMIN']*3600,hc['BPA']); pix=abs(hc['CDELT1'])*3600
log=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); log.append(s)
P(f'H40a cube: beam {bm[0]:.2f}x{bm[1]:.2f} PA {bm[2]:.0f}, pix {pix:.3f}", dv={dv:.1f} km/s, v range {vel.min():.0f}..{vel.max():.0f}; mfs beam {hm["BMAJ"]*3600:.2f}x{hm["BMIN"]*3600:.2f} at {hm["CRVAL3"]/1e9:.2f} GHz')
w=WCS(hc).celestial
def mad(a): a=a[np.isfinite(a)]; return 1.4826*np.median(np.abs(a-np.median(a)))
# H40a moment 0, with a linear baseline from channels outside +-1500 km/s of VSYS (excluding edges contaminated by H40a only)
v=vel-VSYS; base=(np.abs(v)>700)
X=np.vstack([np.ones(base.sum()),v[base]/500.]).T; Pm=np.linalg.pinv(X)
cc=np.tensordot(Pm,cub[base],axes=(1,0)); sub=cub-np.tensordot(np.vstack([np.ones(len(v)),v/500.]).T,cc,axes=(1,0))
def m0(lo,hi): k=(v>=lo)&(v<=hi); return sub[k].sum(0)*dv, k.sum()
# line-free windows for noise: same width windows elsewhere
Wn=(-350,350)
m_line,nch=m0(*Wn)
nulls=[]
for c0 in (-1300,-1000,1000,1300,1600,1900):
    lo,hi=c0-350,c0+350
    if (v>=lo).any() and v[v<=hi].size and v.min()<=lo and v.max()>=hi: nulls.append(m0(lo,hi)[0])
P(f'm0 window {Wn} km/s ({nch} ch); null windows: {len(nulls)}')
# sky -> pixel map to Xb4d grid
H=Cube('Xb4d','HCN'); O=Cube('Xb4d','HCO+'); tb=(1.05,1.05,0.)
# smooth native 0.65x0.56 (beam PA 66) to 1.05 circular via Gaussian convolution (deconvolved kernel, elliptical in sky approx: circularise by sigma per axis in the beam frame)
def smooth_to(img,bmaj,bmin,bpa,target=1.05):
    # kernel = target beam deconvolved from native; native beam ellipse rotated by PA (E of N); build in pixel grid with x=RA
    sa=np.sqrt(max(target**2-bmaj**2,0))/2.3548/pix; sb=np.sqrt(max(target**2-bmin**2,0))/2.3548/pix
    n=int(5*max(sa,sb)); yy,xx=np.mgrid[-n:n+1,-n:n+1].astype(float)
    th=np.deg2rad(bpa)  # PA from north to east: major axis direction (sin th, cos th) in (x=E.., y=N): x axis increases W here (cdelt<0)
    u=-xx*np.sin(th)+yy*np.cos(th)   # coordinate along major axis (x increases to west => east=-x)
    t=xx*np.cos(th)+yy*np.sin(th)
    k=np.exp(-0.5*((u/sa)**2+(t/sb)**2)); k/=k.sum()
    from scipy.signal import fftconvolve
    # result in Jy/beam(target): (native Jy/beam km/s) * Omega_t/Omega_n
    return fftconvolve(img,k,mode='same')*(target**2)/(bmaj*bmin)
def to_grid(img,c):
    # sample on cube-cache grid c (ra,de arrays)
    RA,DE=np.meshgrid(c.ra,c.de); x,y=w.world_to_pixel_values(RA,DE)
    return map_coordinates(img,[y,x],order=1,mode='nearest')
S_line=smooth_to(m_line,*bm); S_nulls=[smooth_to(n,*bm) for n in nulls]
S_cont=smooth_to(mfs,hm['BMAJ']*3600,hm['BMIN']*3600,hm['BPA'])*1e3/1.0   # mfs is Jy/beam -> smooth_to scales by beam area, keep Jy/beam(1.05)->mJy
# for continuum smooth_to multiplies by target^2/(bmaj*bmin): correct for Jy/beam surface-brightness scaling (done)
L=to_grid(S_line,H)*1e3  # mJy/b km/s
Lc=to_grid(S_cont,H)     # mJy/b (3.1mm cont)
sL=np.mean([mad(to_grid(n,H)*1e3) for n in nulls]) if nulls else np.nan
# noise from field far from sources (beyond 6") in the line map using the grid
R=np.hypot(H.X,H.Y); far=R>6
sL_sp=mad(L[far]); sC_sp=mad(Lc[far])
P(f'noise at 1.05": H40a m0 null-window MAD {sL:.1f} (far-field {sL_sp:.1f}) mJy/b km/s; 3mm cont far-field {sC_sp:.3f} mJy/b')
W=(-350,350)
mH=smooth(H.m0(*W),H,tb)*1e3; mO=smooth(O.m0(*W),O,tb)*1e3   # mJy/b km/s
REG={'N':(0.0,2.63),'S':(-0.15,-2.64)}
res={}
for reg,(x0,y0) in REG.items():
    pk=peak_near(O,mO,x0,y0); r=np.hypot(O.X-pk[0],O.Y-pk[1])
    P(f'\n{reg}: HCO+ peak {pk[0]:+.2f},{pk[1]:+.2f}; H40a peak value (1.05") within 0.6": {L[np.hypot(O.X-pk[0],O.Y-pk[1])<0.6].max():.1f} mJy/b km/s; 3mm cont peak {Lc[np.hypot(O.X-pk[0],O.Y-pk[1])<0.6].max():.3f} mJy/b')
    P(' annulus   R(HCN/HCO+)  <H40a> mJy/b km/s (+-)  <3mm cont> mJy/b  H40a/cont(km/s)')
    rows=[]
    for a,b in ((0,0.5),(0.5,1),(1,1.5),(1.5,2),(2,3)):
        m=(r>=a)&(r<b); Rr=mH[m].sum()/mO[m].sum(); l=L[m].mean(); c=Lc[m].mean()
        # annulus mean noise: independent beams ~ area/beam area
        nb=max(1,m.sum()*H.pix**2/(1.1331*1.05**2)); e=sL_sp/np.sqrt(nb)
        rows.append((a,b,Rr,l,e,c)); P(f'  {a:.1f}-{b:.1f}   {Rr:.3f}   {l:7.1f} +- {e:.1f}   {c:.3f}   {l/c if c>0 else np.nan:.0f}')
    res[reg]=rows
    ra=np.array(rows)
    P(f' annulus-level Spearman(H40a,R)={spearmanr(ra[:,3],ra[:,2])[0]:+.2f}; (cont,R)={spearmanr(ra[:,5],ra[:,2])[0]:+.2f} (n=5, descriptive)')
    # pixel level on beam-spaced grid
    step=7  # 1.05"/0.15
    best=[]
    for rmax in (2.0,3.0):
        for ox in range(step):
            for oy in range(step):
                g=np.zeros(L.shape,bool); g[oy::step,ox::step]=True
                s=g&(r<rmax)&(mH>3*70)&(mO>3*70)
                other=[vv for k,vv in REG.items() if k!=reg][0]; ob=peak_near(O,mO,*other); s&=(r<np.hypot(O.X-ob[0],O.Y-ob[1]))
                if s.sum()>=5: best.append((rmax,s.sum(),spearmanr(L[s],(mH/mO)[s])[0],spearmanr(Lc[s],(mH/mO)[s])[0],spearmanr(mO[s],(mH/mO)[s])[0]))
    best=np.array(best)
    for rmax in (2.0,3.0):
        b=best[best[:,0]==rmax]
        if len(b): P(f' r<{rmax}": pixel-level (beam-spaced, mean over {len(b)} offsets, n~{b[:,1].mean():.0f}): Spearman(H40a,R)={b[:,2].mean():+.2f}; (3mm cont,R)={b[:,3].mean():+.2f}; (HCO+,R)={b[:,4].mean():+.2f}')
    # aperture H40a detection in r<0.5 and r<1 at native 0.65 beam (no smoothing): 
# unsmoothed spectra and S/N of H40a in apertures
P('\nH40a native-resolution spectrum/aperture detection:')
for reg,(dx,dy) in REG.items():
    ra=RA0+dx/3600/np.cos(np.deg2rad(DE0)); de=DE0+dy/3600
    x,y=w.world_to_pixel_values(ra,de); x=float(x); y=float(y)
    yy,xx=np.mgrid[0:cub.shape[1],0:cub.shape[2]]
    for rad in (0.35,0.5,1.0):
        ap=np.hypot(xx-x,yy-y)*pix<=rad
        spec=sub[:,ap].mean(1)*1e3  # mean pbcor Jy/beam per pixel (mJy/beam) -- point at native beam
        k=(v>=Wn[0])&(v<=Wn[1]); I=spec[k].sum()*dv
        # noise: std of same-width integrals in line-free parts via channel rms scaling
        rms=np.std(spec[base]); eI=rms*dv*np.sqrt(k.sum())
        # Gaussian fit-free: line peak (3-channel mean) near v in +-150
        P(f' {reg} aperture r<{rad}": integral over {Wn} = {I:.1f} +- {eI:.1f} mJy/b km/s (S/N {I/eI:.1f}); channel rms {rms:.3f} mJy/b; continuum(mfs) {mfs[ap].mean()*1e3:.3f} mJy/b')
        if rad==0.5:
            sm=np.convolve(spec,np.ones(5)/5,'same')
            P(f'   5-channel (120 km/s) smoothed spectrum, -600..600 km/s step 120: '+' '.join(f'{sm[np.argmin(abs(v-q))]:.2f}' for q in range(-600,601,120)))

# ---- spectral index 98 -> 258 GHz at 1.05" and free-free/dust decomposition; expected H40a for pure free-free
z=np.load('/workspace/out/q7/cube_Xb4d_HCN.npz'); x0,y0=int(z['x0']),int(z['y0']); ny,nx=z['cube'].shape[1:]
d1=np.nan_to_num(fits.getdata('/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')[0,0].astype(np.float64))[y0:y0+ny,x0:x0+nx]
C1=smooth(d1,H,tb)*1e3
nu1,nu3=257.8,98.05
P('\n=== spectral index alpha(98-258 GHz) at matched 1.05" (3mm: 2016.1.00965.S mfs; 1.1mm: Xb4d line-free continuum), N and S')
P(' cont noise: 1.1mm far-field %.3f, 3mm %.3f mJy/b; plus 5%% flux-scale error each'%(mad(C1[far]),sC_sp))
ALPHA={}
for reg,(x0_,y0_) in REG.items():
    pk=peak_near(O,mO,x0_,y0_); r=np.hypot(O.X-pk[0],O.Y-pk[1])
    P(f'{reg}:  annulus  S3mm  S1.1mm  alpha(+-)  f_ff(98) [a_ff=-0.1,a_d=3.5(beta=1.5)]  f_ff(258)  | with a_d=4.0 f_ff(98)   H40a/cont expected (pure ff, Te=1e4) vs measured')
    for a,b in ((0,0.5),(0.5,1),(1,1.5),(1.5,2)):
        m=(r>=a)&(r<b); s3=Lc[m].mean(); s1=C1[m].mean(); nb=max(1,m.sum()*H.pix**2/(1.1331*1.05**2))
        e3=np.hypot(sC_sp/np.sqrt(nb),0.05*s3); e1=np.hypot(mad(C1[far])/np.sqrt(nb),0.05*s1)
        al=np.log(s1/s3)/np.log(nu1/nu3); eal=np.hypot(e1/s1,e3/s3)/np.log(nu1/nu3)
        def fff(al_,ad):
            q=(nu1/nu3); # S1/S3=q^al = f*q^-0.1+(1-f)*q^ad
            return (q**ad-q**al_)/(q**ad-q**-0.1)
        f98=fff(al,3.5); f98_lo=fff(al+eal,3.5); f98_hi=fff(al-eal,3.5); f98b=fff(al,4.0)
        f258=f98*(nu1/nu3)**-0.1/ q_ if False else f98*(nu1/nu3)**-0.1*s3/s1
        l=L[m].mean(); el=sL_sp/np.sqrt(nb)
        exp_ff=3.3*f98   # km/s: pure free-free line/continuum ~3.3 km/s at 98 GHz, Te=1e4 K (Gordon & Sorochenko 2002 scaling), x ff fraction of continuum
        P(f'   {a:.1f}-{b:.1f}  {s3:.3f} {s1:.3f}  alpha={al:.2f}+-{eal:.2f}  f_ff98={f98:.2f} [{f98_hi:.2f},{f98_lo:.2f}]  f_ff258={f258:.2f}  | {f98b:.2f}   {exp_ff:.1f} km/s vs {l/s3:.0f}+-{el/s3:.0f} km/s')
        ALPHA[f'{reg}_{a}-{b}']=dict(s3=float(s3),s1=float(s1),alpha=float(al),ealpha=float(eal),f_ff98=float(f98))
json.dump(ALPHA,open('/workspace/out/q10/q10_alpha.json','w'),indent=1)
open('/workspace/out/q10/q10_h40a.log','w').write('\n'.join(log))
json.dump({k:[list(map(float,r)) for r in v_] for k,v_ in res.items()},open('/workspace/out/q10/q10_h40a.json','w'),indent=1)
open('/workspace/out/q10/q10_h40a.log','w').write('\n'.join(log))
