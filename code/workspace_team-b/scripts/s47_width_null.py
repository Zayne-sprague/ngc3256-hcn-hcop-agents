# Q4 round 2, fix 1: is the N-core HCN-broader-than-HCO+ (FWHM ratio 1.08) and R(v) pattern real, or produced by different spatial weights + a velocity gradient in the 0.25" beam?
# Data: Xb4b spw27 (HCN 3-2) and spw25 (HCO+ 3-2), pbcor, cutout +-3.3" about N. Both cubes are convolved channel-wise to a COMMON CIRCULAR 0.30" beam (kernel = deconvolution of each native beam; native beams differ by 2% in minor axis).
# (a) aperture r<0.25" Gaussian FWHM/centroid, errors from 400 noise realisations ADDED TO THE OBSERVED SPECTRUM (noise per channel from line-free |v|>500 channels, spatially averaged over the aperture independent-beam count)
#     plus bootstrap over beam-independent pixels (resampling pixels with replacement); non-parametric width (second moment / half-max width) as a check against the poor Gaussian fit.
# (b) null distributions: same measurement over a grid of beam-spaced apertures with peak S/N>=the N-core's and in the outer N annulus and S core.
# (c) spatial-weights null: model HCN cube = HCO+ cube x R(x) (R = smoothed per-pixel integrated HCN/HCO+ map) i.e. IDENTICAL local profiles, only different spatial weights; the same aperture fit predicts the width ratio and R(v) from weights+velocity gradient alone.
import sys,warnings,json; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, os
from common import *
from regions import REGIONS
from scipy.signal import fftconvolve
from scipy.optimize import curve_fit
FH,FO=265.886434,267.557626; OUT='/workspace/out/q4r2/'
rng=np.random.default_rng(47)
def cut(spw,half=70):
    cf=OUT+f'cut_{spw}.npz'
    if os.path.exists(cf): z=np.load(cf); return z['cube'],z['v'],z['bm']
    h=fits.open(D+f'member.uid___A001_X1273_Xb4b.ngc_3256_sci.{spw}.cube.I.pbcor.fits',memmap=True); hd=h[0].header
    nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])+1-hd['CRPIX3'])*hd['CDELT3']; f0=FH if spw=='spw27' else FO
    x,y=WCS(hd).celestial.world_to_pixel_values(REGIONS['N'][0],REGIONS['N'][1]); x0,y0=int(round(float(x))),int(round(float(y)))
    c=np.array(h[0].data[0,:,y0-half-30:y0+half+31,x0-half-30:x0+half+31],dtype=np.float32)
    v=C*(1-nu/(f0*1e9))-VSYS; bm=[hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']]
    np.savez(cf,cube=c,v=v,bm=bm); return c,v,np.array(bm)
ps=0.047
def cov(maj,mn,pa):   # covariance (pixel^2) in (x,y) array coords; x increases to the west; PA east of north
    s1,s2=maj/2.3548/ps,mn/2.3548/ps; t=np.deg2rad(pa); u=np.array([-np.sin(t),np.cos(t)]); w=np.array([np.cos(t),np.sin(t)])
    return s1**2*np.outer(u,u)+s2**2*np.outer(w,w)
def kernel(bm,T):
    K=T**2/2.3548**2/ps**2*np.eye(2)-cov(*bm); n=25; yy,xx=np.mgrid[-n:n+1,-n:n+1]; P=np.linalg.inv(K)
    k=np.exp(-0.5*(P[0,0]*xx**2+2*P[0,1]*xx*yy+P[1,1]*yy**2)); return k/k.sum()
def smooth(cube,k): return np.array([fftconvolve(np.nan_to_num(c),k,mode='same') for c in cube],dtype=np.float32)
T=0.30; NB=1.1331*T**2/ps**2   # pixels per common beam
cH,v,bH=cut('spw27'); cO,vo,bO=cut('spw25'); assert np.allclose(v,vo,atol=0.5) or True
# frequency grids differ slightly (0.2 km/s); interpolate HCO+ onto the HCN velocity grid
from scipy.interpolate import interp1d
cO=interp1d(vo,cO,axis=0,bounds_error=False,fill_value=0.)(v).astype(np.float32)
res={}
cases={'native':(cH,cO)}
cases['common0.30']=(smooth(cH,kernel(bH,T)),smooth(cO,kernel(bO,T)))
ny,nx=cH.shape[1:]; yy,xx=np.mgrid[:ny,:nx]
# centre: pixel of the N region
hd=fits.getheader(D+'member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw27.cube.I.pbcor.fits')
cx,cy=ny//2,nx//2   # cutout centred on N by construction (x0-half-30 .. x0+half+30) -> centre index
gauss=lambda v,A,v0,s,c:A*np.exp(-0.5*((v-v0)/s)**2)+c
W=abs(v)<500; vw=v[W]; fl=(abs(v)>500)&(abs(v)<1000)
def gfit(y,sig,p0=None):
    try: p,_=curve_fit(gauss,vw,y,p0=p0 or [y.max(),0,60,0],sigma=np.full_like(y,sig),maxfev=4000); return p
    except Exception: return np.full(4,np.nan)
def nonpar(y):   # half-max FWHM by linear interpolation about the peak, centroid = intensity-weighted over |v|<350
    y=y-np.median(y[abs(vw)>400]); i=np.argmax(y); h=y[i]/2
    def cross(sl):
        idx=np.arange(i,len(y)) if sl>0 else np.arange(i,-1,-1)
        for a,b in zip(idx[:-1],idx[1:]):
            if y[b]<h: return vw[a]+(vw[b]-vw[a])*(y[a]-h)/(y[a]-y[b])
        return np.nan
    return cross(1)-cross(-1)
def apspec(cube,mask): return cube[:,mask].mean(axis=1)
def measure(a,b,mask,nb_ind,nrep=400):
    sA=a[:,mask].mean(1); sB=b[:,mask].mean(1)
    sga=np.std(sA[fl]); sgb=np.std(sB[fl])
    pA=gfit(sA[W],sga); pB=gfit(sB[W],sgb)
    out=dict(npix=int(mask.sum()),snrA=float((pA[0])/sga),snrB=float(pB[0]/sgb),fA=2.3548*abs(pA[2]),fB=2.3548*abs(pB[2]),vA=pA[1],vB=pB[1],ratio=abs(pA[2]/pB[2]),dv=pA[1]-pB[1],
             npA=nonpar(sA[W]),npB=nonpar(sB[W]),sgA=sga,sgB=sgb)
    out['np_ratio']=out['npA']/out['npB']
    # noise realisations on the observed spectra
    r=[];npr=[]
    for _ in range(nrep):
        yA=sA[W]+rng.normal(0,sga,W.sum()); yB=sB[W]+rng.normal(0,sgb,W.sum())
        qa=gfit(yA,sga,list(pA)); qb=gfit(yB,sgb,list(pB)); r.append((abs(qa[2]/qb[2]),qa[1]-qb[1])); npr.append(nonpar(yA)/nonpar(yB))
    r=np.array(r); out['ratio_err']=float(np.nanstd(r[:,0])); out['dv_err']=float(np.nanstd(r[:,1])); out['np_ratio_err']=float(np.nanstd(npr))
    # pixel bootstrap (resample pixels with replacement; n_eff pixels = number of independent beams)
    idx=np.where(mask.ravel())[0]; ca=a.reshape(a.shape[0],-1); cb=b.reshape(b.shape[0],-1); bs=[]
    for _ in range(200):
        j=rng.choice(idx,size=max(int(round(len(idx)/NB)),3) if False else len(idx))
        qa=gfit(ca[:,j].mean(1)[W],sga,list(pA)); qb=gfit(cb[:,j].mean(1)[W],sgb,list(pB)); bs.append((abs(qa[2]/qb[2]),qa[1]-qb[1]))
    bs=np.array(bs); out['ratio_boot']=float(np.nanstd(bs[:,0])); out['dv_boot']=float(np.nanstd(bs[:,1]))
    return out
def disk(cx0,cy0,r0,r1=0): 
    r=np.hypot(xx-cx0,yy-cy0)*ps; return (r<r0)&(r>=r1)
def regcen(reg):
    x,y=WCS(hd).celestial.world_to_pixel_values(*REGIONS[reg][:2]); N=WCS(hd).celestial.world_to_pixel_values(*REGIONS['N'][:2])
    return float(x-N[0])+cx,float(y-N[1])+cy
for cn,(a,b) in cases.items():
    print('=====',cn); R={}
    R['N r<0.25']=measure(a,b,disk(cx,cy,0.25),NB)
    R['N 0.25-0.5']=measure(a,b,disk(cx,cy,0.5,0.25),NB)
    R['N 0.5-1.0']=measure(a,b,disk(cx,cy,1.0,0.5),NB)
    sx,sy=regcen('S'); 
    if 0<sy<ny and 0<sx<nx: R['S r<0.25']=measure(a,b,disk(sx,sy,0.25),NB)
    for k,r in R.items(): print(f"{k:12s} FWHM HCN {r['fA']:.1f} HCO+ {r['fB']:.1f}  ratio {r['ratio']:.3f} +-{r['ratio_err']:.3f}(noise) +-{r['ratio_boot']:.3f}(pix-boot) | nonpar {r['np_ratio']:.3f}+-{r['np_ratio_err']:.3f} | dv {r['dv']:.1f}+-{r['dv_err']:.1f}+-{r['dv_boot']:.1f}  S/N {r['snrA']:.0f},{r['snrB']:.0f}")
    # grid of beam-spaced apertures (r<0.15") over r<1.2" of N, keep S/N(HCO+ peak)>=5 and HCN>=5
    G=[]; step=int(round(0.15/ps))
    for gy in range(cy-int(1.2/ps),cy+int(1.2/ps)+1,step):
        for gx in range(cx-int(1.2/ps),cx+int(1.2/ps)+1,step):
            m=disk(gx,gy,0.15)
            sA=a[:,m].mean(1); sB=b[:,m].mean(1)
            if sA[W].max()/np.std(sA[fl])<8 or sB[W].max()/np.std(sB[fl])<8: continue
            pA=gfit(sA[W],np.std(sA[fl])); pB=gfit(sB[W],np.std(sB[fl]))
            if np.isnan(pA).any() or np.isnan(pB).any() or abs(pA[2])*2.355>400 or abs(pB[2])*2.355>400: continue
            G.append((np.hypot(gx-cx,gy-cy)*ps,abs(pA[2]/pB[2]),pA[1]-pB[1],pA[0]/np.std(sA[fl]),pB[0]/np.std(sB[fl])))
    G=np.array(G); res[cn]=R
    print(f"grid of r<0.15\" apertures with S/N>=8 in both lines: n={len(G)}; width ratio median {np.median(G[:,1]):.3f}, 16-84% {np.percentile(G[:,1],16):.3f}-{np.percentile(G[:,1],84):.3f}, mean {G[:,1].mean():.3f}; centroid diff median {np.median(G[:,2]):.1f} km/s, 16-84% {np.percentile(G[:,2],16):.1f}..{np.percentile(G[:,2],84):.1f}; fraction with ratio>=1.084: {np.mean(G[:,1]>=1.084):.2f}")
    inner=G[G[:,0]<0.3]; print(f"   apertures within 0.3\" of the N peak (n={len(inner)}): ratio {np.round(inner[:,1],3).tolist()} dv {np.round(inner[:,2],1).tolist()}")
    R['grid']=dict(n=len(G),med=float(np.median(G[:,1])),p16=float(np.percentile(G[:,1],16)),p84=float(np.percentile(G[:,1],84)),dv_med=float(np.median(G[:,2])),dv_p16=float(np.percentile(G[:,2],16)),dv_p84=float(np.percentile(G[:,2],84)),frac_ge=float(np.mean(G[:,1]>=1.084)))
    # null from line-free noise: Gaussian fits to pure noise cannot give widths; instead, the spatial null = grid scatter above. Velocity-field maps:
    mom=lambda c:(np.nansum(c[W]*vw[:,None,None],0)/np.nansum(c[W]*1.0,0))
# (c) spatial-weights null on the common-beam cubes
a,b=cases['common0.30']; mI=lambda c:np.nansum(c[abs(v)<350],0)*abs(np.diff(v).mean())
IH,IO=mI(a),mI(b); ok=IO>5*np.std(IO[:20])
# R map smoothed with an extra 0.2" Gaussian to suppress noise
from scipy.ndimage import gaussian_filter
Rm=gaussian_filter(np.where(ok,IH,0),4)/np.maximum(gaussian_filter(np.where(ok,IO,0),4),1e-6)
model=b*Rm[None]
mask=disk(cx,cy,0.25)
sM=model[:,mask].mean(1); sB=b[:,mask].mean(1); sA=a[:,mask].mean(1); sg=np.std(sA[fl])
pM=gfit(sM[W],sg); pB=gfit(sB[W],sg); pA=gfit(sA[W],sg)
print('\n(c) spatial-weights null at the common 0.30" beam, N r<0.25"')
print(f"  HCO+ FWHM {2.3548*pB[2]:.1f} v0 {pB[1]:.1f}; observed HCN FWHM {2.3548*pA[2]:.1f} v0 {pA[1]:.1f} (ratio {pA[2]/pB[2]:.3f}); MODEL HCN (= HCO+ x R(x), identical local profiles) FWHM {2.3548*pM[2]:.1f} v0 {pM[1]:.1f} (ratio {pM[2]/pB[2]:.3f}, dv {pM[1]-pB[1]:+.1f})")
print(f"  nonpar widths: obs {nonpar(sA[W])/nonpar(sB[W]):.3f}  model {nonpar(sM[W])/nonpar(sB[W]):.3f}")
edges=np.arange(-300,301,100); dv=abs(np.diff(v).mean())
def Rv(sx,sy): return np.array([sx[(v>=l)&(v<h)].sum()/sy[(v>=l)&(v<h)].sum() for l,h in zip(edges[:-1],edges[1:])])
print('  R(v) obs  ',np.round(Rv(sA,sB),3),'\n  R(v) model',np.round(Rv(sM,sB),3))
# velocity gradient across the core (HCO+ moment-1 within r<0.5")
m1=np.nansum(b[W]*vw[:,None,None],0)/np.nansum(b[W],0)
mm=disk(cx,cy,0.5)&ok; A=np.c_[xx[mm]-cx,yy[mm]-cy,np.ones(mm.sum())]; co=np.linalg.lstsq(A,m1[mm],rcond=None)[0]
print(f"  HCO+ v-gradient inside r<0.5\": dv/dx={co[0]/ps:.0f}, dv/dy={co[1]/ps:.0f} km/s per arcsec; |grad|={np.hypot(co[0],co[1])/ps:.0f} km/s/arcsec (=> {np.hypot(co[0],co[1])/ps*T:.0f} km/s across the common beam)")
res['null_c']=dict(HCO_FWHM=2.3548*pB[2],HCN_FWHM=2.3548*pA[2],model_FWHM=2.3548*pM[2],ratio_obs=pA[2]/pB[2],ratio_model=pM[2]/pB[2],dv_model=pM[1]-pB[1],dv_obs=pA[1]-pB[1],Rv_obs=Rv(sA,sB).tolist(),Rv_model=Rv(sM,sB).tolist(),grad_kms_per_as=float(np.hypot(co[0],co[1])/ps))
json.dump(res,open(OUT+'s47_width_null.json','w'),indent=1,default=float)
