# Q5 helpers: joint matched-filter fits of blended lines with a shared (HCO+-derived) line shape; fixed-position spectra from s90_q5_extract.py
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import curve_fit
C=299792.458
REFLINE={'Xb47':('spw18',356.73424),'Xb4f':('spw16',267.55763),'Xb4d':('spw25',267.55763),'Xb4b':('spw25',267.55763)}
HCNLINE={'Xb47':('spw16',354.50547),'Xb4f':('spw18',265.88618),'Xb4d':('spw27',265.88618),'Xb4b':('spw27',265.88618)}
def spec(ous,spw,reg,kind='ap'):
    z=np.load(f'/workspace/out/q5/spec_{ous}_{spw}.npz'); return z['nu'],z[reg+'_'+kind]*1e3   # mJy/beam (pbcor)
def vel(nu,f0): return C*(1-np.asarray(nu)/(f0*1e9))
def lag1(res):
    d=res-res.mean(); r=np.sum(d[1:]*d[:-1])/np.sum(d*d); r=min(max(r,0),0.9); return (1+r)/(1-r)
def gauss(v,A,v0,s): return A*np.exp(-0.5*((v-v0)/s)**2)
def line_shape(ous,reg,kind='ap',mode='profile',tname='ref'):
    """Return (shape(v_rel) unit-integral [1/(km/s)], vcen [km/s radio, w.r.t. its own rest freq], fwhm, snr). v_rel = velocity from the profile centroid.
    mode 'profile': measured profile (3-channel boxcar, baseline-subtracted, clipped >=0 and truncated at +-350 km/s from the centroid, +-450 if gaussian fails);
    mode 'gauss'  : Gaussian fitted to the reference-line spectrum (centre, FWHM free, bounds 60-400 km/s)."""
    spw,f0=(REFLINE if tname=='ref' else HCNLINE)[ous]
    nu,y=spec(ous,spw,reg,kind); v=vel(nu,f0); o=np.argsort(v); v,y=v[o],y[o]
    # initial centre: brightest 5-channel smoothed in 2400-3200 km/s
    m=(v>2300)&(v<3300); ys=np.convolve(y,np.ones(5)/5,'same'); v0=v[m][np.argmax(ys[m])]
    fl=(abs(v-v0)>500)&(abs(v-v0)<1100); pf=np.polyfit(v[fl],y[fl],1); yb=y-np.polyval(pf,v); rms=np.std(yb[fl]); 
    try:
        p,cov=curve_fit(gauss,v[abs(v-v0)<500],yb[abs(v-v0)<500],p0=[ys[m].max(),v0,70],bounds=([0,v0-200,25],[np.inf,v0+200,200]))
    except Exception: p=np.array([ys[m].max(),v0,70.]); 
    fwhm=2.355*p[2]; vc=p[1]
    snr=p[0]/rms
    if mode=='gauss':
        s=p[2]; vg=np.arange(-1500,1500.1,1.); t=np.exp(-0.5*(vg/s)**2); t/=t.sum()*(vg[1]-vg[0]); return (lambda x,vg=vg,t=t: np.interp(x,vg,t,left=0,right=0)),vc,fwhm,snr
    vg=np.arange(-1500,1500.1,1.); ysm=np.convolve(yb,np.ones(3)/3,'same'); prof=np.interp(vg+vc,v,ysm,left=0,right=0); prof=np.clip(prof,0,None)
    prof[abs(vg)>350]=0; cen=(prof*vg).sum()/prof.sum(); prof=np.interp(vg+cen,vg,prof,left=0,right=0); prof/=prof.sum()*(vg[1]-vg[0])
    return (lambda x,vg=vg,t=prof: np.interp(x,vg,t,left=0,right=0)),vc+cen,fwhm,snr
def joint_fit(ous,spw,reg,lines,shape,vcen,window=1100,order=1,kind='ap',center=None,extra=None,fixed=None):
    """Linear LSQ of shape-templated lines (+poly baseline) to the spectrum in |v - v(center line)| < window. lines: list of (name,f0 GHz) or (name,f0,weight-group).
    extra: list of (name,[(f,w),...]) composite components (fixed relative weights). Returns dict name->(W,sigma) [mJy/beam km/s... in Jy/beam km/s], rms, infl, chi."""
    nu,y=spec(ous,spw,reg,kind); f_c=center or lines[0][1]
    vc=vel(nu,f_c)-vcen; m=abs(vc)<window; v=vc[m]; yy=y[m]; nuu=nu[m]
    cols=[];names=[]
    for L in lines:
        t=shape(vel(nuu,L[1])-vcen)
        if fixed and L[0] in fixed: yy=yy-fixed[L[0]]*1e3*t; continue     # amplitude held at a given integral (Jy/beam km/s)
        cols.append(t); names.append(L[0])
    for nm,comp in (extra or []):
        cols.append(sum(w*shape(vel(nuu,f)-vcen) for f,w in comp)/sum(w for f,w in comp)); names.append(nm)
    x=(v/window)
    for k in range(order+1): cols.append(x**k)
    X=np.array(cols).T; beta,*_=np.linalg.lstsq(X,yy,rcond=None); res=yy-X@beta
    # robust noise: MAD of residual outside the lines (|v|>450) else all
    dof=len(yy)-X.shape[1]; rms=np.sqrt(np.sum(res**2)/dof); infl=lag1(res)
    cov=np.linalg.inv(X.T@X)*rms**2*infl
    out={n:(beta[i]*1e-3,np.sqrt(cov[i,i])*1e-3) for i,n in enumerate(names)}
    k=len(names); out['_rms']=rms; out['_infl']=infl; out['_rss']=np.sum(res**2)
    out['_corr']={ (names[i],names[j]):cov[i,j]/np.sqrt(cov[i,i]*cov[j,j]) for i in range(k) for j in range(i+1,k)}
    return out
def slide_sigma(ous,spw,reg,f_line,shape,vcen,window=500,kind='ap',excl=500,step=25):
    """Empirical sigma of the single-line amplitude: same template slid across the whole spw (|v - vline|>excl), linear baseline in +-window. MAD*1.4826 and std (Jy/beam km/s)."""
    nu,y=spec(ous,spw,reg,kind); v=vel(nu,f_line)-vcen; amps=[]
    for c in np.arange(v.min()+window*0.7,v.max()-window*0.7,step):
        if abs(c)<excl: continue
        m=abs(v-c)<window
        if m.sum()<100: continue
        t=shape(v[m]-c); x=(v[m]-c)/window; X=np.array([t,np.ones(m.sum()),x]).T
        b,*_=np.linalg.lstsq(X,y[m],rcond=None); amps.append(b[0]*1e-3)
    a=np.array(amps); return 1.4826*np.median(abs(a-np.median(a))),a.std(),len(a),np.median(a)

def upper_limit(ous,spw,reg,lines,name,shape,vcen,nsig=3,amax=3.0,**kw):
    """One-sided limit on amplitude `name` (>=0) from the profile chi2 (others free), chi2 scaled by rms^2*infl. Returns UL (Jy/beam km/s) at nsig, and the best a>=0."""
    base=joint_fit(ous,spw,reg,lines,shape,vcen,**kw); rms=base['_rms']; infl=base['_infl']
    sig=base[name][1]; grid=np.linspace(0,max(amax,12*sig),400); ch=[]
    for a in grid:
        r=joint_fit(ous,spw,reg,lines,shape,vcen,fixed={name:a},**kw)
        ch.append(r['_rss']/(rms**2*infl))
    ch=np.array(ch); k=ch.argmin(); d=(ch-ch[k])
    ok=np.where(d>=nsig**2)[0]; ok=ok[ok>k]
    return (grid[ok[0]] if len(ok) else np.nan),grid[k]
