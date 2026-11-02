# CH3OH 2(0)-1(-1)E 254.0154 GHz: moment-0 maps (Xb4d 0.95", Xb4f 6x3.7"), HCO+ 3-2 in the same window, matched-beam
# comparison with the Xb55 J=1-3 E blend, and an LTE excitation test. Window: CH3OH-frame velocity 2783 +- 150 km/s
# (= HCO+ 3-2 v0 at the Bridge, +-1.3 FWHM). Baseline: per-pixel linear fit to line-free channels.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
from regions import REGIONS
C=299792.458; REST=254.015377; HCO=267.557626; V0=2783.; HW=150.
def load(ous,spw):
    h=fits.open(f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.spw{spw}.cube.I.pbcor.fits',memmap=True)
    H=h[0].header; d=np.squeeze(h[0].data).astype(float)
    nu=(H['CRVAL3']+(np.arange(H['NAXIS3'])-H['CRPIX3']+1)*H['CDELT3'])/1e9
    return d,H,nu
def mom0(d,H,nu,rest,v0,hw,free=((1500,2150),(3450,4100))):
    v=C*(1-nu/rest); dv=abs(C*H['CDELT3']/1e9/rest)
    sel=np.zeros(len(v),bool)
    for a,b in free: sel|=(v>a)&(v<b)
    win=np.abs(v-v0)<hw
    ny,nx=d.shape[1:]; X=np.vander(v[sel]-3000,2)
    Y=d[sel].reshape(sel.sum(),-1); coef=np.linalg.lstsq(X,Y,rcond=None)[0]
    base=(np.vander(v[win]-3000,2)@coef).reshape(win.sum(),ny,nx)
    m=((d[win]-base).sum(0))*dv
    # empirical noise: same-width windows slid over line-free channels (baseline-subtracted)
    res=(Y-X@coef).reshape(-1,ny,nx); n=win.sum()
    sig=np.sqrt(np.nanmean(res**2,0))*np.sqrt(n)*dv   # white-noise estimate per pixel
    return m,sig,dv
def beam(H): return H['BMAJ']*3600,H['BMIN']*3600,H['BPA']
def kern_sigma_pix(bA,bB,pix):
    # circularised deconvolution (areas) -> kernel FWHM, Gaussian sigma in pixels; elliptical handled by two axes below
    return None
out={}
for ous,(sp254,spH) in {'Xb4d':('29','25'),'Xb4f':('20','16')}.items():
    dA,HA,nuA=load(ous,sp254); dB,HB,nuB=load(ous,spH)
    mA,sA,dvA=mom0(dA,HA,nuA,REST,V0,HW); mB,sB,dvB=mom0(dB,HB,nuB,HCO,V0-9,HW)   # HCO+ window centred at its own v0 (CH3OH frame +9)
    w=WCS(HA).celestial; pix=abs(HA['CDELT2'])*3600
    print(f'== {ous}: beam254 {beam(HA)[0]:.2f}x{beam(HA)[1]:.2f}; pixel {pix:.3f}"; median per-pixel sigma(I254) {np.nanmedian(sA):.4f} Jy/b km/s')
    # sigma map from empirical scatter of the moment-0 map outside the central region
    ny,nx=mA.shape; yy,xx=np.mgrid[:ny,:nx]; cx,cy=w.all_world2pix(156.96346,-43.90461,0)
    r_=np.hypot(xx-cx,yy-cy)*pix; far=(r_>(4 if ous=='Xb4d' else 7))&(r_<(11 if ous=='Xb4d' else 14))&np.isfinite(mA)&(mA!=0)
    mad=1.4826*np.nanmedian(np.abs(mA[far]-np.nanmedian(mA[far]))); print(f'   annulus MAD of mom0 = {mad:.4f} (white-noise est {np.nanmedian(sA[far]):.4f}); n={far.sum()}')
    near=np.hypot(xx-cx,yy-cy)*pix<12
    # peaks
    from scipy.ndimage import maximum_filter
    pk=(mA==maximum_filter(mA,size=int(round(1.0/pix))*2+1))&near&(mA>4*mad)
    rows=[]
    for y,x in zip(*np.where(pk)):
        ra,dec=w.all_pix2world(x,y,0); 
        near_reg=min(REGIONS,key=lambda k:np.hypot((REGIONS[k][0]-ra)*np.cos(np.radians(dec)),REGIONS[k][1]-dec))
        dist=np.hypot((REGIONS[near_reg][0]-ra)*np.cos(np.radians(dec)),REGIONS[near_reg][1]-dec)*3600
        rows.append((mA[y,x],ra,dec,near_reg,dist,mB[y,x]))
    rows.sort(reverse=True)
    for r in rows[:6]: print(f'   peak {r[0]:.3f} Jy/b km/s ({r[0]/mad:.1f} sigma_MAD) at {r[1]:.5f},{r[2]:.5f} near {r[3]} ({r[4]:.1f}")  HCO+3-2 {r[5]:.2f}  ratio(Jy)={r[0]/r[5]:.4f}')
    # region values
    ν254=REST/(1-V0/C); nuH=HCO/(1-V0/C)
    reg_out={}
    for k in ('N','S','Bridge','Sblue','W','E1','E2'):
        x,y=w.all_world2pix(REGIONS[k][0],REGIONS[k][1],0); x=int(round(float(x)));y=int(round(float(y)))
        if not(0<=x<nx and 0<=y<ny): continue
        a,b_=mA[y,x],mB[y,x]
        # brightness-temperature-integral ratio: K km/s ratio = (S254/nu254^2)/(SH/nuH^2) in matched beam
        kr=(a/ν254**2)/(b_/nuH**2) if b_ else np.nan
        reg_out[k]=dict(I254=a,err_mad=mad,snr=a/mad,IHCO32=b_,ratio_Jy=a/b_,ratio_K=kr,err_ratio=mad/b_)
        print(f'   {k:6s} I(254)={a:+.3f}+-{mad:.3f} ({a/mad:+.1f} sig)  I(HCO+3-2)={b_:.2f}  CH3OH/HCO+ (Jy/b km/s)={a/b_:.4f}+-{mad/b_:.4f}  (K km/s) {kr:.4f}')
    out[ous]=reg_out
    np.savez(f'/workspace/out/s21_ch3oh254_mom0_{ous}.npz',ch3oh=mA,hco=mB,mad=mad,wcs=w.to_header_string())
    if ous=='Xb4f':
        # matched beam: smooth Xb4f to the Xb55 beam (covariance deconvolution, pixel units), compare with Xb55 J=1-3 blend
        b55=(9.357,6.284,-88.45); b4=beam(HA)
        def cov(b,pixs):
            a=b[0]/2.3548/pixs;c=b[1]/2.3548/pixs; t=np.radians(b[2]+90)  # PA east of north -> use x=east? (RA axis flipped): handle by angle from +x
            R=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]]); return R@np.diag([a*a,c*c])@R.T
        # PA measured from North toward East; image x increases toward West, y toward North. major-axis direction vector in (x,y): (-sin PA, cos PA)
        def cov2(b):
            a=b[0]/2.3548/pix;c=b[1]/2.3548/pix; t=np.radians(b[2]); u=np.array([-np.sin(t),np.cos(t)]); v=np.array([np.cos(t),np.sin(t)])
            return a*a*np.outer(u,u)+c*c*np.outer(v,v)
        K=cov2(b55)-cov2(b4); ev,evec=np.linalg.eigh(K); print('   kernel eigenvalues (pix^2):',ev)
        ev=np.clip(ev,1e-6,None)
        # apply by rotating coordinates: use scipy affine via explicit convolution kernel
        n=41; g=np.mgrid[-n//2+1:n//2+1,-n//2+1:n//2+1]; P=np.stack([g[1].ravel(),g[0].ravel()]); Ki=np.linalg.inv(K+1e-6*np.eye(2))
        ker=np.exp(-0.5*np.einsum('ij,ik,kj->j',P,Ki,P)).reshape(n,n); ker/=ker.sum()
        from scipy.signal import fftconvolve
        fac=(b55[0]*b55[1])/(b4[0]*b4[1])   # Jy/beam scales with beam area for a smooth map convolved & renormalised
        sm=lambda m: fftconvolve(np.nan_to_num(m),ker,mode='same')*fac
        smA=sm(mA); smB=sm(mB)
        sm_out={}
        for k in ('N','S','Bridge'):
            x,y=w.all_world2pix(REGIONS[k][0],REGIONS[k][1],0); x=int(round(float(x)));y=int(round(float(y)))
            sm_out[k]=dict(I254_55=smA[y,x],IHCO32_55=smB[y,x]); print(f'   smoothed to Xb55 beam: {k} I254={smA[y,x]:.3f}  HCO+3-2={smB[y,x]:.1f} Jy/b km/s')
        # also noise of smoothed map
        sA_far=1.4826*np.median(np.abs(smA[far]-np.median(smA[far]))); print('   smoothed far-field MAD',sA_far)
        out['smoothed']=sm_out; out['smoothed_mad']=sA_far
json.dump(out,open('/workspace/out/s21_ch3oh254.json','w'),indent=1,default=float)
