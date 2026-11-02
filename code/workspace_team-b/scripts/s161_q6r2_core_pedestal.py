"""Q6 r2 fix 2: core+pedestal beam-convolved fits of the Xb4b moment-0 maps (pedestal width fixed from Xb4d), core FWHM, core flux
fraction, radial profiles. python3 -I s161_q6r2_core_pedestal.py"""
import sys, json, numpy as np, warnings
sys.path.insert(0,'/workspace/scripts'); warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter
from regions import REGIONS
C=299792.458; V0=2780.; HW=350.
LINES={'HCO+3-2':('spw25',267.5576),'HCN3-2':('spw27',265.8864)}
REG={k:REGIONS[k] for k in ('N','S','E1','E2')}
P=lambda obs,spw:'/data/member.uid___A001_X1273_%s.ngc_3256_sci.%s.cube.I.pbcor.fits'%(obs,spw)
def m0map(obs,spw,f0,ra,de,half):
    h=fits.open(P(obs,spw),memmap=True); hd=h[0].header; d=h[0].data
    while d.ndim>3: d=d[0]
    nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']; v=C*(1-nu/(f0*1e9)); dv=abs(np.median(np.diff(v)))
    pix=abs(hd['CDELT2'])*3600; w=WCS(hd).celestial; cx,cy=[float(t) for t in w.world_to_pixel_values(ra,de)]; n=int(half/pix); ix,iy=int(round(cx)),int(round(cy))
    win=np.abs(v-V0)<HW; lf=np.abs(v-V0)>600
    cut=np.array(d[:,iy-n:iy+n+1,ix-n:ix+n+1],dtype=np.float64); base=np.nanmedian(cut[lf],axis=0); sig=np.nanmedian(np.nanstd(cut[lf],axis=0))
    m0=np.nansum(cut[win]-base,axis=0)*dv; s0=sig*np.sqrt(win.sum())*dv
    return m0,pix,(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),(cx-(ix-n),cy-(iy-n)),s0
def cov(maj,minr,pa_deg):
    s1,s2=maj/2.3548,minr/2.3548; pa=np.deg2rad(pa_deg); u=np.array([-np.sin(pa),np.cos(pa)]); w=np.array([np.cos(pa),np.sin(pa)])
    return s1**2*np.outer(u,u)+s2**2*np.outer(w,w)
def gauss(X,Y,x0,y0,S):   # peak-normalised Gaussian with covariance S (arcsec^2)
    Si=np.linalg.inv(S); dx=X-x0; dy=Y-y0; return np.exp(-0.5*(Si[0,0]*dx*dx+2*Si[0,1]*dx*dy+Si[1,1]*dy*dy))
def amp_conv(S_int,Sb):   # peak factor of an integral-normalised-by-peak intrinsic Gaussian after convolution with the beam
    return np.sqrt(np.linalg.det(S_int)/np.linalg.det(S_int+Sb))
def model(p,X,Y,Sb,Wped):
    # p: Ac (intrinsic peak, Jy/beam km/s), x0,y0, maj,min,pa(deg) intrinsic core FWHM, Ap (intrinsic peak), xp,yp
    Ac,x0,y0,a,b,pa,Ap,xp,yp=p
    Sc=cov(a,b,pa); Sp=cov(Wped,Wped,0)
    return Ac*amp_conv(Sc,Sb)*gauss(X,Y,x0,y0,Sc+Sb)+Ap*amp_conv(Sp,Sb)*gauss(X,Y,xp,yp,Sp+Sb)
def fit(m,X,Y,Sb,Wped,sig,p0=None,box=None):
    ok=np.isfinite(m)&(np.hypot(X,Y)<box)
    pk=np.nanmax(m)
    k=np.argmax(np.where(np.hypot(X,Y)<0.4,m,-np.inf)); y0,x0=Y.flat[k],X.flat[k]
    if p0 is None:
        p0=[0.5*pk*2,x0,y0,0.4,0.3,30.,0.4*pk,x0,y0]
    lo=[0,x0-0.3,y0-0.3,0.02,0.02,-90,0,x0-0.6,y0-0.6]; hi=[10*pk,x0+0.3,y0+0.3,2.0,2.0,90,10*pk,x0+0.6,y0+0.6]
    p0=np.clip(p0,np.array(lo)+1e-6,np.array(hi)-1e-6)
    f=lambda p:((model(p,X[ok],Y[ok],Sb,Wped)-m[ok])/sig)
    r=least_squares(f,p0,bounds=(lo,hi),x_scale=[pk,0.1,0.1,0.2,0.2,30,pk,0.3,0.3])
    return r.x,float(np.sum(r.fun**2)),ok.sum()
def make_noise(shape,bm,pix,sig,rng):
    s=gaussian_filter(rng.normal(size=shape),bm[1]/2.3548/pix); return s/s.std()*sig
def main():
    rng=np.random.default_rng(5); out={}
    for ln,(spw,f0) in LINES.items():
        out[ln]={}
        # pedestal width from Xb4d: single elliptical Gaussian + constant per nucleus, deconvolved
        for rn,(ra,de,_) in REG.items():
            m,pix,bm,c,s0=m0map('Xb4b',spw,f0,ra,de,4.0)
            n=m.shape[0]; yy,xx=np.mgrid[:n,:n]; X=-(xx-c[0])*pix; Y=(yy-c[1])*pix    # X east-positive? use +x=east -> -(x) ; consistent w/ cov(PA) with x west: flip
            X=(xx-c[0])*pix   # x increases westward as in pixel space; cov() uses the same convention (u=(-sinPA,cosPA))
            Sb=cov(bm[0],bm[1],bm[2]); npb=1.1331*bm[0]*bm[1]/pix**2; res={}
            box=2.5 if rn in('N','S') else 1.5
            for Wped in (1.2,1.8,2.5):
                p,chi2,nn=fit(m,X,Y,Sb,Wped,s0*np.sqrt(npb)**0 ,box=box)
                # reduced chi2 inflated by independent-beam count
                chi2r=chi2/(nn/npb)/ (1) if False else chi2/ (nn-9)
                gm=np.sqrt(max(p[3]*p[4],1e-9)); 
                # core/pedestal fluxes (Jy km/s): integral of intrinsic Gaussian = A * pi/(4ln2) * a*b ; in Jy/beam km/s units convert with beam area
                Om=1.1331*bm[0]*bm[1]; Fc=p[0]*(np.pi/(4*np.log(2)))*p[3]*p[4]/Om; Fp=p[6]*(np.pi/(4*np.log(2)))*Wped**2/Om
                # statistical errors: refit on model+correlated noise
                gs=[];maj=[];mn=[]
                mod=model(p,X,Y,Sb,Wped)
                for _ in range(40):
                    q,_,_=fit(mod+make_noise(m.shape,bm,pix,s0,rng),X,Y,Sb,Wped,s0,p0=list(p),box=box); gs.append(np.sqrt(q[3]*q[4])); maj.append(max(q[3],q[4]))
                res[str(Wped)]=dict(core_a=float(p[3]),core_b=float(p[4]),core_pa=float(p[5]),core_gm=float(gm),gm_p16_84=np.percentile(gs,[16,84]).tolist(),
                    core_peak_int=float(p[0]),ped_peak_int=float(p[6]),Fcore=float(Fc),Fped=float(Fp),core_frac_total=float(Fc/(Fc+Fp)),chi2r=float(chi2/(nn/npb)),
                    ped_off=[float(p[7]-p[1]),float(p[8]-p[2])],core_obs_peak_frac=float(amp_conv(cov(p[3],p[4],p[5]),Sb)*p[0]/(amp_conv(cov(p[3],p[4],p[5]),Sb)*p[0]+amp_conv(cov(Wped,Wped,0),Sb)*p[6])))
            # radial profile (peak-centred, annuli 0.1") of data and W=1.8 model
            Wf=1.8; p,_,_=fit(m,X,Y,Sb,Wf,s0,box=box); mod=model(p,X,Y,Sb,Wf); r=np.hypot(X-p[1],Y-p[2]); edges=np.arange(0,2.5,0.125)
            prof=[[float(0.5*(e0+e1)),float(np.nanmean(m[(r>=e0)&(r<e1)])),float(np.nanmean(mod[(r>=e0)&(r<e1)]))] for e0,e1 in zip(edges[:-1],edges[1:])]
            res['profile']=prof; res['beam']=bm; res['sig0']=s0
            out[ln][rn]=res
            print(ln,rn,' '.join('W%s: core %.2fx%.2f gm %.2f [%.2f,%.2f] cf %.2f chi2r %.1f'%(w,r_['core_a'],r_['core_b'],r_['core_gm'],*r_['gm_p16_84'],r_['core_frac_total'],r_['chi2r']) for w,r_ in res.items() if w[0].isdigit()),flush=True)
    json.dump(out,open('/workspace/out/q6r2/s161.json','w'),indent=1)
main()
