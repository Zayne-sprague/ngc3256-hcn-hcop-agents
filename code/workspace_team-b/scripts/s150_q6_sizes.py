"""Q6: deconvolved sizes, line widths, Tb and beam filling of dense-gas structures at 0.24" (Xb4b).
2D elliptical Gaussian + constant fits to moment-0 (HCO+ 3-2 spw25, HCN 3-2 spw27) and to the 4-spw continuum,
beam-deconvolved, errors scaled by correlated-noise (Condon 1997: Npix per beam). Run: python3 -I s150_q6_sizes.py"""
import sys, json, numpy as np, warnings
sys.path.insert(0,'/workspace/scripts')
warnings.filterwarnings('ignore')
from astropy.io import fits
from scipy.optimize import curve_fit
from regions import REGIONS
D='/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.'
C=299792.458; V0=2780.; HW=350.
DIST=40.0  # Mpc (Sakamoto+2014 / Pereira-Santaella+2024: 5" = 850 pc => 170 pc/")
PC=194.0
LINES={'HCO+3-2':('spw25',267.5576),'HCN3-2':('spw27',265.8864)}
REG=dict(REGIONS); REG['SE']=(156.96684,-43.90635,'SE continuum source')
def gfun(xy,A,x0,y0,a,b,pa,c):   # a,b FWHM (arcsec) ; pa rad from +y(N) toward +x(E)? use math angle from x axis
    x,y=xy; ct,st=np.cos(pa),np.sin(pa)
    u=(x-x0)*ct+(y-y0)*st; v=-(x-x0)*st+(y-y0)*ct
    s=4*np.log(2)
    return A*np.exp(-s*(u/a)**2-s*(v/b)**2)+c
def deconv(a,b,pa,ba,bb,bpa):
    """Gaussian deconvolution (Condon 1997 / AIPS DECON), sizes FWHM, angles from x axis [rad]. returns maj,min,pa (<0 sizes => 0)"""
    a2,b2=a*a,b*b; ba2,bb2=ba*ba,bb*bb
    ca=a2*np.cos(2*pa)-b2*0 ; # build via second-moment matrices
    def M(a,b,t): 
        R=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]]); return R@np.diag([a*a,b*b])@R.T
    S=M(a,b,pa)-M(ba,bb,bpa)
    w,V=np.linalg.eigh(S)
    mj=np.sqrt(max(w[1],0)); mn=np.sqrt(max(w[0],0)); t=np.arctan2(V[1,1],V[0,1])
    return mj,mn,t
def getplane(fn,x0,y0,half,kind,v0=V0,hw=HW,nu0=None):
    h=fits.open(fn,memmap=True); hd=h[0].header
    pix=abs(hd['CDELT2'])*3600.; ra0,de0=hd['CRVAL1'],hd['CRVAL2']; cx,cy=hd['CRPIX1']-1,hd['CRPIX2']-1
    return h,hd,pix,ra0,de0,cx,cy
def pos2pix(hd,ra,de):
    from astropy.wcs import WCS
    w=WCS(hd).celestial; x,y=w.world_to_pixel_values(ra,de); return float(x),float(y)
def fit_region(img,hd,pix,ra,de,beam,box,sigma,label):
    cx,cy=pos2pix(hd,ra,de); n=int(round(box/pix))
    ix,iy=int(round(cx)),int(round(cy))
    sub=img[iy-n:iy+n+1,ix-n:ix+n+1]
    # RA axis increases to the left: x_arcsec(East) = -(pixel offset)*pix
    X,Y=np.meshgrid(-(np.arange(-n,n+1))*pix,(np.arange(-n,n+1))*pix)  # X = East offset, Y = North offset
    ok=np.isfinite(sub)
    # start at the max within 0.3" of the nominal position
    r=np.hypot(X,Y); m=np.where(r<0.4,sub,-np.inf); k=np.nanargmax(m); yk,xk=np.unravel_index(k,sub.shape)
    p0=[sub[yk,xk],X[yk,xk],Y[yk,xk],beam[0]*1.2,beam[1]*1.2,0.5,0.0]
    lo=[0,-box,-box,0.3*beam[1],0.3*beam[1],-np.pi,-np.inf]; hi=[np.inf,box,box,box*2,box*2,np.pi,np.inf]
    try:
        p,cov=curve_fit(gfun,(X[ok],Y[ok]),sub[ok],p0=p0,bounds=(lo,hi),sigma=np.full(ok.sum(),sigma),absolute_sigma=True,maxfev=20000)
    except Exception as e:
        return None
    res=sub[ok]-gfun((X[ok],Y[ok]),*p)
    npb=1.1331*beam[0]*beam[1]/pix**2
    cov=cov*npb    # correlated-noise inflation of the covariance (Condon 1997)
    e=np.sqrt(np.diag(cov))
    return dict(label=label,p=p,e=e,npb=npb,chi=float(np.std(res)/sigma),cx=cx,cy=cy,X0=p[1],Y0=p[2])
def analyse(p,e,beam,rng):
    """beam=(maj,min,pa_math). returns deconvolved maj/min (arcsec) with MC errors, flag unresolved"""
    A,x0,y0,a,b,pa,c=p; ea,eb=e[3],e[4]
    if b>a: a,b=b,a; ea,eb=eb,ea; pa=pa+np.pi/2
    pa_m=(pa)%np.pi
    mj,mn,t=deconv(a,b,pa_m,*beam)
    ms=[]
    for _ in range(400):
        aa=max(rng.normal(a,ea),0.05); bb=max(rng.normal(b,eb),0.05)
        if bb>aa: aa,bb=bb,aa
        ms.append(deconv(aa,bb,pa_m,*beam)[:2])
    ms=np.array(ms)
    return mj,mn,ms
def main():
    rng=np.random.default_rng(1); out={}
    # ---------------- continuum
    cont=D+'spw25_27_29_31.cont.I.pbcor.fits'
    h=fits.open(cont,memmap=True); hd=h[0].header; img=np.squeeze(h[0].data).astype(float)
    pix=abs(hd['CDELT2'])*3600.; bm=(hd['BMAJ']*3600,hd['BMIN']*3600,np.deg2rad(90-hd['BPA']))
    # sigma from MAD of the map far from sources
    yy,xx=np.mgrid[:img.shape[0],:img.shape[1]]; rr=np.hypot(xx-img.shape[1]/2,yy-img.shape[0]/2)*pix
    sig_c=1.4826*np.nanmedian(np.abs(img[(rr>8)&(rr<14)]-np.nanmedian(img[(rr>8)&(rr<14)])))
    print('continuum beam',bm[0],bm[1],np.rad2deg(bm[2]),'rms mJy/beam',sig_c*1e3)
    out['cont']={'beam':[bm[0],bm[1],float(np.rad2deg(bm[2]))],'rms_mJy':sig_c*1e3,'fits':{}}
    nu_c=257.760
    for name,(ra,de,_) in REG.items():
        for box in (0.6,1.0,1.5):
            r=fit_region(img,hd,pix,ra,de,bm,box,sig_c,name)
            if r is None: continue
            mj,mn,ms=analyse(r['p'],r['e'],bm,rng)
            A=r['p'][0]; 
            out['cont']['fits'].setdefault(name,{})[str(box)]=dict(Apk_mJy=A*1e3,eA=r['e'][0]*1e3,snr=A/sig_c,fwhm_maj=float(max(r['p'][3],r['p'][4])),fwhm_min=float(min(r['p'][3],r['p'][4])),dec_maj=mj,dec_min=mn,dec_maj_p16_84=np.percentile(ms[:,0],[16,84]).tolist(),dec_min_p16_84=np.percentile(ms[:,1],[16,84]).tolist(),dec_maj_p95=float(np.percentile(ms[:,0],95)),dec_min_p95=float(np.percentile(ms[:,1],95)),bkg_mJy=r['p'][6]*1e3,offs=[r['X0'],r['Y0']],Sint_mJy=float(A*1e3*r['p'][3]*r['p'][4]/(bm[0]*bm[1])))
    # ---------------- lines
    out['lines']={}
    for ln,(spw,f0) in LINES.items():
        fn=D+spw+'.cube.I.pbcor.fits'; h=fits.open(fn,memmap=True); hd=h[0].header
        data=h[0].data; 
        while data.ndim>3: data=data[0]
        nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']; v=C*(1-nu/(f0*1e9))
        dv=abs(np.median(np.diff(v))); pix=abs(hd['CDELT2'])*3600.
        bm=(hd['BMAJ']*3600,hd['BMIN']*3600,np.deg2rad(90-hd['BPA']))
        win=np.abs(v-V0)<HW; lf=(np.abs(v-V0)>600)
        idx=np.where(win)[0]; i0,i1=idx.min(),idx.max()+1
        print(ln,'chan dv',dv,'nwin',win.sum(),'v range',v[win].min(),v[win].max(),'beam',bm[0],bm[1])
        res={'beam':[bm[0],bm[1],float(np.rad2deg(bm[2]))],'dv':dv,'fits':{},'spec':{}}
        nu_l=nu0=f0/(1+V0/C*0)  # placeholder
        for name,(ra,de,_) in REG.items():
            cx,cy=pos2pix(hd,ra,de); pad=int(2.5/pix); ix,iy=int(round(cx)),int(round(cy))
            if ix-pad<0 or iy-pad<0 or ix+pad>=data.shape[2] or iy+pad>=data.shape[1]: continue
            ys,xs=slice(iy-pad,iy+pad+1),slice(ix-pad,ix+pad+1)
            cut=np.array(data[:,ys,xs],dtype=float)
            m0=np.nansum(cut[i0:i1],axis=0)*dv                      # Jy/beam km/s
            lfc=np.nanstd(cut[lf],axis=0) if lf.sum()>10 else None
            sig_ch=np.nanmedian(lfc) ; sig0=sig_ch*np.sqrt(win.sum())*dv
            # embed in a big array with the same geometry for fit_region
            big=np.full((2*pad+1,2*pad+1),np.nan); big[:]=m0
            class H: pass
            hh=dict(hd); 
            cxl,cyl=cx-(ix-pad),cy-(iy-pad)
            for box in (0.6,1.0,1.5):
                nn=int(round(box/pix)); 
                r=None
                sub_hd=None
                # local fit (own coordinates): reuse fit_region with a fake wcs through pixel position
                from types import SimpleNamespace
                class W: 
                    pass
                r=fit_local(big,cxl,cyl,pix,bm,box,sig0,name)
                if r is None: continue
                mj,mn,ms=analyse(r['p'],r['e'],bm,rng)
                A=r['p'][0]
                res['fits'].setdefault(name,{})[str(box)]=dict(Apk=A,eA=r['e'][0],snr=A/sig0,fwhm_maj=float(max(r['p'][3],r['p'][4])),fwhm_min=float(min(r['p'][3],r['p'][4])),e_maj=float(r['e'][3]),dec_maj=mj,dec_min=mn,dec_maj_p16_84=np.percentile(ms[:,0],[16,84]).tolist(),dec_min_p16_84=np.percentile(ms[:,1],[16,84]).tolist(),dec_maj_p95=float(np.percentile(ms[:,0],95)),dec_min_p95=float(np.percentile(ms[:,1],95)),bkg=r['p'][6],offs=[r['X0'],r['Y0']],Sint=float(A*r['p'][3]*r['p'][4]/(bm[0]*bm[1])),pa_deg=float((90-np.rad2deg(r['p'][5]))%180),chi=r['chi'])
            # spectra: peak pixel of m0 near region (fit centre) and 0.25"-radius aperture mean
            kk=np.argmax(np.where(np.hypot((np.arange(2*pad+1)[None,:]-cxl),(np.arange(2*pad+1)[:,None]-cyl))*pix<0.4,m0,-np.inf)); py,px=np.unravel_index(kk,m0.shape)
            sp_pk=cut[:,py,px]
            yy2,xx2=np.mgrid[:2*pad+1,:2*pad+1]; ap=np.hypot(xx2-px,yy2-py)*pix<0.25
            sp_ap=np.nanmean(cut[:,ap],axis=1)
            res['spec'][name]=dict(pix_offset_arcsec=[float(-(px-cxl)*pix),float((py-cyl)*pix)],sig_ch=float(sig_ch),pk=sp_pk.tolist(),ap=sp_ap.tolist(),v=v.tolist() if name=='N' else None)
            np.save(f'/workspace/out/q6/m0_{ln}_{name}.npy',m0)
        out['lines'][ln]=res
    json.dump(out,open('/workspace/out/q6/s150.json','w'),default=float)
def fit_local(sub,cx,cy,pix,beam,box,sigma,label):
    n=int(round(box/pix)); ix,iy=int(round(cx)),int(round(cy))
    s=sub[iy-n:iy+n+1,ix-n:ix+n+1]
    X,Y=np.meshgrid(-(np.arange(-n,n+1)+ix-cx)*pix,(np.arange(-n,n+1)+iy-cy)*pix)
    ok=np.isfinite(s); r=np.hypot(X,Y); m=np.where(r<0.4,s,-np.inf); k=np.nanargmax(m); yk,xk=np.unravel_index(k,s.shape)
    p0=[s[yk,xk],X[yk,xk],Y[yk,xk],beam[0]*1.2,beam[1]*1.2,0.5,0.0]
    lo=[0,-box,-box,0.3*beam[1],0.3*beam[1],-np.pi,-np.inf]; hi=[np.inf,box,box,box*2,box*2,np.pi,np.inf]
    try:
        p,cov=curve_fit(gfun,(X[ok],Y[ok]),s[ok],p0=p0,bounds=(lo,hi),sigma=np.full(ok.sum(),sigma),absolute_sigma=True,maxfev=20000)
    except Exception as e:
        return None
    resid=s[ok]-gfun((X[ok],Y[ok]),*p); npb=1.1331*beam[0]*beam[1]/pix**2
    e=np.sqrt(np.diag(cov*npb))
    return dict(label=label,p=p,e=e,chi=float(np.std(resid)/sigma),X0=p[1],Y0=p[2])
if __name__=='__main__': main()
