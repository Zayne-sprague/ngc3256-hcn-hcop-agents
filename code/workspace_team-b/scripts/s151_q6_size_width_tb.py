"""Q6 part 2: non-parametric (half-max area) sizes, compact+extended 2-Gaussian fits, line widths (peak pixel / 0.25" aperture,
beam-smearing correction from the local velocity gradient), brightness temperatures. python3 -I s151_q6_size_width_tb.py"""
import sys, json, numpy as np, warnings
sys.path.insert(0,'/workspace/scripts'); warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import curve_fit
from scipy.ndimage import label, gaussian_filter
from regions import REGIONS
D='/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.'
C=299792.458; V0=2780.; HW=350.; PC=194.0
LINES={'HCO+3-2':('spw25',267.5576),'HCN3-2':('spw27',265.8864)}
REG=dict(REGIONS); REG['SE']=(156.96684,-43.90635,'SE')
H=6.62607e-27; KB=1.380649e-16
def tb_rj(S_jy_beam,nu_ghz,bmaj,bmin):   # K ; beam in arcsec
    return 1.222e6*S_jy_beam*1e3/ (nu_ghz**2*bmaj*bmin) /1e3*1.0 if False else 1.222e3*(S_jy_beam*1e3)/(nu_ghz**2*bmaj*bmin)
def tb_planck(trj,nu_ghz):
    x=H*nu_ghz*1e9/KB
    return x/np.log1p(x/trj) if trj>0 else np.nan
def halfmax(m,pix,cx,cy,rmax,frac=0.5):
    ys,xs=np.indices(m.shape); r=np.hypot(xs-cx,ys-cy)*pix
    sm=gaussian_filter(m,0.8)
    k=np.argmax(np.where(r<0.4,sm,-np.inf)); py,px=np.unravel_index(k,m.shape)
    pk=sm[py,px]; msk=(sm>frac*pk)&(np.hypot(xs-px,ys-py)*pix<rmax)
    lab,_=label(msk); l=lab[py,px]; a=(lab==l).sum()*pix**2
    return a,pk,(py,px)
def make_noise(shape,beam,pix,sig,rng):
    s=gaussian_filter(rng.normal(size=shape),(beam[1]/2.355/pix))   # circular smoothing with minor sigma
    return s/ s.std()*sig
def g2(xy,A,x0,y0,a,b,pa,c,A2,a2,b2):
    x,y=xy; s=4*np.log(2)
    def one(A,a,b):
        ct,st=np.cos(pa),np.sin(pa); u=(x-x0)*ct+(y-y0)*st; v=-(x-x0)*st+(y-y0)*ct
        return A*np.exp(-s*(u/a)**2-s*(v/b)**2)
    return one(A,a,b)+one(A2,a2,b2)+c
def main():
    rng=np.random.default_rng(2); out={}
    # continuum half-max
    fn=D+'spw25_27_29_31.cont.I.pbcor.fits'; h=fits.open(fn,memmap=True); hd=h[0].header; img=np.squeeze(h[0].data).astype(float)
    pix=abs(hd['CDELT2'])*3600; bm=(hd['BMAJ']*3600,hd['BMIN']*3600); w=WCS(hd).celestial
    yy,xx=np.mgrid[:img.shape[0],:img.shape[1]]; rr=np.hypot(xx-img.shape[1]/2,yy-img.shape[0]/2)*pix
    sig=1.4826*np.nanmedian(np.abs(img[(rr>8)&(rr<14)]-np.nanmedian(img[(rr>8)&(rr<14)])))
    out['cont']={}
    for nme,(ra,de,_) in REG.items():
        cx,cy=[float(t) for t in w.world_to_pixel_values(ra,de)]; n=int(3/pix); ix,iy=int(round(cx)),int(round(cy))
        sub=img[iy-n:iy+n+1,ix-n:ix+n+1]; cxl,cyl=cx-(ix-n),cy-(iy-n)
        out['cont'][nme]=hm_stats(sub,pix,cxl,cyl,bm,sig,rng,rmax=1.2 if nme in('N','Bridge') else 1.5)
    # lines
    out['lines']={}
    for ln,(spw,f0) in LINES.items():
        h=fits.open(D+spw+'.cube.I.pbcor.fits',memmap=True); hd=h[0].header; data=h[0].data
        while data.ndim>3: data=data[0]
        nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']; v=C*(1-nu/(f0*1e9)); dv=abs(np.median(np.diff(v)))
        pix=abs(hd['CDELT2'])*3600; bm=(hd['BMAJ']*3600,hd['BMIN']*3600); w=WCS(hd).celestial
        win=np.abs(v-V0)<HW; lf=np.abs(v-V0)>600; i0,i1=np.where(win)[0][[0,-1]]; i1+=1
        nu_obs_mid=f0*(1-V0/C); out['lines'][ln]={}
        for nme,(ra,de,_) in REG.items():
            cx,cy=[float(t) for t in w.world_to_pixel_values(ra,de)]; pad=int(2.5/pix); ix,iy=int(round(cx)),int(round(cy))
            cut=np.array(data[:,iy-pad:iy+pad+1,ix-pad:ix+pad+1],dtype=float); cxl,cyl=cx-(ix-pad),cy-(iy-pad)
            base=np.nanmedian(cut[lf],axis=0); cut=cut-0  # cubes are already continuum subtracted; line-free median checked below
            sig_ch=np.nanmedian(np.nanstd(cut[lf],axis=0)); 
            m0=np.nansum(cut[i0:i1]-base,axis=0)*dv; sig0=sig_ch*np.sqrt(win.sum())*dv
            st=hm_stats(m0,pix,cxl,cyl,bm,sig0,rng,rmax=1.2 if nme in('N','Bridge') else 1.5)
            # peak pixel (of smoothed m0 within 0.4" of nominal)
            ys,xs=np.indices(m0.shape); r=np.hypot(xs-cxl,ys-cyl)*pix; sm=gaussian_filter(m0,0.8)
            k=np.argmax(np.where(r<0.4,sm,-np.inf)); py,px=np.unravel_index(k,m0.shape)
            spk=cut[:,py,px]-base[py,px]; ap=np.hypot(xs-px,ys-py)*pix<0.25; sap=np.nanmean(cut[:,ap]-base[ap],axis=1)
            sp={}
            for key,s in (('peakpix',spk),('ap0.25',sap)):
                sp[key]=spec_props(v,s,sig_ch/(1 if key=='peakpix' else np.sqrt(ap.sum()/ (1.1331*bm[0]*bm[1]/pix**2)+1)),nu_obs_mid,bm,dv)
            # velocity gradient near the peak from moment-1 (5 sigma pixels, window)
            I=np.nansum((cut[i0:i1]-base)*v[i0:i1,None,None],axis=0)*dv
            m1=np.where(m0>5*sig0,I/np.where(m0>5*sig0,m0,1),np.nan)
            sel=np.isfinite(m1)&(np.hypot(xs-px,ys-py)*pix<0.4)
            grad=None
            if sel.sum()>=12:
                A=np.c_[(xs[sel]-px)*pix,-(ys[sel]-py)*pix,np.ones(sel.sum())]; cf=np.linalg.lstsq(A,m1[sel],rcond=None)[0]
                G=float(np.hypot(cf[0],cf[1])); grad=dict(G_kms_per_arcsec=G,theta_deg=float(np.rad2deg(np.arctan2(cf[0],cf[1]))%360),npix=int(sel.sum()),v_mean=float(cf[2]))
                bsm=np.sqrt(bm[0]*bm[1])
                for key in sp:
                    s2=sp[key]['sigma_mom2']
                    sg=G*bsm/2.355; sp[key]['sigma_grad']=sg
                    sp[key]['sigma_int_mom2']=float(np.sqrt(max(s2**2-sg**2-(dv/2.355)**2,0)))
            out['lines'][ln][nme]=dict(size=st,spec=sp,grad=grad,sig_ch_mJy=sig_ch*1e3,sig0=sig0,off=[float((px-cxl)*pix*-1),float((py-cyl)*pix)])
        out['lines'][ln]['_beam']=bm; out['lines'][ln]['_dv']=dv
    json.dump(out,open('/workspace/out/q6/s151.json','w'),default=float)
def cross(v,sm,j,frac):
    half=frac*sm[j]; lo=j
    while lo>0 and sm[lo]>half: lo-=1
    hi=j
    while hi<len(sm)-1 and sm[hi]>half: hi+=1
    vlo=np.interp(half,[sm[lo],sm[lo+1]],[v[lo],v[lo+1]]) if sm[lo+1]!=sm[lo] else v[lo]
    vhi=np.interp(half,[sm[hi],sm[hi-1]],[v[hi],v[hi-1]]) if sm[hi-1]!=sm[hi] else v[hi]
    return vlo,vhi
def spec_props(v,s,sig,nu_obs_ghz,bm,dv):
    sm=np.convolve(s,np.ones(3)/3,'same'); win=np.abs(v-V0)<HW
    j=np.where(win)[0][np.argmax(sm[win])]; Ipk=sm[j]       # 3-channel boxcar peak (limits the noise bias to ~0.3 sigma)
    vlo,vhi=cross(v,sm,j,0.5); fwhm=abs(vhi-vlo)
    v2l,v2h=cross(v,sm,j,0.2); w20=abs(v2h-v2l)
    # noise-robust 2nd moment: contiguous channels around the peak above 2 sigma_3ch
    s3=sig/np.sqrt(3); lo=j
    while lo>0 and sm[lo]>2*s3: lo-=1
    hi=j
    while hi<len(sm)-1 and sm[hi]>2*s3: hi+=1
    sel=slice(lo+1,hi); w_=s[sel]; vv=v[sel]; cen=(w_*vv).sum()/w_.sum(); sig2=np.sqrt((w_*(vv-cen)**2).sum()/w_.sum())
    # single Gaussian fit within +-300 km/s of vsys
    m=np.abs(v-V0)<300
    try:
        f=lambda x,A,x0,sg: A*np.exp(-0.5*((x-x0)/sg)**2)
        p,cv=curve_fit(f,v[m],s[m],p0=[Ipk,v[j],40.],sigma=np.full(m.sum(),sig),absolute_sigma=True)
        g=dict(A_mJy=float(p[0]*1e3),v0=float(p[1]),fwhm=float(abs(p[2])*2.355),efwhm=float(np.sqrt(cv[2,2])*2.355),chi2r=float(np.sum(((s[m]-f(v[m],*p))/sig)**2)/(m.sum()-3)))
    except Exception: g=None
    trj=tb_rj(Ipk,nu_obs_ghz,*bm)
    return dict(peak_mJy=float(Ipk*1e3),v_peak=float(v[j]),fwhm_halfmax=float(fwhm),w20=float(w20),vlo=float(vlo),vhi=float(vhi),centroid=float(cen),sigma_mom2=float(sig2),n_mom2_chan=int(hi-lo-1),gauss1=g,
                Tb_RJ=float(trj),eTb=float(tb_rj(s3,nu_obs_ghz,*bm)),Tb_planck=float(tb_planck(trj,nu_obs_ghz)),snr=float(Ipk/s3),I_Kkms=float(tb_rj(np.nansum(s[win])*dv,nu_obs_ghz,*bm)))
def hm_stats(m,pix,cx,cy,bm,sig,rng,rmax=1.5):
    a,pk,(py,px)=halfmax(m,pix,cx,cy,rmax)
    res=[]
    for _ in range(150):
        mm=m+make_noise(m.shape,bm,pix,sig,rng); aa,pp,_=halfmax(mm,pix,cx,cy,rmax); res.append(aa)
    res=np.array(res); b2=bm[0]*bm[1]
    deq=lambda A: 2*np.sqrt(A/np.pi)
    dd=np.sqrt(np.clip(deq(res)**2-b2,0,None)); d0=np.sqrt(max(deq(a)**2-b2,0))
    out=dict(peak_sm=float(pk),Ahm_arcsec2=float(a),Deq_obs=float(deq(a)),Ddec=float(d0),Ddec_p16_84=np.percentile(dd,[16,84]).tolist(),Ddec_p95=float(np.percentile(dd,95)),Deq_boot=float(deq(res).mean()))
    # Gaussian + constant elliptical single-component fit baseline vs 2-component (compact+extended, same centre/orientation)
    n=min(int(1.5/pix),(m.shape[0]-1)//2-1); ix,iy=int(round(cx)),int(round(cy)); s=m[iy-n:iy+n+1,ix-n:ix+n+1]
    X,Y=np.meshgrid(-(np.arange(-n,n+1)+ix-cx)*pix,(np.arange(-n,n+1)+iy-cy)*pix); ok=np.isfinite(s)
    npb=1.1331*bm[0]*bm[1]/pix**2
    try:
        p0=[0.7*pk,0,0,bm[0]*1.3,bm[1]*1.3,0.5,0,0.3*pk,1.2,1.0]
        lo=[0,-0.3,-0.3,0.8*bm[1],0.8*bm[1],-np.pi,-np.inf,0,0.6,0.6]; hi=[np.inf,0.3,0.3,0.7,0.7,np.pi,np.inf,np.inf,4,4]
        # shift start to the actual peak
        ky=np.argmax(np.where(np.hypot(X,Y)<0.4,s,-np.inf)); yk,xk=np.unravel_index(ky,s.shape); p0[1],p0[2]=X[yk,xk],Y[yk,xk]
        lo[1],lo[2],hi[1],hi[2]=p0[1]-0.3,p0[2]-0.3,p0[1]+0.3,p0[2]+0.3
        p,cov=curve_fit(g2,(X[ok],Y[ok]),s[ok],p0=p0,bounds=(lo,hi),sigma=np.full(ok.sum(),sig),absolute_sigma=True,maxfev=40000)
        e=np.sqrt(np.diag(cov*npb)); r=s[ok]-g2((X[ok],Y[ok]),*p)
        out['g2']=dict(Acomp=float(p[0]),eA=float(e[0]),a=float(p[3]),b=float(p[4]),ea=float(e[3]),eb=float(e[4]),Aext=float(p[7]),aext=float(p[8]),bext=float(p[9]),chi=float(np.std(r)/sig),
          core_dec_maj=float(np.sqrt(max(max(p[3],p[4])**2-bm[0]**2,0))),core_dec_min=float(np.sqrt(max(min(p[3],p[4])**2-bm[1]**2,0))))
    except Exception as ex:
        out['g2']=None
    return out
if __name__=='__main__': main()
