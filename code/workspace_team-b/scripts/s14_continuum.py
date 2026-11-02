# Continuum: (A) blind >5 sigma source finder + 2D Gaussian fits with error budget on the line-free pbcor cont images;
# (B) all OUS convolved to common beams; matched peak and aperture flux densities, spectral indices with flux-scale errors.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
from scipy import ndimage, optimize
from scipy.signal import fftconvolve
CONT={'Xb55':'*Xb55*16_18_20_22.cont.I.pbcor.fits','Xb47':'*Xb47*16_18_20_22.cont.I.pbcor.fits','Xb4f':'*Xb4f*16_18_20_22.cont.I.pbcor.fits',
      'Xb4d':'*Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits','Xb4b':'*Xb4b*25_27_29_31.cont.I.pbcor.fits'}
CAL={'Xb55':0.10,'Xb47':0.10,'Xb4f':0.05,'Xb4d':0.05,'Xb4b':0.05}
PH=(156.96346,-43.90461); FW={'Xb55':33.0,'Xb47':16.6,'Xb4f':22.0,'Xb4d':22.0,'Xb4b':22.0}
log=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); log.append(s)
def local_rms(d,blk):
    ny,nx=d.shape; out=np.full(d.shape,np.nan)
    ys=range(0,ny,blk); xs=range(0,nx,blk); R=np.full((len(ys),len(xs)),np.nan)
    for i,y in enumerate(ys):
        for j,x in enumerate(xs):
            c=d[max(0,y-blk//2):y+blk+blk//2,max(0,x-blk//2):x+blk+blk//2]; c=c[np.isfinite(c)]
            if c.size>50:
                m=np.median(c); s=1.4826*np.median(abs(c-m)); c2=c[abs(c-m)<4*s]; R[i,j]=np.std(c2) if c2.size>30 else s
    R=ndimage.generic_filter(R,lambda a: np.nanmedian(a),size=3,mode='nearest')
    Y,X=np.mgrid[:ny,:nx]; from scipy.interpolate import RegularGridInterpolator
    f=RegularGridInterpolator((np.array(list(ys))+blk/2,np.array(list(xs))+blk/2),np.nan_to_num(R,nan=np.nanmedian(R)),bounds_error=False,fill_value=None)
    return f(np.stack([Y.ravel(),X.ravel()],-1)).reshape(d.shape)
def g2(p,xx,yy):
    A,x0,y0,sx,sy,th,c=p; ct,st=np.cos(th),np.sin(th); u=(xx-x0)*ct+(yy-y0)*st; v=-(xx-x0)*st+(yy-y0)*ct
    return A*np.exp(-0.5*((u/sx)**2+(v/sy)**2))+c
def fit_source(d,x,y,bpix,sigma):
    r=int(max(4,2.2*bpix)); x0,x1,y0,y1=int(x)-r,int(x)+r+1,int(y)-r,int(y)+r+1
    if x0<0 or y0<0 or x1>d.shape[1] or y1>d.shape[0]: return None
    cut=d[y0:y1,x0:x1]; yy,xx=np.mgrid[y0:y1,x0:x1]; s0=bpix/2.355
    p0=[cut.max(),x,y,s0,s0,0,0]
    try:
        p,cv=optimize.curve_fit(lambda xy,*p: g2(p,xy[0],xy[1]).ravel(),(xx,yy),cut.ravel(),p0=p0,sigma=np.full(cut.size,sigma),absolute_sigma=True,maxfev=5000)
        return p,np.sqrt(np.diag(cv))
    except Exception: return None
rows=[]; imgs={}
for ous,pat in CONT.items():
    fn=files(pat)[0]; d,h=load2d(fn); w=WCS(h).celestial; ps=abs(h['CDELT1'])*3600; bmaj,bmin,bpa=h['BMAJ']*3600,h['BMIN']*3600,h['BPA']
    nu=h['CRVAL3']/1e9; imgs[ous]=(d,h,w,ps,bmaj,bmin,bpa,nu)
    blk=int(max(12,6*bmaj/ps)); rms=local_rms(d,blk); sn=d/rms
    bpix=bmaj/ps    # FWHM in pixels (major)
    mx=ndimage.maximum_filter(np.nan_to_num(d,nan=-1),size=int(max(5,1.5*bmaj/ps)))
    pk=np.argwhere((d==mx)&(sn>5)); neg=np.argwhere((d==ndimage.minimum_filter(np.nan_to_num(d,nan=1),size=int(max(5,1.5*bmaj/ps))))&(sn<-5))
    P(f'## {ous} cont {nu:.1f} GHz beam {bmaj:.2f}x{bmin:.2f}: median local rms {np.nanmedian(rms)*1e3:.3f} mJy/beam; >5sigma maxima {len(pk)}; < -5sigma minima {len(neg)} (false-positive check)')
    for (iy,ix) in sorted(pk.tolist(),key=lambda p:-sn[p[0],p[1]]):
        f=fit_source(d,ix,iy,bpix,rms[iy,ix]); ra,dec=w.pixel_to_world_values(ix,iy)
        r=dict(ous=ous,nu_GHz=nu,ra0=float(ra),dec0=float(dec),peak_mJy=d[iy,ix]*1e3,local_rms_mJy=rms[iy,ix]*1e3,snr=sn[iy,ix])
        if f is not None:
            p,e=f; A,x0,y0,sx,sy,th,c=p; raf,decf=w.pixel_to_world_values(x0,y0)
            fx,fy=sorted([abs(sx)*2.355*ps,abs(sy)*2.355*ps],reverse=True)
            Sint=A*2*np.pi*abs(sx*sy)*ps**2/(np.pi*bmaj*bmin/(4*np.log(2)))*1e3   # mJy: A[Jy/beam]*Omega_fit/Omega_beam
            # deconvolved size (FWHM, arcsec): simple quadrature (major/minor separately, ignoring PA mismatch)
            dmaj=np.sqrt(max(fx**2-bmaj**2,0)); dmin=np.sqrt(max(fy**2-bmin**2,0))
            snr=A/rms[iy,ix]
            sig_pos=0.42*np.sqrt(bmaj*bmin)/snr    # Condon (1997) relative position error per axis, arcsec
            r.update(ra=float(raf),dec=float(decf),fit_peak_mJy=A*1e3,fit_peak_err_mJy=e[0]*1e3,fit_maj=fx,fit_min=fy,decon_maj=dmaj,decon_min=dmin,S_int_mJy=Sint,
                     S_int_err_mJy=Sint*np.hypot(e[0]/A,0.1)  ,sig_pos_arcsec=sig_pos,sig_pos_abs_arcsec=np.hypot(sig_pos,0.1*np.sqrt(bmaj*bmin)))
        rows.append(r)
cat=pd.DataFrame(rows)
cat=cat[cat.ra.notna()].copy()
cat['fit_flag']=np.where((cat.S_int_mJy<=8*cat.fit_peak_mJy)&(cat.decon_maj<3.0)&(cat.fit_peak_mJy>0.5*cat.peak_mJy),'ok','extended/blended: fit unreliable')
REG={'N':(156.96346,-43.90388),'S':(156.96340,-43.90533),'E1':(156.96608,-43.90350),'E2':(156.96548,-43.90540),'SE':(156.96685,-43.90635),'W':(156.96114,-43.90590)}
def assoc(r):
    best=None
    for k,(ra,dec) in REG.items():
        d=np.hypot((r.ra-ra)*np.cos(np.deg2rad(dec)),r.dec-dec)*3600
        if best is None or d<best[1]: best=(k,d)
    return best[0] if best[1]<(2.0 if best[0]=='N' else 1.0) else ('N-complex' if np.hypot((r.ra-REG['N'][0])*np.cos(np.deg2rad(REG['N'][1])),r.dec-REG['N'][1])*3600<2.5 else 'other')
cat['association']=cat.apply(assoc,axis=1)
cat.to_csv('/workspace/out/s14_source_catalogue.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',200)
for ous in CONT:
    c=cat[cat.ous==ous]
    if len(c): P(c[['ous','ra','dec','snr','peak_mJy','fit_peak_mJy','S_int_mJy','decon_maj','decon_min','sig_pos_arcsec']].round(5 if False else 4).head(25).to_string())
# ---- (B) common beam
def cov(a,b,pa):
    th=np.deg2rad(pa); e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
def convolve_to(d,ps,b0,tgt):
    C2=cov(*tgt)-cov(*b0); w_,U=np.linalg.eigh(C2)
    if w_.min()<-1e-3*w_.max(): return None
    C2=U@np.diag(np.clip(w_,1e-6,None))@U.T
    n=int(8*np.sqrt(w_.max())/ps)|1; n=min(n,min(d.shape)-1|1)
    ky,kx=np.mgrid[:n,:n]-n//2; q=np.stack([(-kx*ps).ravel(),(ky*ps).ravel()]); Ci=np.linalg.inv(C2)
    ker=np.exp(-0.5*np.einsum('ij,ij->j',q,Ci@q)).reshape(n,n); ker/=ker.sum()
    return fftconvolve(np.nan_to_num(d),ker,mode='same')*(tgt[0]*tgt[1])/(b0[0]*b0[1])
SRC={'N':NUC['N'],'S':NUC['S'],'E1':(156.96608,-43.90350)}
res=[]
for tag,tgt in {'B5beam_8.95x5.98':(8.95,5.98,-88.8),'B6lowres_5.84x3.52':(5.84,3.52,82.4)}.items():
    P(f'## common beam {tag}')
    for ous,(d,h,w,ps,bmaj,bmin,bpa,nu) in imgs.items():
        c=d.copy() if (abs(bmaj-tgt[0])<0.03*tgt[0] and abs(bmin-tgt[1])<0.03*tgt[1]) else convolve_to(d,ps,(bmaj,bmin,bpa),tgt)
        if c is None: P(f'  {ous}: native beam not smaller than target, skipped'); continue
        # local noise: robust rms of the convolved map in an off-source annulus (r 12-20")
        Y,X=np.mgrid[:d.shape[0],:d.shape[1]]; mid=w.world_to_pixel_values(*[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)])
        rr=np.hypot(X-mid[0],Y-mid[1])*ps; ann=(rr>9)&(rr<16)&np.isfinite(d); rm=1.4826*np.median(abs(c[ann]-np.median(c[ann])))
        row=dict(beam=tag,ous=ous,nu_GHz=nu,rms_mJy=rm*1e3)
        for k,(ra,dec) in SRC.items():
            x,y=w.world_to_pixel_values(ra,dec); m=np.hypot(X-x,Y-y)*ps<=1.0
            row[k+'_peak_mJy']=np.max(np.where(m,c,-1e9))*1e3
        # integrated flux in a matched aperture of radius 6" around the N-S midpoint from the ORIGINAL image
        bomega=np.pi*bmaj*bmin/(4*np.log(2))
        for R in (3.0,6.0,10.0):
            for k in ('N','S'):
                x,y=w.world_to_pixel_values(*NUC[k]); ap=np.hypot(X-x,Y-y)*ps<=R
                row[f'{k}_S_r{R:.0f}_mJy']=np.nansum(d[ap])*ps**2/bomega*1e3
        for R in (8.0,12.0):
            ap=(rr<=R)&np.isfinite(d); row[f'tot_r{R:.0f}_mJy']=np.nansum(c[ap])*ps**2/(np.pi*tgt[0]*tgt[1]/(4*np.log(2)))*1e3
        res.append(row)
        P('  ',{k:(round(v,3) if isinstance(v,float) else v) for k,v in row.items() if k not in('beam',)})
R=pd.DataFrame(res); R.to_csv('/workspace/out/s14_common_beam.csv',index=False)
# spectral indices at the B5 beam and B6low beam (peak flux density in matched beam), errors: noise + flux scale (independent per OUS)
P('## spectral indices (matched beam, peak flux densities; error = noise + flux-scale in quadrature)')
al=[]
for tag in R.beam.unique():
    r=R[R.beam==tag].set_index('ous')
    for k in ('N','S','E1'):
        pairs=[('Xb55','Xb4f'),('Xb55','Xb4d'),('Xb55','Xb4b'),('Xb4f','Xb47'),('Xb4d','Xb47'),('Xb4b','Xb47'),('Xb55','Xb47')]
        for a,b in pairs:
            if a in r.index and b in r.index:
                sa,sb=r.loc[a,k+'_peak_mJy'],r.loc[b,k+'_peak_mJy']; ea=np.hypot(r.loc[a,'rms_mJy'],CAL[a]*sa); eb=np.hypot(r.loc[b,'rms_mJy'],CAL[b]*sb)
                if sa>0 and sb>0:
                    L=np.log(r.loc[b,'nu_GHz']/r.loc[a,'nu_GHz']); alpha=np.log(sb/sa)/L; ea_=np.hypot(ea/sa,eb/sb)/abs(L)
                    al.append(dict(beam=tag,src=k,a=a,b=b,nu_a=r.loc[a,'nu_GHz'],nu_b=r.loc[b,'nu_GHz'],S_a=sa,S_b=sb,alpha=alpha,alpha_err=ea_))
for tag in R.beam.unique():
    r=R[R.beam==tag].set_index('ous')
    for Rr in (8,12):
        for a,b in [('Xb55','Xb4f'),('Xb55','Xb4d'),('Xb55','Xb4b'),('Xb4f','Xb47'),('Xb4d','Xb47'),('Xb4b','Xb47'),('Xb55','Xb47')]:
            if a in r.index and b in r.index:
                sa,sb=r.loc[a,f'tot_r{Rr}_mJy'],r.loc[b,f'tot_r{Rr}_mJy']; ea=CAL[a]*sa; eb=CAL[b]*sb
                L=np.log(r.loc[b,'nu_GHz']/r.loc[a,'nu_GHz']); al.append(dict(beam=tag,src=f'total r<{Rr}"',a=a,b=b,nu_a=r.loc[a,'nu_GHz'],nu_b=r.loc[b,'nu_GHz'],S_a=sa,S_b=sb,alpha=np.log(sb/sa)/L,alpha_err=np.hypot(ea/sa,eb/sb)/abs(L)))
A=pd.DataFrame(al); A.to_csv('/workspace/out/s14_alpha.csv',index=False); P(A.round(2).to_string())
open('/workspace/out/s14_continuum.txt','w').write('\n'.join(log))
