"""Q56 r2: systematic error budget of R32/10 (HCN, HCO+) at N and S, and uv-matching test.
Maps: 1-0 (X961 12m+ACA, native 1.77") and 3-2 (uv-complete 12m+7m 'comb' cubes, circular 2.5"/3.5") -> RJ Tb cubes at circular theta -> moment-0 maps (K km/s) with per-pixel
baseline, resampled (bilinear) to one common 0.1" grid centred between N and S.  Variants: window (+-100/150/350), baseline (const|linear|none-free-range), position (nominal pixel, +-0.25 theta
jackknife grid, HCO+-3-2 peak within 0.5"), estimator (pixel vs aperture-mean r=0.5 theta, r=theta [Q34 style]), uv filter (image-plane high-pass of all maps at uvmin = 0, 7.7, 13 klambda).
Output out/q56/s344.txt/json.  usage: python3 -I s344_q56_systematics.py"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
import s320_q52_chem_ratios as m
from s180_q29_load import to_tb_conv, CEN as CEN0
TBC={}; C=299792.458; O='/workspace/out/q34/'
F0={'HCN':265.88618,'HCO+':267.55762}; TAG={'HCN':'HCN32','HCO+':'HCOp32'}
GR=0.1; HALF=9.0; NG=int(2*HALF/GR)+1; CEN=m.REGIONS['N'][:2]
cen=(0.5*(m.REGIONS['N'][0]+m.REGIONS['S'][0]),0.5*(m.REGIONS['N'][1]+m.REGIONS['S'][1]))
cosd=np.cos(np.deg2rad(cen[1]))
gx,gy=np.meshgrid(np.arange(NG)*GR-HALF,np.arange(NG)*GR-HALF)    # offsets arcsec E (x), N (y)
ra=cen[0]-gx/3600/cosd; de=cen[1]+gy/3600   # RA increases to the east = -x
def kw(v,dv,lo,hi): return np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv
def mom(Tb,v,dv,w,base,free):
    Tb=np.nan_to_num(Tb)
    if base=='const': b=np.median(Tb[free],axis=0)[None]
    elif base=='lin':
        x=v[free]; X=np.vstack([np.ones_like(x),x]).T; coef=np.linalg.lstsq(X,Tb[free].reshape(free.sum(),-1),rcond=None)[0]
        b=(np.outer(v,coef[1])+coef[0]).reshape(Tb.shape)
    else: b=0
    return np.tensordot(kw(v,dv,-w,w),Tb-b,axes=(0,0))*dv
def resample(M,wc):
    px,py=wc.world_to_pixel_values(ra,de); return map_coordinates(M,[py,px],order=1,mode='nearest')
def maps(J,ln,th,pos_v0,w,base):
    """returns moment-0 map on common grid (K km/s) for the two nuclei velocities: dict rn->map"""
    out={}
    if J=='32':
        h=fits.open(O+f'{TAG[ln]}_{str(th).replace(".","p")}_comb.pbcor.fits'); d=h[0].data.astype(float); hd=h[0].header; wc=WCS(hd).celestial; f0=F0[ln]; bm=hd['BMAJ']*3600
        nu=(hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3'])/1e9
        Tb=np.nan_to_num(d)*1.222e6/(f0**2*bm**2)
        for rn,v0 in pos_v0.items():
            v=C*(1-nu/f0)-v0; free=abs(v)>380; dv=abs(np.mean(np.diff(v)))
            out[rn]=resample(mom(Tb,v,dv,w,'const' if base!='none' else 'none',free),wc)
    else:
        path,f0=m.LINES['10'][ln]; cut,nu,pix,bm,pxy=m.getcube(path,th)
        Tb=TBC.setdefault((path,th),to_tb_conv(cut,nu,pix,bm,th,float(np.mean(nu)))); h=fits.open(path,memmap=True); wc0=WCS(h[0].header).celestial
        # cutout origin
        cx,cy=[float(t) for t in wc0.world_to_pixel_values(*CEN0)]; n=int(max(8.0,3*th)/pix); x0=max(int(round(cx))-n,0); y0=max(int(round(cy))-n,0)
        wc=wc0.deepcopy(); wc.wcs.crpix=[wc.wcs.crpix[0]-x0,wc.wcs.crpix[1]-y0]
        oth={'HCN':'HCO+','HCO+':'HCN'}[ln]; g0=m.LINES['10'][oth][1]
        for rn,v0 in pos_v0.items():
            v=m.vel(nu,f0,v0); vg=m.vel(nu,g0,v0); o=np.argsort(v); vv=v[o]; T=Tb[o]
            free=(abs(vv)>600)&(abs(vv)<1200)&~((vg[o]>-650)&(vg[o]<650)); dv=abs(np.mean(np.diff(vv)))
            out[rn]=resample(mom(T,vv,dv,w,base,free),wc)
    return out
def uvfilter(M,uvmin):
    if uvmin==0: return M
    # apodised image-plane high pass: remove |uv|<uvmin (klambda) via a Gaussian-edged mask; map padded x2
    n=NG; P=np.pad(M-0,(n//2,n//2),mode='edge'); F=np.fft.fft2(P); u=np.fft.fftfreq(P.shape[0],d=GR/206265.)/1e3
    U=np.hypot(*np.meshgrid(u,u)); H=1-np.exp(-0.5*(U/(uvmin/2.5))**2)*(U<uvmin*1.5)  # soft edge, ~full removal below uvmin/2
    return np.fft.ifft2(F*H).real[n//2:n//2+n,n//2:n//2+n]
def samp(M,x,y,rad,th):
    if rad==0: return float(map_coordinates(M,[[(y+HALF)/GR],[(x+HALF)/GR]],order=1)[0])
    msk=np.hypot(gx-x,gy-y)<=rad; return float(M[msk].mean())
def offs(rn):
    p=m.REGIONS[rn]; return (-(p[0]-cen[0])*3600*cosd,(p[1]-cen[1])*3600)
if __name__=='__main__':
    pos_v0={'N':2785.,'S':2750.}; res={}; out=[]
    def P(t): print(t,flush=True); out.append(t)
    for th in (2.5,3.5):
      for w in (100,150,350):
        for base in ('lin','const'):
            M={(J,ln):maps(J,ln,th,pos_v0,w,base) for J in ('10','32') for ln in ('HCN','HCO+')}
            for uv in (0,7.7,13.0):
                MF={k:{rn:uvfilter(a,uv) for rn,a in v.items()} for k,v in M.items()}
                for rn in 'NS':
                    x0,y0=offs(rn)
                    # HCO+ 3-2 peak within 0.5"
                    mm=np.where(np.hypot(gx-x0,gy-y0)<0.5,MF[('32','HCO+')][rn],-np.inf); j=np.argmax(mm); xp,yp=gx.flat[j],gy.flat[j]
                    for lab,(dx,dy,rad) in {'nom':(0,0,0),'peak':(xp-x0,yp-y0,0),'ap0.5':(0,0,0.5*th),'ap1':(0,0,th)}.items():
                        for ln in ('HCN','HCO+'):
                            i10=samp(MF[('10',ln)][rn],x0+dx,y0+dy,rad,th); i32=samp(MF[('32',ln)][rn],x0+dx,y0+dy,rad,th)
                            res[f'{th}|{w}|{base}|{uv}|{rn}|{lab}|{ln}']=[i10,i32]
                    # jackknife grid +-0.25 theta around nominal
                    for ix,iy in [(a,b) for a in (-1,0,1) for b in (-1,0,1)]:
                        for ln in ('HCN','HCO+'):
                            dx,dy=0.25*th*ix,0.25*th*iy
                            res[f'{th}|{w}|{base}|{uv}|{rn}|jk{ix}{iy}|{ln}']=[samp(MF[('10',ln)][rn],x0+dx,y0+dy,0,th),samp(MF[('32',ln)][rn],x0+dx,y0+dy,0,th)]
    json.dump(res,open('/workspace/out/q56/s344.json','w'),default=float)
    print('saved',len(res))
