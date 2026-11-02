# Q5h: injection-recovery for (i) deconvolved sizes and (ii) line-vs-continuum centroid offsets, same estimator as the real fits.
# Image plane, Xb4b native (0.246x0.237"). Noise = real line-free sums (HCN/HCO+) and real off-source continuum patches.
# python3 -I scripts/q5h_inject.py -> out/q5h/q5h.log, q5h_size_inj.csv, q5h_size_real.csv, q5h_offsets.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q5lib import *
import pandas as pd, warnings; warnings.filterwarnings('ignore')
from multiprocessing import Pool
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as ndi
OH='/workspace/out/q5h/'; log=None
def pr(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True)
    if log: log.write(s+'\n'); log.flush()
d=load(); dx=d['dx']; dy=d['dy']; pix=float(d['pix'])
mH,_=win_map(d,'H',*VW); mO,_=win_map(d,'O',*VW)
IMG={'HCN':mH,'HCO+':mO}; BM={'HCN':d['beamH'],'HCO+':d['beamO']}; BC={k:beamcov(*BM[k]) for k in BM}
NUL={'HCN':null_maps_scaled(d,'H',len(win_map(d,'H',*VW)[0])*0+int((VW[1]-VW[0])/float(d['dvH']))),
     'HCO+':null_maps_scaled(d,'O',int((VW[1]-VW[0])/float(d['dvO'])))}
# continuum crop on the same grid as the cube cache
hd=fits.getheader('/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')
RA0,DE0=156.96346,-43.90461
cx,cy=[float(a) for a in WCS(hd).celestial.world_to_pixel_values(RA0,DE0)]; hp=int(6./pix); x0=int(cx)-hp; y0=int(cy)-hp
cfull=np.nan_to_num(fits.getdata('/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits').squeeze())
CONT=cfull[y0:y0+dx.shape[0],x0:x0+dx.shape[1]]; BM['cont']=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); BC['cont']=beamcov(*BM['cont'])
IMG['cont']=CONT
# continuum noise patches: annulus 6.5-11" from centre (finite, source free?) -> list of 255x255 not possible; use the full-field noise std and a patch bank of the central-size blocks
Y,X=np.mgrid[0:840,0:840]; rr=np.hypot(X-cx,Y-cy)*pix
ann=(rr>7)&(rr<11.5)&(cfull!=0)
pr('cont beam',BM['cont'],'crop peaks',CONT.max(),'annulus std (mJy/b)',cfull[ann].std()*1e3,'N',ann.sum())
SMIN=0.01
def fitone(img,x,y,rad,bc,bm,amp0):
    mask=np.hypot(dx-x,dy-y)<rad
    z=np.where(mask,ndi.gaussian_filter(img,2),-1e9); iy,ix=np.unravel_index(np.argmax(z),z.shape)
    p=[img[iy,ix]*1.3,dx[iy,ix],dy[iy,ix],0.1,0.1,0.,0,0,0]
    lo=[0,x-rad,y-rad,SMIN,SMIN,-400,-9,-9,-9]; hi=[50,x+rad,y+rad,1.0,1.0,400,9,9,9]
    # scale the amplitude bound/units: normalise image by amp0
    r=fit(img/amp0,dx,dy,mask,np.array(p,float)*np.r_[1/amp0,1,1,1,1,1,1/amp0,1/amp0,1/amp0],np.array(lo)*np.r_[1/amp0,1,1,1,1,1,1,1,1]*np.r_[1,1,1,1,1,1,1/amp0,1/amp0,1/amp0],np.array(hi),bc,1,True)
    par=r.x; fa,fb=decon(par); return par[1],par[2],fa,fb,np.sqrt(fa*fb)
def inject(args):
    line,fw,amp,rad,seed=args; rng=np.random.default_rng(seed)
    bc=BC[line]
    # noise map: line-free sum, random block; for the continuum: random shift of a source-free patch of the annulus is not available -> white-ish smoothed noise with same std
    if line=='cont':
        nm=cont_noise(rng)
    else:
        nm=NUL[line][rng.integers(len(NUL[line]))]
    while True:
        x=rng.uniform(-3.5,3.5); y=rng.uniform(-3.5,3.5)
        if min(np.hypot(x-0,y-2.6),np.hypot(x+.15,y+2.6),np.hypot(x-.3,y-.15))>1.3: break
    # model() peak = F*sqrt(det B/det C): choose F so that the beam-convolved peak equals amp
    C=cov(fw/2.3548,fw/2.3548,0)+bc; F=amp*np.sqrt(np.linalg.det(C)/np.linalg.det(bc))
    src=model(np.array([F,x,y,fw/2.3548,fw/2.3548,0.]),dx,dy,bc,1,False)
    img=nm+src
    try: xf,yf,fa,fb,fg=fitone(img,x,y,rad,bc,BM[line],amp)
    except Exception: return (np.nan,)*5
    return (xf-x,yf-y,fa,fb,fg)
def cont_noise(rng):
    # generate correlated noise with the real off-source rms and the clean-beam autocorrelation: filter white noise by the beam, rescale to the annulus std
    sig=cfull[ann].std(); w=rng.normal(size=dx.shape)
    k=beamcov(*BM['cont']); sx=np.sqrt(k[0,0])/pix; sy=np.sqrt(k[1,1])/pix
    # elliptical gaussian filter through a rotated frame is overkill: use circular filter of the geometric-mean sigma
    s=np.sqrt(np.sqrt(np.linalg.det(k)))/pix; g=ndi.gaussian_filter(w,s/np.sqrt(2)); return g/g.std()*sig
if __name__=='__main__':
    log=open(OH+'q5h.log','w')
    pr('null maps HCN',len(NUL['HCN']),'HCO+',len(NUL['HCO+']))
    nstd={k:np.mean([m[np.hypot(dx,dy)<4].std() for m in NUL[k]]) for k in NUL}; pr('null m0 rms (Jy/b km/s) within 4":',nstd,' cont (Jy/b):',cfull[ann].std())
    NT=150; rows=[]
    sizes=[0.0,0.1,0.2,0.3,0.4,0.5,0.7,1.0]
    jobs=[]
    amps={'HCN':[0.32,0.85],'HCO+':[0.35,0.85]}
    with Pool(16) as pool:
        for line in ('HCN','HCO+'):
            for amp in amps[line]:
                for rad in (0.35,0.5):
                    for fw in sizes:
                        res=pool.map(inject,[(line,fw,amp,rad,i+1000*int(fw*10)+7) for i in range(NT)])
                        a=np.array(res); ok=np.isfinite(a[:,0])
                        for q,nm in ((0,'dx'),(1,'dy')): pass
                        rows.append(dict(line=line,amp=amp,rad=rad,fwhm_true=fw,n=ok.sum(),sd_x=a[ok,0].std(),sd_y=a[ok,1].std(),bias_r=np.median(np.hypot(a[ok,0],a[ok,1])),
                            g16=np.percentile(a[ok,4],2.5),g50=np.median(a[ok,4]),g84=np.percentile(a[ok,4],97.5),gstd=a[ok,4].std(),frac_floor=np.mean(a[ok,3]<0.03)))
                        pr(rows[-1])
    pd.DataFrame(rows).to_csv(OH+'q5h_size_inj.csv',index=False)
