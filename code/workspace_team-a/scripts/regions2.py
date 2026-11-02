# Region definitions from moment-0 maps of several tracers and two resolutions + astrometric error budget.
# Maps: Xb4d and Xb4b: HCO+ 3-2 (spw25), HCN 3-2 (spw27), CH3OH 2(0)-1(-1) E (spw29, weak tracer), 1.1mm continuum (cont image).  Peaks (S/N>thr, r<12") are fitted with a circular 2-D Gaussian
# (free centre, FWHM, amplitude, offset); centroid errors from the fit (noise = per-pixel line-free rms) are compared with theta/(2 S/N); pointing/phase-transfer floor taken as 10% of the beam.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
import pandas as pd
RA0,DE0=156.96346,-43.90461
def mom(ous,spw,f0,vc=2780.,vw=350.):
    h=fits.open(f(ous,spw),memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C; dv=abs(np.median(np.diff(v)))
    sel=np.where((v>=vc-vw)&(v<=vc+vw))[0]; cube=np.array(h[0].data[0,sel.min():sel.max()+1],float)
    # noise: line-free channels = the first 20 channels of the cube (far from every line)
    off=np.array(h[0].data[0,0:20],float); sig=np.nanstd(off,axis=0)
    return hd,np.nansum(cube,axis=0)*dv,sig*dv*np.sqrt(len(sel))
def cont(ous):
    p=D+f'member.uid___A001_X1273_{ous}.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'; h=fits.open(p); hd=h[0].header; im=np.array(h[0].data).squeeze()
    # noise: robust rms of pixels at r>14"
    w=WCS(hd).celestial; Y,X=np.mgrid[:im.shape[0],:im.shape[1]]; ra,de=w.pixel_to_world_values(X,Y); r=np.hypot((ra-RA0)*np.cos(np.deg2rad(DE0))*3600,(de-DE0)*3600)
    far=(r>14)&(r<20); sg=1.4826*np.nanmedian(np.abs(im[far]-np.nanmedian(im[far])))
    return hd,im*1e3,np.full(im.shape,sg*1e3)
def g2(xy,a,xc,yc,fw,b):
    return a*np.exp(-4*np.log(2)*((xy[0]-xc)**2+(xy[1]-yc)**2)/fw**2)+b
def analyse(hd,m0,e0,thr,label,unit):
    w=WCS(hd).celestial; ny,nx=m0.shape; Y,X=np.mgrid[:ny,:nx]; ra,de=w.pixel_to_world_values(X,Y)
    dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600; r=np.hypot(dx,dy); pix=abs(hd['CDELT2'])*3600; bm=np.sqrt(hd['BMAJ']*hd['BMIN'])*3600
    e1=ndi.uniform_filter(np.where(np.isfinite(e0),e0,np.nanmedian(e0)),max(3,int(bm/pix)|1)); sm=ndi.gaussian_filter(np.nan_to_num(m0),0.5*bm/pix/2.355)
    snr=sm/e1; ms=int(max(1.2*bm,0.5)/pix)|1
    lm=(sm==ndi.maximum_filter(sm,size=ms))&(snr>thr)&(r<12)
    rows=[]
    for y,x in sorted(zip(*np.where(lm)),key=lambda t:-snr[t]):
        box=int(1.6*bm/pix); ys=slice(max(0,y-box),y+box+1); xs=slice(max(0,x-box),x+box+1)
        Z=m0[ys,xs]; E=e0[ys,xs]; ok=np.isfinite(Z)&np.isfinite(E)&(E>0)
        try:
            p,cv=curve_fit(g2,(dx[ys,xs][ok],dy[ys,xs][ok]),Z[ok],p0=[m0[y,x],dx[y,x],dy[y,x],bm*1.2,0],sigma=E[ok],bounds=([0,dx[y,x]-bm,dy[y,x]-bm,0.5*bm,-np.inf],[np.inf,dx[y,x]+bm,dy[y,x]+bm,6*bm,np.inf]),maxfev=20000)
            chi=np.sum(((Z[ok]-g2((dx[ys,xs][ok],dy[ys,xs][ok]),*p))/E[ok])**2)/(ok.sum()-5); e=np.sqrt(np.diag(cv))*np.sqrt(max(chi,1))
        except Exception as ex: continue
        pk=p[0]/np.nanmedian(E[ok]); fl=np.sqrt(0.1**2*bm**2)
        rows.append(dict(map=label,dRA=p[1],dDec=p[2],e_fit=float(np.hypot(e[1],e[2])/np.sqrt(2)),sigma_theory=bm/(2*pk),floor_0p1beam=0.1*bm,peak=p[0],unit=unit,SN=float(pk),fwhm_fit=p[3],beam=bm,
                         fwhm_deconv=float(np.sqrt(max(p[3]**2-bm**2,0))),r=float(np.hypot(p[1],p[2]))))
    return rows
if __name__=='__main__':
    allr=[]
    for ous in ('Xb4d','Xb4b'):
        for spw,nm,f0,vw,thr in (('spw25','HCO+ 3-2',267.557626,350,8),('spw27','HCN 3-2',265.886434,350,8),('spw29','CH3OH 2(0)-1(-1)E',254.015377,130,4)):
            hd,m0,e0=mom(ous,spw,f0,vw=vw); allr+=analyse(hd,m0,e0,thr,f'{ous} {nm}','Jy/b km/s')
            np.savez(f'/workspace/out/m0b_{ous}_{spw}.npz',m0=m0,e0=e0)
            print(ous,nm,'done',flush=True)
        hd,im,e=cont(ous); allr+=analyse(hd,im,e,8,f'{ous} 1.1mm cont','mJy/b')
    df=pd.DataFrame(allr); df.round(3).to_csv('/workspace/out/peaks_astrometry.csv',index=False)
    pd.set_option('display.width',250,'display.max_rows',300); print(df.round(2).to_string())
