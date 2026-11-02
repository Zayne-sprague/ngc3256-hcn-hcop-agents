# Q5d: beam-deconvolved 2D Gaussian sizes + centroids of N, S, Bridge in Xb4b HCN, HCO+ 3-2 moment-0 (+-350 km/s) and the 1.1 mm continuum; line-continuum offsets.
# python3 -I scripts/q5d_fits.py -> out/q5/q5d.log, q5d_fits.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q5lib import *
from q1_lib import offsets, RA0, DE0
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
d=load(); dx=d['dx']; dy=d['dy']; pix=float(d['pix']); log=open(O+'q5d.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
# continuum, same crop as the cubes (+-6", same pixel grid)
h=fits.open('/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'); hd=h[0].header
w=WCS(hd).celestial; cx,cy=[float(a) for a in w.world_to_pixel_values(RA0,DE0)]; hp=int(6./pix)
x0,y0=max(0,int(cx)-hp),max(0,int(cy)-hp); cont=np.squeeze(h[0].data)[y0:int(cy)+hp+1,x0:int(cx)+hp+1].astype(float)
cdx,cdy=offsets(hd,*cont.shape,x0,y0); cbeam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
pr('continuum crop',cont.shape,'dx range',cdx.min(),cdx.max(),'beam',cbeam,'; cube grid',dx.shape,dx.min(),dx.max())
assert cont.shape==dx.shape and np.allclose(cdx,dx,atol=1e-3)
mH,nH=win_map(d,'H',*VW); mO,nO=win_map(d,'O',*VW)
# empirical thermal noise per pixel (line-free sums) for the weights; continuum noise from an emission-free corner
sH=np.std(null_maps_scaled(d,'H',nH)[0][60:190,60:190]); sO=np.std(null_maps_scaled(d,'O',nO)[0][60:190,60:190])
sC=1.4826*np.median(np.abs(cont[:60,:60]-np.median(cont[:60,:60])))
pr(f'thermal noise per pixel: HCN {sH:.4f}, HCO+ {sO:.4f} Jy/b km/s; continuum {sC*1e6:.0f} uJy/b')
MAPS={'HCN':(mH,sH,(d['beamH'])),'HCO+':(mO,sO,(d['beamO'])),'cont':(cont*1e3,sC*1e3,cbeam)}   # continuum in mJy/beam
SRC={'S':(-0.15,-2.55),'N':(0.0,2.6),'Bridge':(0.3,0.15)}
rows=[]
def run(name,src,radius,ncomp=1,bg=True,tag=''):
    img,sg,bm=MAPS[name]; x,y=SRC[src]; Bc=beamcov(*bm)
    X=dx; Y=dy; mask=np.hypot(dx-x,dy-y)<radius
    # start at local maximum
    z=np.where(mask,ndi_smooth(img),-1e9); iy,ix=np.unravel_index(np.argmax(z),z.shape); pk=img[iy,ix]
    F0=pk*1.3; p=[F0,dx[iy,ix],dy[iy,ix],0.12,0.10,0.]; lo=[0,x-radius,y-radius,0.005,0.005,-400]; hi=[1e3,x+radius,y+radius,2.0,2.0,400]
    for i in range(1,ncomp):
        p+=[0.3*F0,x,y,0.4,0.3,0.]; lo+=[0,x-radius,y-radius,0.005,0.005,-400]; hi+=[1e3,x+radius,y+radius,3.0,3.0,400]
    if bg: p+=[0,0,0]; lo+=[-1e3,-1e3,-1e3]; hi+=[1e3,1e3,1e3]
    r=fit(img,X,Y,mask,np.array(p),np.array(lo),np.array(hi),Bc,ncomp,bg)
    par=r.x; res=img-model(par,X,Y,Bc,ncomp,bg); rms=np.std(res[mask]); 
    # positional/size errors: Condon (1997) for the first component with S/N = peak/rms_res (rms of residual includes real structure -> conservative)
    out=[]
    for ic in range(ncomp):
        F,xc,yc,sa,sb,pa=par[6*ic:6*ic+6]; fa,fb=decon(par[6*ic:6*ic+6]); pa=(pa if sa>=sb else pa+90)%180
        bmaj,bmin=bm[0],bm[1]; cm,cn=np.sqrt(fa**2+bmaj**2),np.sqrt(fb**2+bmin**2)
        pkc=model(np.r_[par[6*ic:6*ic+6]],X,Y,Bc,1,False).max()      # peak of the convolved component (map units)
        k=(cm*cn/(4*bmaj*bmin)); snr=pkc/rms
        rho=lambda a,b: np.sqrt(k*(1+bmaj**2/cm**2)**a*(1+bmin**2/cn**2)**b)*snr   # Condon 1997 eqs 21 (noise = fit residual rms: conservative, includes real structure)
        ra_,rm_,rn_=rho(1.5,1.5),rho(2.5,0.5),rho(0.5,2.5)
        e_pos=cm/2.3548*np.sqrt(2)/ra_; e_maj=np.sqrt(2)/rm_*cm; e_min=np.sqrt(2)/rn_*cn
        # deconvolved-size error: d(fa)=cm/fa*d(cm)
        efa=cm/fa*e_maj if fa>0.03 else np.nan; efb=cn/fb*e_min if fb>0.03 else np.nan
        out.append(dict(src=src,map=name,tag=tag,radius=radius,ncomp=ncomp,comp=ic,flux=F,x=xc,y=yc,fwhm_maj=fa,fwhm_min=fb,pa=pa,peak=pkc,res_rms=rms,snr=snr,e_pos=e_pos,e_maj=efa,e_min=efb,unres=(fa<0.03 or fb<0.03)))
    rows.extend(out); return out
from scipy import ndimage as ndi
def ndi_smooth(a): return ndi.gaussian_filter(a,2)
for name in ('HCN','HCO+','cont'):
    for src,rads in (('S',(0.6,0.8,1.0,1.4)),('Bridge',(0.6,0.8,1.0)),('N',(0.6,0.8,1.0,1.5,2.0))):
        for r_ in rads:
            for bg in (True,):
                try: run(name,src,r_,1,bg)
                except Exception as e: print('fail',name,src,r_,e)
    for src,rads in (('S',(1.4,2.0)),('N',(2.0,2.5))):
        for r_ in rads:
            try: run(name,src,r_,2,True,tag='2comp')
            except Exception as e: print('fail2',name,src,r_,e)
df=pd.DataFrame(rows); df.to_csv(O+'q5d_fits.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
pr(df[['src','map','radius','ncomp','comp','flux','x','y','fwhm_maj','fwhm_min','pa','peak','snr','e_pos','e_maj','e_min','res_rms']].round(3).to_string())
