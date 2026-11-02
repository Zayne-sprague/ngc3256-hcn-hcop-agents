# Continuum source inventory: peaks, 2D-Gaussian fits (deconvolved), line contamination mfs vs line-free cont
import sys,warnings,json,glob; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from lineset import RA0,DE0
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
import pandas as pd
COS=np.cos(np.deg2rad(DE0))
def img(path):
    h=fits.open(path); hd=h[0].header; im=np.squeeze(h[0].data).astype(float); return hd,im,WCS(hd).celestial
def offs(w,x,y):
    ra,de=w.pixel_to_world_values(x,y); return (ra-RA0)*COS*3600,(de-DE0)*3600
def rms_est(hd,im,w,rin,rout):
    ny,nx=im.shape; Y,X=np.mgrid[:ny,:nx]; dx,dy=offs(w,X,Y); r=np.hypot(dx,dy); m=(r>rin)&(r<rout)&np.isfinite(im)
    v=im[m]; 
    for _ in range(5):
        s=1.4826*np.median(np.abs(v-np.median(v))); v=v[np.abs(v-np.median(v))<4*s]
    return s,np.median(v)
def g2(xy,a,x0,y0,sx,sy,th,c):
    x,y=xy; ct,st=np.cos(th),np.sin(th); xr=(x-x0)*ct+(y-y0)*st; yr=-(x-x0)*st+(y-y0)*ct
    return a*np.exp(-0.5*((xr/sx)**2+(yr/sy)**2))+c
def fitsrc(hd,im,w,ra,de,box_arcsec):
    pix=abs(hd['CDELT2'])*3600; x0,y0=[float(a) for a in w.world_to_pixel_values(ra,de)]
    r=int(box_arcsec/pix); xi,yi=int(round(x0)),int(round(y0))
    ys,xs=np.mgrid[max(yi-r,0):yi+r+1,max(xi-r,0):xi+r+1]; z=im[ys,xs]
    if not np.all(np.isfinite(z)) or z.size<20: return None
    # start: peak
    iy,ix=np.unravel_index(np.argmax(z),z.shape)
    bmaj,bmin=hd['BMAJ']*3600/pix,hd['BMIN']*3600/pix
    p0=[z.max(),xs[iy,ix],ys[iy,ix],bmaj/2.355,bmin/2.355,0.,0.]
    try:
        pp,cv=curve_fit(g2,(xs.ravel(),ys.ravel()),z.ravel(),p0=p0,maxfev=20000)
    except Exception: return None
    e=np.sqrt(np.diag(cv)); a,gx,gy,sx,sy,th,c=pp
    fx,fy=abs(sx)*2.355*pix,abs(sy)*2.355*pix; a_,b_=max(fx,fy),min(fx,fy)
    bm,bn=hd['BMAJ']*3600,hd['BMIN']*3600
    S=a*(fx*fy)/(bm*bn)      # Jy (integrated), peak Jy/beam * (source gauss area / beam area), with convolved sizes
    dec=lambda f,b: np.sqrt(f**2-b**2) if f>b else 0.
    ra_f,de_f=[float(q) for q in w.pixel_to_world_values(gx,gy)]
    return dict(peak_fit=a,e_peak=e[0],ra=ra_f,dec=de_f,fwhm_maj_conv=a_,fwhm_min_conv=b_,S_int=S,decon_maj=dec(a_,bm),decon_min=dec(b_,bn),offset=c)
SRC={'N':(0.0,2.55),'S':(-0.15,-2.55),'NEclump':(6.75,3.90),'Eclump':(5.25,-2.85),'SEcont':(8.8,-5.1),'Bridge':(0.3,0.15)}
if __name__=='__main__':
    rows=[]
    files=sorted(glob.glob(D+'*cont*pbcor*.fits')+glob.glob(D+'*mfs.I.pbcor*.fits'))
    for p in files:
        key=os.path.basename(p).replace('member.uid___A001_X1273_','').replace('.ngc_3256_sci','').replace('_ngc_3256_sci','_tap').replace('.fits','')
        ous=key.split('.')[0].split('_')[0]
        hd,im,w=img(p); pix=abs(hd['CDELT2'])*3600
        sg,med=rms_est(hd,im,w,{'Xb55':20,'Xb47':12,'Xb4f':14,'Xb4d':10,'Xb4b':10,'Xb4d_tap':10}[ous if ous in('Xb55','Xb47','Xb4f','Xb4d','Xb4b') else 'Xb4d_tap'],
                       {'Xb55':45,'Xb47':25,'Xb4f':30,'Xb4d':18,'Xb4b':18}.get(ous,18))
        bm,bn=hd['BMAJ']*3600,hd['BMIN']*3600
        for sn,(dx,dy) in SRC.items():
            ra=RA0+dx/3600/COS; de=DE0+dy/3600
            x0,y0=[float(a) for a in w.world_to_pixel_values(ra,de)]
            r=int(max(1.0,0.4*bm)/pix); xi,yi=int(round(x0)),int(round(y0))
            if not(r<xi<im.shape[1]-r and r<yi<im.shape[0]-r): continue
            st=im[yi-r:yi+r+1,xi-r:xi+r+1]; pk=np.nanmax(st); iy,ix=np.unravel_index(np.nanargmax(st),st.shape)
            row=dict(file=key,band=hd.get('CRVAL3',0),src=sn,peak_Jybeam=pk,rms=sg,snr=pk/sg,beam=f'{bm:.2f}x{bn:.2f}',nu_GHz=hd['CRVAL3']/1e9 if hd.get('CTYPE3','').startswith('FREQ') else np.nan)
            if pk/sg>5:
                ft=fitsrc(hd,im,w,*[float(q) for q in w.pixel_to_world_values(xi-r+ix,yi-r+iy)],box_arcsec=max(1.5,1.2*bm))
                if ft: row.update(ft)
            rows.append(row)
    df=pd.DataFrame(rows); df.to_csv('/workspace/out/continuum_table.csv',index=False)
    pd.set_option('display.width',250,'display.max_rows',900,'display.max_columns',30)
    print(df[['file','src','peak_Jybeam','rms','snr','beam','peak_fit','S_int','decon_maj','decon_min']].round(5).to_string())
