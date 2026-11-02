# Spatial empirical-noise test for tentative lines: template-matched integrated intensity I_T at every aperture position of the field
# (same template, same velocity window, per-pixel baseline from line-free channels) -> distribution away from the nuclei; compare with the value at N/S.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
from scipy import ndimage
import json
exec(open('/workspace/scripts/s12_line_table.py').read().split("rows=[]\nfor lab,f0,lst,note in LINES:")[0])
TESTS=[('HOC+ 2-1','N','Xb55','spw18'),('HOC+ 2-1','S','Xb55','spw18'),('CH3OH J-1-J0 E blend','N','Xb55','spw22'),('CH3OH J-1-J0 E blend','S','Xb55','spw22'),
       ('SO 5_6-4_5','S','Xb4f','spw22'),('SO 5_6-4_5','S','Xb4d','spw31'),('NS 11/2-9/2 f (253.968; "253.9 bump")','S','Xb4f','spw20'),('NS 11/2-9/2 f (253.968; "253.9 bump")','S','Xb4d','spw29'),
       ('HC3N 28-27','N','Xb4f','spw20'),('HC3N 28-27','N','Xb4d','spw29'),('SiS 14-13','S','Xb4f','spw20'),('CS 7-6','N','Xb47','spw22'),('H41beta RRL','N','Xb55','spw16')]
LD={l[0]:l for l in LINES}
out=[]
def disk(r):
    n=int(np.ceil(r))*2+1; yy,xx=np.mgrid[:n,:n]-n//2; k=(np.hypot(xx,yy)<=r).astype(float); return k/k.sum()
for lab,reg,ous,spw in TESTS:
    f0=LD[lab][1]; T=make_template(ous,reg,'ap'); half=0.5*(T['hi']-T['lo']); vc=0.5*(T['hi']+T['lo'])
    fn=files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]
    h,hd,w,nu=cube_info(fn); v=vel(nu,f0); ps=abs(hd['CDELT1'])*3600; dv=abs(np.median(np.diff(v)))
    cube=np.array(h[0].data[0],dtype=np.float32)
    wc=WCS(hd).celestial; cx,cy=wc.world_to_pixel_values(*[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)])
    R=10.0 if hd['BMAJ']*3600<2 else 30.0; r=int(R/ps); x0,x1=max(0,int(cx-r)),min(hd['NAXIS1'],int(cx+r)); y0,y1=max(0,int(cy-r)),min(hd['NAXIS2'],int(cy+r))
    cube=cube[:,y0:y1,x0:x1]
    free=np.ones(len(v),bool)
    for fk in [f0]+STRONGK.get((ous,spw),[])+EXTRA.get((ous,spw),[]):
        nuk=fk*1e9*(1-vc/C); vk=C*(1-nuk/(f0*1e9)); free&=abs(v-vk)>half+150
    base=np.nanmedian(cube[free],axis=0); cube=cube-base
    bmaj,bmin=hd['BMAJ']*3600,hd['BMIN']*3600; apr=0.5*np.sqrt(bmaj*bmin)/ps; K=disk(apr)
    t=np.interp(v,T['v'],T['t'],left=0,right=0); t=t/t.max()
    A=np.tensordot(t,np.nan_to_num(cube),axes=(0,0))/np.sum(t**2)*np.sum(t)*dv    # per-pixel matched integrated intensity
    A=ndimage.convolve(A,K,mode='constant')
    yy,xx=np.mgrid[:A.shape[0],:A.shape[1]]
    def pos(k):
        px,py=wc.world_to_pixel_values(REGIONS[k][0],REGIONS[k][1]); return px-x0,py-y0
    excl=np.zeros(A.shape,bool)
    for k in ('N','S','Bridge','E1','E2','Sblue','W'):
        px,py=pos(k); excl|=np.hypot(xx-px,yy-py)*ps<1.2*np.sqrt(bmaj*bmin)
    edge=np.hypot(xx-A.shape[1]/2,yy-A.shape[0]/2)*ps<0.9*R
    sel=edge&~excl; x=A[sel]
    px,py=pos(reg); val=A[int(round(py)),int(round(px))]
    sig=1.4826*np.median(abs(x-np.median(x)))
    nind=sel.sum()*ps**2/(np.pi*(apr*ps)**2)
    o=dict(line=lab,region=reg,ous=ous,value=float(val),field_median=float(np.median(x)),field_sigma_MAD=float(sig),field_std=float(np.std(x)),
           snr_vs_field=float((val-np.median(x))/sig),frac_pixels_above=float((x>=val).mean()),n_pix=int(sel.sum()),n_indep=float(nind),max_field=float(x.max()),p99=float(np.percentile(x,99)))
    out.append(o); print({k:(round(v_,4) if isinstance(v_,float) else v_) for k,v_ in o.items()},flush=True)
    del cube
json.dump(out,open('/workspace/out/s12b_spatial_scatter.json','w'),indent=1)
