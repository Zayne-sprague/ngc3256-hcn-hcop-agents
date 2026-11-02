# Objective region definition from Xb4d (0.95") HCO+ 3-2 and HCN 3-2 moment-0 maps (vsys±350 km/s) and blue-wing map
import sys,warnings,json; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy import ndimage as ndi
RA0,DE0=156.96346,-43.90461
VSYS=2780.
def maps(spw,f0,lo,hi):
    h=fits.open(f('Xb4d',spw),memmap=True); hd=h[0].header; v=(1-freqaxis(hd)/1e9/f0)*C; dv=abs(np.median(np.diff(v)))
    cube=np.array(h[0].data[0],float); sig=np.nanstd(cube[:15],axis=0)
    m=(v>=VSYS+lo)&(v<=VSYS+hi); return hd,np.nansum(cube[m],axis=0)*dv,sig*dv*np.sqrt(m.sum())
def peaks(hd,m0,e0,thr,fp_arcsec=1.5,rmax=14.):
    w=WCS(hd).celestial; pix=abs(hd['CDELT2'])*3600
    ny,nx=m0.shape; Y,X=np.mgrid[:ny,:nx]; ra,de=w.pixel_to_world_values(X,Y)
    dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600; r=np.hypot(dx,dy)
    # noise: local-pbcor-aware, e0 map smoothed
    e1=np.where(np.isfinite(e0),e0,np.nanmedian(e0)); snr=np.nan_to_num(m0)/ndi.uniform_filter(e1,9)
    lm=(m0==ndi.maximum_filter(np.nan_to_num(m0),size=int(fp_arcsec/pix)))&(snr>thr)&(r<rmax)
    ys,xs=np.where(lm); out=[(float(snr[y,x]),float(m0[y,x]),float(ra[y,x]),float(de[y,x]),float(dx[y,x]),float(dy[y,x])) for y,x in zip(ys,xs)]
    return sorted(out,reverse=True),snr,dx,dy
if __name__=='__main__':
    res={}
    for spw,nm,f0 in [('spw25','HCO+ 3-2',267.557626),('spw27','HCN 3-2',265.886434)]:
        hd,m0,e0=maps(spw,f0,-350,350); pk,_,_,_=peaks(hd,m0,e0,8)
        print(nm,'core peaks (S/N>8):')
        for p in pk: print(f'   S/N={p[0]:.0f} I={p[1]:.2f} Jy/b km/s  RA={p[2]:.5f} Dec={p[3]:.5f}  dRA={p[4]:+.2f}" dDec={p[5]:+.2f}"')
        res[nm]=pk
        hd,m0,e0=maps(spw,f0,-450,-150); pk,_,_,_=peaks(hd,m0,e0,5)
        print(nm,'blue wing (-450..-150) peaks (S/N>5):')
        for p in pk: print(f'   S/N={p[0]:.0f} I={p[1]:.2f}  dRA={p[4]:+.2f}" dDec={p[5]:+.2f}"')
        res[nm+' blue']=pk
        hd,m0,e0=maps(spw,f0,150,450); pk,_,_,_=peaks(hd,m0,e0,5)
        print(nm,'red wing (150..450) peaks (S/N>5):')
        for p in pk: print(f'   S/N={p[0]:.0f} I={p[1]:.2f}  dRA={p[4]:+.2f}" dDec={p[5]:+.2f}"')
    json.dump(res,open('/workspace/out/xb4d_peaks.json','w'),indent=1)
