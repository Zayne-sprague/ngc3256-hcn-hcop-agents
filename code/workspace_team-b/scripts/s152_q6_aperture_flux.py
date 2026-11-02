"""Q6: aperture fluxes (Jy km/s, Jy) at N and S in Xb4b/Xb4d/Xb4f: compact fraction and recovery. python3 -I s152_q6_aperture_flux.py"""
import sys,json,numpy as np,warnings
sys.path.insert(0,'/workspace/scripts'); warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from regions import REGIONS
C=299792.458;V0=2780.;HW=350.
OUS={'Xb4b':('spw25','spw27','spw25_27_29_31'),'Xb4d':('spw25','spw27','spw25_27_29_31'),'Xb4f':('spw16','spw18','spw16_18_20_22')}
F0={'HCO+3-2':267.5576,'HCN3-2':265.8864}
out={}
for ous,(s_hco,s_hcn,sc) in OUS.items():
    for kind,fn,f0 in (('HCO+3-2',f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.{s_hco}.cube.I.pbcor.fits',F0['HCO+3-2']),('HCN3-2',f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.{s_hcn}.cube.I.pbcor.fits',F0['HCN3-2']),('cont',f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.{sc}.cont.I.pbcor.fits',None)):
        h=fits.open(fn,memmap=True); hd=h[0].header; d=h[0].data
        while d.ndim>(3 if f0 else 2): d=d[0]
        pix=abs(hd['CDELT2'])*3600; bm=(hd['BMAJ']*3600,hd['BMIN']*3600); npb=1.1331*bm[0]*bm[1]/pix**2
        w=WCS(hd).celestial
        if f0:
            nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']; v=C*(1-nu/(f0*1e9)); dv=abs(np.median(np.diff(v)))
            win=np.abs(v-V0)<HW; lf=np.abs(v-V0)>600; i0,i1=np.where(win)[0][[0,-1]]; i1+=1
        for nme in ('N','S'):
            ra,de,_=REGIONS[nme]; cx,cy=[float(t) for t in w.world_to_pixel_values(ra,de)]
            R=int(5/pix)+1; ix,iy=int(round(cx)),int(round(cy)); ys,xs=np.indices((2*R+1,2*R+1))
            if f0:
                cut=np.array(d[:,iy-R:iy+R+1,ix-R:ix+R+1],dtype=float); m=(np.nansum(cut[i0:i1],axis=0)-win.sum()*np.nanmedian(cut[lf],axis=0))*dv
            else: m=np.array(d[iy-R:iy+R+1,ix-R:ix+R+1],dtype=float)
            # centre on local peak within 0.4"
            r0=np.hypot(xs-(cx-(ix-R)),ys-(cy-(iy-R)))*pix; k=np.argmax(np.where(r0<0.4,m,-np.inf)); py,px=np.unravel_index(k,m.shape)
            r=np.hypot(xs-px,ys-py)*pix
            for rad in (0.5,1.0,2.0,4.0):
                if rad<0.6*bm[1]: continue    # aperture smaller than the beam is not a flux
                out.setdefault(kind,{}).setdefault(nme,{}).setdefault(ous,{})[str(rad)]=float(np.nansum(m[r<rad])/npb)
            out[kind][nme][ous]['peak']=float(m[py,px]); out[kind][nme][ous]['beam']=bm
json.dump(out,open('/workspace/out/q6/s152.json','w'),indent=1)
for k,v in out.items():
    for n,vv in v.items():
        print(k,n,{o:{a:(round(b,3) if not isinstance(b,list) else None) for a,b in x.items() if a!='beam'} for o,x in vv.items()})
