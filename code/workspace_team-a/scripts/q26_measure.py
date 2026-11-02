# Q26: annulus HCN/HCO+ from the common-uv re-imaged cubes (CASA, scripts/q26_image.py): u>6/8 kl, 8" taper, 10" circular restoring beam.
# PB-corrected within PB>0.5, annuli (0-4,4-8,8-12") about the N-S midpoint (-0.075,0) offset (as Q11), window VSYS+-350 km/s.
import numpy as np, json, glob
from astropy.io import fits
from astropy.wcs import WCS
D='/workspace/out/q26/casa/'; RA0,DE0=156.96346,-43.90461; C=299792.458; VS=2780.; HW=350.
ANN=[(0,4),(4,8),(8,12)]; XM,YM=-0.075,0.
REST={('HCN','2-1'):177.2611,('HCO+','2-1'):178.3751,('HCN','3-2'):265.8864,('HCO+','3-2'):267.5576}
def meas(ous,sp,J,u,hw=HW,vs=VS,pbl=0.5):
    nm=f'{ous}_{sp.replace("+","p")}_{J}_u{u}'
    h=fits.open(D+nm+'.image.pbcor.fits'); hd=h[0].header; im=np.array(h[0].data,float)
    pb=np.array(fits.getdata(D+nm+'.pb.fits'),float); ps=fits.getdata(D+nm+'.psf.fits')
    if pb.ndim==3: pb=pb[0]
    fr=hd['CRVAL3']+hd['CDELT3']*np.arange(im.shape[0]); v=(1-fr/1e9/REST[(sp,J)])*C; dv=abs(np.median(np.diff(v)))
    w=WCS(hd).celestial; Y,X=np.mgrid[:im.shape[1],:im.shape[2]]; ra,de=w.pixel_to_world_values(X,Y)
    dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600; r=np.hypot(dx-XM,dy-YM)
    pix=abs(hd['CDELT2'])*3600; omp=np.pi*100/(4*np.log(2))/pix**2
    good=np.isfinite(im).all(0)&(pb>pbl)
    rows=[]
    spec={}
    for a,b in ANN:
        m=(r>=a)&(r<b)&good; frac=m.sum()/((r>=a)&(r<b)).sum()
        s=np.nansum(np.where(m,im,0),axis=(1,2))/omp     # Jy per channel (beam-integrated), annulus sum
        spec[(a,b)]=s
        sel=(abs(v-vs)<=hw); nul=(abs(v-vs)>600)&(np.abs(s)>0)
        # nulls: only channels with data (PSF not blank) and beyond 600 km/s; noise of an N-channel window = std(ch)*sqrt(N)
        f=s[sel].sum()*dv; N=sel.sum(); sg=np.std(s[nul])*np.sqrt(N)*dv if nul.sum()>5 else np.nan
        rows.append((f,sg,frac,nul.sum()))
    return rows,v,spec
if __name__=='__main__':
    out={}
    SET=[('Xb55','2-1'),('Xb4f','3-2'),('Xb4d','3-2')]
    for u in (6,8):
        for ous,J in SET:
            H,vh,_=meas(ous,'HCN',J,u); O,vo,_=meas(ous,'HCO+',J,u)
            for i,(a,b) in enumerate(ANN):
                (fh,eh,fr,_),(fo,eo,_,_)=H[i],O[i]; R=fh/fo; eR=abs(R)*np.hypot(eh/fh,eo/fo)
                print(f'u>{u}kl {ous} {J} {a}-{b}": HCN {fh:7.2f}+-{eh:.2f} HCO+ {fo:7.2f}+-{eo:.2f} Jy km/s  R={R:.3f}+-{eR:.3f}  PB>0.5 area frac {fr:.2f}')
                out[f'u{u}|{ous}|{J}|{a}-{b}']=dict(HCN=fh,eHCN=eh,HCO=fo,eHCO=eo,R=R,eR=eR,frac=fr)
    json.dump(out,open('/workspace/out/q26/q26_ratios.json','w'),indent=1,default=float)
