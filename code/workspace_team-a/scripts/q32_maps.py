# Q32 shared: load 0.18" CO(2-1) cube, build native-resolution maps (K RJ) around S and N for several velocity measures
import numpy as np, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
C=299792.458; F0=230.538e9
P='/archive/files/2015.1.00714.S/member.uid___A001_X2d6_Xae.ari_l.NGC_3256_sci.spw1_229288MHz.12m.cube.I.pbcor.fits'
RA0,DE0=156.96346,-43.90461
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
WIN={'S':(2711,100),'N':(2806,155)}
def load(HWas=6.0):
    h=fits.open(P,memmap=True); hd=h[0].header; cube=h[0].data[0]
    nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/F0)*C
    bmaj,bmin=hd['BMAJ']*3600,hd['BMIN']*3600; pix=abs(hd['CDELT2'])*3600
    w=WCS(hd).celestial; HW=int(round(HWas/pix)); K=1.222e3/((nu/1e9)**2*bmaj*bmin)*1e3
    out={}
    for nm,(x0,y0) in CORES.items():
        ra=RA0+x0/3600/np.cos(np.deg2rad(DE0)); de=DE0+y0/3600
        ix,iy=[int(round(float(a))) for a in w.world_to_pixel_values(ra,de)]
        sub=np.array(cube[:,iy-HW:iy+HW+1,ix-HW:ix+HW+1],float)*K[:,None,None]
        v0,fw=WIN[nm]; ks=np.where((v>v0-fw/2)&(v<v0+fw/2)&(v<2844))[0]
        kp=np.where((v>2300)&(v<2844))[0]; kpk=kp[np.argmax(sub[kp][:,HW,HW])]
        off=np.where(v<2300)[0]
        maps={'peak':sub[kpk],'3ch':sub[kpk-1:kpk+2].mean(0),'HCNwin':sub[ks].mean(0)}
        noise={k:float(np.std(sub[off][:,:,:].reshape(len(off),-1).mean(1)*0+sub[off][:,HW,HW])) for k in ['x']}
        rms1=float(np.std(sub[off][:,HW,HW])); 
        out[nm]=dict(maps=maps,HW=HW,rms1=rms1,vpk=float(v[kpk]),nch_win=len(ks))
    return out,(bmaj,bmin,pix)
