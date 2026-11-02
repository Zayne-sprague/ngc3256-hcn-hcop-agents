"""Q62: 20" aperture (midpoint) line fluxes [Jy km/s] of the CO cubes vs Sakamoto+2014 (CO 1-0/2-1/3-2 = 1.0e3 / 3.0e3 (4.0e3 at 3" res) / 5.7e3 Jy km/s, no missing-flux correction).
Native resolution pbcor cubes, |v-V0|<400 km/s, linear baseline from 700<|v-V0|<1500.  usage: python3 -I s373_q62_fluxcheck.py"""
import sys,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s180_q29_load as L
from astropy.io import fits
C=L.C; V0=2780.; P='/archive/files/'
CUBES={'CO10_X77':(P+'2016.2.00094.S/member.uid___A001_X1262_X77.NGC3256_sci.spw16.cube.I.pbcor.fits',115.2712),
 'CO10_Xd5':(P+'2016.2.00042.S/member.uid___A001_X1262_Xd5.ngc3256_sci.spw16.cube.I.pbcor.fits',115.2712),
 'CO10_X4f':(P+'2018.1.00223.S/member.uid___A001_X133f_X4f.NGC3256_sci.spw16.cube.I.pbcor.fits',115.2712),
 'CO21_7m':(P+'2015.1.00714.S/member.uid___A001_X2d6_Xb2.NGC_3256_sci.spw0_227531MHz.mos.7m.cube.I.pbcor.fits',230.538),
 'CO32_X67':(P+'2016.2.00094.S/member.uid___A001_X1262_X67.NGC3256_sci.spw20.cube.I.pbcor.fits',345.796),
 'CO32_X85a':(P+'2018.1.00223.S/member.uid___A001_X133d_X85a.NGC3256_sci.spw20.cube.I.pbcor.fits',345.796)}
CEN=(156.96343,-43.90461)
for tag,(path,f0) in CUBES.items():
    cut,nu,pix,bm,pxy=L.load(path,half=40.)
    v=C*(1-nu/f0); o=np.argsort(v); v=v[o]; cut=np.nan_to_num(cut[o])
    x0,y0=pxy(CEN); ny,nx=cut.shape[1:]; yy,xx=np.mgrid[:ny,:nx]; ap=np.hypot(xx-x0,yy-y0)*pix<=10.
    spec=cut[:,ap].sum(1)/(1.1331*bm[0]*bm[1]/pix**2)   # Jy per beam-area sum -> Jy
    dv=abs(np.gradient(v)); fr=(abs(v-V0)>700)&(abs(v-V0)<1500); w=abs(v-V0)<400
    p=np.polyfit(v[fr],spec[fr],1); s=spec-np.polyval(p,v)
    print(f'{tag:10s} beam {bm[0]:.1f}x{bm[1]:.1f}  S(20" ap)={ (s[w]*dv[w]).sum():8.0f} Jy km/s  (nfree {fr.sum()})')
