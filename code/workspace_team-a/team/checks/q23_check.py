import numpy as np,warnings;warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
P='/archive/files/2015.1.00714.S/member.uid___A001_X2d6_Xae.ari_l.NGC_3256_sci.spw1_229288MHz.12m.cube.I.pbcor.fits'
h=fits.open(P,memmap=True);hd=h[0].header;c=h[0].data[0]
nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3'];v=(1-nu/230.538e9)*299792.458
ba,bb=hd['BMAJ']*3600,hd['BMIN']*3600;pix=abs(hd['CDELT2'])*3600
w=WCS(hd).celestial;R0,D0=156.96346,-43.90461
k=np.where((v>2300)&(v<2844))[0]
for nm,(x,y) in {'N':(-0.06,2.65),'S':(-0.15,-2.59)}.items():
  ix,iy=w.world_to_pixel_values(R0+x/3600/np.cos(np.deg2rad(D0)),D0+y/3600);ix=int(round(float(ix)));iy=int(round(float(iy)))
  sub=np.array(c[k][:,iy-100:iy+101,ix-100:ix+101],float)
  for f in (0.26,0.57):
    s=np.sqrt(f**2-ba*bb)/2.3548/pix
    S=np.array([gaussian_filter(a,s) for a in sub])*f**2/(ba*bb)
    T=1.222e3*S*1e3/(228.4**2*f*f)   # RJ
    # position null test: shifts of 0.05,0.1"
    out=[]
    for d in (0,2,3,5):
      out.append(round(float(np.max([T[:,100+dy,100+dx].max() for dx in (-d,d,0) for dy in (-d,d,0)])),1))
    print(nm,f,'RJ peak at core, max over ring r=0,0.06,0.09,0.15\":',out,'spectrum peak v',v[k][T[:,100,100].argmax()])
    sp=T[:,100,100];print('   FWHM-ish: channels >half max',(sp>sp.max()/2).sum()*2.54,'km/s')
