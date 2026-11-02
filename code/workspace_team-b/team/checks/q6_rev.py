import numpy as np, astropy.constants as c
from astropy.io import fits
from scipy.ndimage import gaussian_filter
def tb(S,h,nu):
    om=np.pi/(4*np.log(2))*h['BMAJ']*h['BMIN']*(np.pi/180)**2
    return S*1e-26/om*c.c.value**2/(2*c.k_B.value*nu**2)
fb='/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25.cube.I.pbcor.fits'
fd='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw25.cube.I.pbcor.fits'
hb=fits.getheader(fb);hd=fits.getheader(fd)
mb=fits.open(fb,memmap=True)[0].data[0];md=fits.open(fd,memmap=True)[0].data[0]
nu=265.1e9
# S peak pixel location from Xb4b
sub=np.array(mb[141-3:141+4]).mean(0)
i=np.unravel_index(np.nanargmax(sub),sub.shape);print('S pk px',i,'Tb 7ch-mean',tb(sub[i]*1,hb,nu))
iN=None
# smooth to Xb4d resolution (circular approx geometric-mean kernel)
bd=np.sqrt(hd['BMAJ']*hd['BMIN'])*3600;bb=np.sqrt(hb['BMAJ']*hb['BMIN'])*3600
k=np.sqrt(bd**2-bb**2)/2.3548/0.047
print('kernel px',k)
ch=range(135,148)
sm=np.array([gaussian_filter(np.nan_to_num(mb[j]),k)*(bd/bb)**2 for j in ch])  # per-beam normalisation: Jy/beam scales by beam area ratio
# Xb4d peak near same sky pos
import astropy.wcs as w
wb=w.WCS(hb).celestial;wd=w.WCS(hd).celestial
ra,de=wb.pixel_to_world_values(i[1],i[0]);xd,yd=wd.world_to_pixel_values(ra,de);xd,yd=int(round(float(xd))),int(round(float(yd)))
print('Xb4d pix',xd,yd)
for j,s in zip(ch,sm):
    reg=md[j,yd-2:yd+3,xd-2:xd+3]
    print(j,'Xb4b-smoothed peak(Jy/bm_d) at S: %.4f'%(s[i]*1e0/1), ' Xb4d: %.4f'%np.nanmax(reg))
