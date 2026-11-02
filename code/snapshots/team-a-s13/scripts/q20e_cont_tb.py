# Q20e: 1.2 mm continuum brightness temperature at 0.24" (Xb4b cont) at N and S peaks: Td lower bound; compare with CO(2-1) Tpk (Brunetti+2021 Table 2).
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
f='/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'
h=fits.open(f)[0]; im=np.squeeze(h.data); w=WCS(h.header).celestial
bmaj,bmin=h.header['BMAJ']*3600,h.header['BMIN']*3600; nu=h.header.get('CRVAL3',257.76e9)
nu=257.76e9
om=np.pi*bmaj*bmin/(4*np.log(2))*(np.pi/180/3600)**2
print('beam',bmaj,bmin)
def TB(S): # Jy/beam -> RJ T_b ; also Planck-J
    return 1.222e6*S/(nu/1e9)**2/(bmaj*bmin)
import pandas as pd
reg=pd.read_csv('/workspace/out/regions.csv'); print(reg.to_string()[:800])
for name in ('N','S'):
    r=reg[reg.iloc[:,0].astype(str)==name].iloc[0]
    print(name,r.to_dict())
print()
pb=None
for name in ('N','S'):
    r=reg[reg.region==name].iloc[0]
    x,y=w.all_world2pix(r.RA_deg,r.Dec_deg,0); x,y=int(round(float(x))),int(round(float(y)))
    cut=im[y-6:y+7,x-6:x+7]; S=np.nanmax(cut)
    print(name,'peak cont Jy/beam %.5f  RJ Tb %.2f K'%(S,TB(S)))
    hv=0.04799243*257.76  # h nu/k  K
    Tpl=hv/np.log(1+hv/ (TB(S)))  # Planck-equivalent (blackbody, tau>>1) : J(T)=Tb -> T
    print('   thick-blackbody Td lower bound (J(Td)-J(Tbg)=Tb): %.1f K'%(hv/np.log(1+hv/(TB(S)+hv/np.expm1(hv/2.73)))))
