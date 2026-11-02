import numpy as np, glob, os
from astropy.io import fits
from astropy.wcs import WCS
C=299792.458
D='/data/'
def f(ous,spw,kind='cube.I.pbcor',sep='.'):
    return D+f'member.uid___A001_X1273_{ous}{sep}ngc_3256_sci.{spw}.{kind}.fits'
def load_cube(path):
    h=fits.open(path,memmap=True)
    hd=h[0].header
    d=h[0].data
    return hd,d
def freqaxis(hd):
    n=hd['NAXIS3']; return hd['CRVAL3']+(np.arange(n)+1-hd['CRPIX3'])*hd['CDELT3']
