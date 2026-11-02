# Free Gaussian fit of HOC+ 2-1 (178.972 GHz) in the Xb55 spw18 beam-aperture/peak-pixel spectra at N and S
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy.optimize import curve_fit
z=np.load('/workspace/out/regspec/Xb55.spw18.cube.pbcor.npz'); nu=z['nu']
g=lambda v,a,v0,w,c: a*np.exp(-4*np.log(2)*((v-v0)/w)**2)+c
for r in ('N','S'):
    for st in ('ap','pix'):
        y=z[f'{r}_{st}']; v=vrad(nu,178.97214e9)-VSYS; m=abs(v)<700; fl=m&(abs(v)>400); base=np.median(y[fl]); rms=np.std(y[fl])
        p,cv=curve_fit(g,v[m],y[m]-base,p0=[0.003,0,150,0],sigma=np.full(m.sum(),rms),absolute_sigma=True,bounds=([0,-300,30,-0.002],[0.02,300,600,0.002]))
        print(r,st,'A=%.2f+-%.2f mJy v0=%.0f+-%.0f FWHM=%.0f+-%.0f rms=%.2f'%(p[0]*1e3,np.sqrt(cv[0,0])*1e3,p[1],np.sqrt(cv[1,1]),p[2],np.sqrt(cv[2,2]),rms*1e3))
