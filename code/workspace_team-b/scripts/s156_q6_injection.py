"""Q6: injection test of the half-max deconvolved size (s151.hm_stats) with a Gaussian disk on a noisy map (beam 0.245x0.235", pix 0.047"). python3 -I s156_q6_injection.py"""
import sys,numpy as np; sys.path.insert(0,'/workspace/scripts')
import s151_q6_size_width_tb as S
from scipy.ndimage import gaussian_filter
rng=np.random.default_rng(5); pix=0.047; bm=(0.245,0.235); n=107; c=n//2
ys,xs=np.indices((n,n)); r2=((xs-c)**2+(ys-c)**2)*pix**2
for snr in (50,20,8):
    for fwhm in (0.0,0.15,0.3,0.5,0.8):
        # true sky Gaussian FWHM fwhm convolved with beam geometric mean -> observed FWHM sqrt(fwhm^2+bg^2)
        bg=np.sqrt(bm[0]*bm[1]); obs=np.sqrt(fwhm**2+bg**2); sky=np.exp(-4*np.log(2)*r2/obs**2)
        sig=1/snr; got=[]
        for k in range(20):
            m=sky+S.make_noise(sky.shape,bm,pix,sig,rng)
            a,pk,_=S.halfmax(m,pix,c,c,1.5); got.append(np.sqrt(max(4*a/np.pi-bg**2,0)))
        print('SNR %3d true %.2f -> recovered %.3f +- %.3f'%(snr,fwhm,np.mean(got),np.std(got)))
