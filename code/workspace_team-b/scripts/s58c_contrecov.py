# Q12 round 2, step c: continuum as an independent probe of what each array keeps (continuum is smooth and line-free). 10" circular beam, circular aperture sums about the
# N-S midpoint (Jy, = sum/beam area), PB-corrected maps, apertures r<4,8,12".  Ratios Xb4b/Xb4d/Xb4f (same Band 6 band) and the ACA spectral indices Xb55->Xb4f->Xb47.
# usage: python3 -I s58c_contrecov.py -> out/q12b/cont_recovery.txt/.json
import sys,json; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from common import *
import s56_lib as L
from s30_matchedJ_ratios import transfer
CEN=L.CEN; B10=(10.,10.,0.); log=[]
def P(s=''): print(s,flush=True); log.append(s)
NU={'Xb55':'175.5','Xb4f':'265.7','Xb47':'351.5'}
def conv2(img,ps,b0):
    n=img.shape[0]; Pn=int(2**np.ceil(np.log2(2*n))); T=transfer((Pn,Pn),ps,b0,B10)[:,:Pn//2+1]; fac=B10[0]*B10[1]/(b0[0]*b0[1])
    a=np.zeros((Pn,Pn),np.float32); a[:n,:n]=np.nan_to_num(img); return np.fft.irfft2(np.fft.rfft2(a)*T,s=(Pn,Pn))[:n,:n]*fac
RES={}
for ous in ('Xb55','Xb4f','Xb47','Xb4d','Xb4b'):
    d=np.load(f'/workspace/out/q12b/cont_{ous}.npz'); ps=float(d['ps']); b0=tuple(d['beam']); img=d['img']
    fn=files(f'*{ous}.ngc_3256_sci.spw*_*.cont.I.pbcor.fits')[0]; hd=fits.getheader(fn); wc=WCS(hd).celestial
    k=1
    if ous=='Xb4b': k=4; n2=img.shape[0]//4*4; img=img[:n2,:n2].reshape(n2//4,4,n2//4,4).mean(axis=(1,3)); ps*=4
    if ous=='Xb4d': pass
    if ous=='Xb4b':
        # crop +-19" around the midpoint to bound the FFT size
        cx,cy=wc.world_to_pixel_values(*CEN); cx/=k; cy/=k; h=int(19/ps); x0=int(cx-h); y0=int(cy-h); img=img[y0:y0+2*h,x0:x0+2*h]; cx-=x0; cy-=y0
    else: cx,cy=wc.world_to_pixel_values(*CEN)
    m=conv2(img,ps,b0); n=m.shape[0]; Y,X=np.mgrid[:n,:n]; rr=np.hypot(X-cx,Y-cy)*ps; area=1.1331*100/ps**2
    RES[ous]={R:float(m[rr<=R].sum()/area) for R in (4,8,12)}; RES[ous]['peak']=float(m[rr<4].max())
    P(f'{ous} cont native beam {b0[0]:.2f}x{b0[1]:.2f}": 10" aperture flux r<4/8/12" = '+' / '.join(f'{RES[ous][R]*1e3:.1f}' for R in (4,8,12))+f' mJy ; peak {RES[ous]["peak"]*1e3:.1f} mJy/10"-beam')
P('\n## Band-6 12-m / ACA continuum at 10" (same band: pure array effect), and the line comparison (HCO+ 3-2 r<4/8/12": Xb4d/Xb4f 0.81/0.65/0.50; Xb4b/Xb4f 0.69/0.47/-0.01)')
for a,b in (('Xb4d','Xb4f'),('Xb4b','Xb4f'),('Xb4b','Xb4d')):
    P(f'cont {a}/{b} r<4/8/12": '+' / '.join(f'{RES[a][R]/RES[b][R]:.2f}' for R in (4,8,12)))
P('\n## ACA spectral indices of the continuum in the 10" aperture (S ~ nu^alpha); a loss growing with frequency would make alpha(4f->47) < alpha(55->4f) (curvature negative)')
nu={'Xb55':170.74,'Xb4f':257.71,'Xb47':346.53}   # CRVAL3 of the cont products (mean of the 4 spw), GHz
for R in (4,8,12):
    a1=np.log(RES['Xb4f'][R]/RES['Xb55'][R])/np.log(nu['Xb4f']/nu['Xb55']); a2=np.log(RES['Xb47'][R]/RES['Xb4f'][R])/np.log(nu['Xb47']/nu['Xb4f'])
    e=0.1*np.sqrt(2)/np.log(257.71/170.74)   # 10% calibration per band
    P(f'r<{R}": alpha(55->4f)={a1:.2f}  alpha(4f->47)={a2:.2f}  curvature={a2-a1:+.2f} (+-{0.1*np.sqrt(2)/np.log(346.53/257.71):.2f}..{e:.2f} from 10% calibration/band)')
json.dump(RES,open('/workspace/out/q12b/cont_recovery.json','w')); open('/workspace/out/q12b/cont_recovery.txt','w').write('\n'.join(log))
