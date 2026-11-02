# Q11b: uv-matched test. Xb55 (177 GHz) 2-1 10" maps are high-passed in the Fourier domain at the minimum uv radius of the Band-6 3-2 data
# (Bmin/lambda at 265/267 GHz; 7-m ACA Bmin~8.7 m, 12-m Xb4d Bmin~15 m), so that the 2-1 map lacks the same angular scales as 3-2. Compare 2-1(filtered) vs 3-2 annuli,
# and window / baseline-order sensitivity of the 3-2 annulus ratios.
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
tb=(10.,10.,0.); ANN=[(0,4),(4,8),(8,12)]
def sm(d,arr=None):
    k=kernel(d['hd'],tb,d['beam']); a=np.nan_to_num(d['m0'] if arr is None else arr); return fftconvolve(a,k,mode='same')*OM(tb)/OM(d['beam'])
def hp(img,pix,umin_kl,npad=512):
    n=img.shape; P=np.zeros((npad,npad)); P[:n[0],:n[1]]=img
    F=np.fft.fft2(P); u=np.fft.fftfreq(npad,d=pix/206265.)/1e3
    U=np.hypot(*np.meshgrid(u,u)); F[U<umin_kl]=0
    return np.real(np.fft.ifft2(F))[:n[0],:n[1]]
def ann(d,s,fl=(30,45)):
    r=np.hypot(d['dx']+0.075,d['dy']); omp=OM(tb)/d['pix']**2
    z=np.nanmean(s[(r>=fl[0])&(r<fl[1])]) if fl else 0
    return [(s[(r>=a)&(r<b)]-z).sum()/omp for a,b in ANN]
c=299792458.
H=load('HCN 2-1',nullspec=False); O=load('HCO+ 2-1',nullspec=False)
sH,sO=sm(H),sm(O)
print('Xb55 2-1 unfiltered      HCN',np.round(ann(H,sH,None),2),'HCO+',np.round(ann(O,sO,None),2),'R',np.round(np.array(ann(H,sH,None))/np.array(ann(O,sO,None)),3))
for B in (8.7,15.0):
  umin=B/(c/265.9e9)/1e3
  for lab,fl in (('no zero',None),('far-field',(30,45))):
    hH,hO=hp(sH,H['pix'],umin),hp(sO,O['pix'],umin)
    aH,aO=ann(H,hH,fl),ann(O,hO,fl)
    print(f'Xb55 2-1 high-passed u>{umin:.1f}kl (B={B} m) [{lab}]: HCN',np.round(aH,2),'HCO+',np.round(aO,2),'R',np.round(np.array(aH)/np.array(aO),3))
# window / baseline sensitivity of the 3-2 annulus ratios (Xb4f), 10"
print('--- 3-2 Xb4f window sensitivity (no zero level / far-field 30-45)')
for hw in (250.,350.,450.):
  r={}
  for k in ('HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'):
    d=load(k,nullspec=False,window=(VSYS-hw,VSYS+hw)); r[k]=(d,sm(d))
  for lab,fl in (('none',None),('ff',(30,45))):
    aH=np.array(ann(*r['HCN 3-2@Xb4f'],fl)); aO=np.array(ann(*r['HCO+ 3-2@Xb4f'],fl))
    print(f' +-{hw:.0f} [{lab}] HCN',np.round(aH,2),'HCO+',np.round(aO,2),'R',np.round(aH/aO,3))
