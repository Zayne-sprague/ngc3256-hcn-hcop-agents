# Q11 like-for-like uv-matched test: SAME Fourier band-pass (cosine taper u1->u2) applied to the 10" maps of BOTH 2-1 (Xb55) and 3-2 (Xb4f, Xb4d), annuli where both survive positive;
# fraction of flux kept per annulus/line; also edge-aware normalised convolution restricted to a PB>0.5 mask (Gaussian PB, FWHM 22" 12m-like / 37" 7m-like at 266 GHz, 55" at 177 GHz).
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
tb=(10.,10.,0.); ANN=[(0,4),(4,8),(6,9),(8,12)]
def sm(d,pbfw=None,pbmin=0.5):
    k=kernel(d['hd'],tb,d['beam']); a=np.nan_to_num(d['m0']); r2=d['dx']**2+d['dy']**2
    if pbfw is None: return fftconvolve(a,k,mode='same')*OM(tb)/OM(d['beam'])
    m=(np.exp(-4*np.log(2)*r2/pbfw**2)>=pbmin).astype(float)           # normalised convolution: divide by convolved mask
    num=fftconvolve(a*m,k,mode='same'); den=fftconvolve(m,k,mode='same')
    return np.where(den>0.2*k.sum(),num/np.maximum(den,1e-9)*k.sum(),np.nan)*OM(tb)/OM(d['beam'])
def bp(img,pix,u1,u2,npad=512):
    n=img.shape; P=np.zeros((npad,npad)); P[:n[0],:n[1]]=img
    F=np.fft.fft2(P); u=np.fft.fftfreq(npad,d=pix/206265.)/1e3; U=np.hypot(*np.meshgrid(u,u))
    t=np.clip((U-u1)/max(u2-u1,1e-9),0,1); F*=0.5*(1-np.cos(np.pi*t))
    return np.real(np.fft.ifft2(F))[:n[0],:n[1]]
def ann(d,s):
    r=np.hypot(d['dx']+0.075,d['dy']); omp=OM(tb)/d['pix']**2
    return np.array([np.nansum(s[(r>=a)&(r<b)])/omp for a,b in ANN])
M={'2-1 Xb55':('HCN 2-1','HCO+ 2-1'),'3-2 Xb4f':('HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'),'3-2 Xb4d':('HCN 3-2@Xb4d','HCO+ 3-2@Xb4d')}
D={m:[load(k,nullspec=False) for k in v] for m,v in M.items()}
print('annuli',ANN)
print('== no filter')
S={m:[sm(d) for d in v] for m,v in D.items()}; F0={}
for m in M:
    h,o=[ann(d,s) for d,s in zip(D[m],S[m])]; F0[m]=(h,o); print(f'{m:9s} HCN {h.round(2)} HCO+ {o.round(2)} R {(h/o).round(3)}')
for (u1,u2) in ((0,0),(6,10),(8,12),(10,14)):
    print(f'== common band-pass cosine taper u1={u1} u2={u2} kl (applied to every map)')
    for m in M:
        fh=[ann(d,bp(s,d['pix'],u1,u2)) for d,s in zip(D[m],S[m])]; h,o=fh
        ok=(h>0)&(o>0); R=np.where(ok,h/np.where(o==0,1,o),np.nan)
        print(f'{m:9s} HCN {h.round(2)} HCO+ {o.round(2)} R {R.round(3)}  kept HCN {(h/F0[m][0]).round(2)} HCO+ {(o/F0[m][1]).round(2)}')
print('== PB>0.5 masked, edge-aware normalised convolution (unfiltered)')
for m,pb in (('2-1 Xb55',55.),('3-2 Xb4f',37.),('3-2 Xb4f',22.),('3-2 Xb4d',22.)):
    h,o=[ann(d,sm(d,pb)) for d in D[m]]; print(f'{m:9s} PBfwhm {pb:.0f}" HCN {h.round(2)} HCO+ {o.round(2)} R {(h/o).round(3)}')
