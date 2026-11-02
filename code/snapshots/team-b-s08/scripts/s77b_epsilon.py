# Q2 round 2, fix 2(ii): size of the dirty-vs-clean beam scale factor epsilon for UNCLEANED (residual) emission, by synthetic uv sampling.
# The wing is outside the clean mask (s75c), so its Jy/beam is Jy per DIRTY beam, but flux = sum/(clean beam area). We simulate: uv density D(u,v) (elliptical, hole at short spacings)
# tuned so that the main lobe of the dirty beam has the Xb4d beam (0.93x0.74", PA -86; also Xb4b 0.25x0.23 and the 1" circular smoothing of both).  A circular-Gaussian source (FWHM 0.5-3") is
# observed: dirty image = FT(D*V). Measured flux in an r=1.5" aperture = sum(dirty)/beam_area(clean, fitted to the main lobe) vs the same aperture on (source (x) clean beam) (what a perfect clean gives, in the same units)
# and vs the true source flux in the same aperture.   epsilon_eff = F_dirty_ap / F_clean_ap.  uv models: (G) Gaussian density; (P) 1/r radial density (typical of a filled array) with outer taper.
# Output out/q2r2/s77b.txt/json. Assumption: no real uv coverage is available, so this is an order-of-magnitude bound on the direction/size of the effect (state as such).
import json, numpy as np
from scipy.optimize import curve_fit
OUT='/workspace/out/q2r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
PS=0.05; N=1024; a=PS/206265.0
fx=np.fft.fftfreq(N,a)   # cycles/rad  (u in wavelengths)
U,V=np.meshgrid(fx,fx)
def dens(model,sig_u,sig_v,pa_deg,umin,shape=1.0):
    th=np.deg2rad(pa_deg); # elongation along PA: beam major axis at PA -> uv minor axis. rotate uv by PA
    ur=U*np.cos(th)-V*np.sin(th); vr=U*np.sin(th)+V*np.cos(th)   # ur along RA-ish rotated
    r=np.hypot(ur/sig_u,vr/sig_v)
    if model=='G': D=np.exp(-0.5*r**2)
    else: D=1.0/np.maximum(r,0.05)*np.exp(-0.5*(r/2.2)**2)         # ~1/r radial density, tapered
    D=np.where(np.hypot(U,V)<umin,0,D); return D
def beam_fit(D):
    B=np.fft.fftshift(np.fft.ifft2(D).real); B/=B.max(); c=N//2; h=int(4.0/PS); sub=B[c-h:c+h+1,c-h:c+h+1]
    Y,X=np.mgrid[-h:h+1,-h:h+1]*PS
    def g(p,A,sx,sy,t): x,y=p; u=x*np.cos(t)+y*np.sin(t); w=-x*np.sin(t)+y*np.cos(t); return (A*np.exp(-0.5*((u/sx)**2+(w/sy)**2))).ravel()
    m=np.hypot(X,Y)<0.9
    # fit main lobe only: use points inside the first null region (|x|<1.2*fwhm); approximate by r<0.9"
    p,_=curve_fit(lambda xy,A,sx,sy,t: g((xy[0],xy[1]),A,sx,sy,t)[m.ravel()] if False else g((xy[0],xy[1]),A,sx,sy,t),(X.ravel(),Y.ravel()),sub.ravel(),p0=[1,0.3,0.3,0.0],maxfev=20000)
    return B,abs(p[1])*2.355,abs(p[2])*2.355
def tune(model,fwhm_maj,fwhm_min,pa,umin):
    sig=lambda f: 2.355/(2*np.pi*(f/206265.0))*1.0
    s_u,s_v=sig(fwhm_min),sig(fwhm_maj); best=None
    for k in np.linspace(0.5,2.5,17):
        D=dens(model,s_u*k,s_v*k,pa,umin); B,f1,f2=beam_fit(D); fa,fb=max(f1,f2),min(f1,f2)
        e=(fa/fwhm_maj-1)**2+(fb/fwhm_min-1)**2
        if best is None or e<best[0]: best=(e,k,D,B,fa,fb)
    return best
def run(label,model,maj,minr,pa,umin):
    e,k,D,B,fa,fb=tune(model,maj,minr,pa,umin); Y,X=np.mgrid[-N//2:N//2,-N//2:N//2]*PS
    cb=np.exp(-0.5*(((X*np.cos(np.deg2rad(pa))+Y*np.sin(np.deg2rad(pa))))**2/(fa/2.355)**2*0+0))  # placeholder not used
    # clean beam = Gaussian with fitted FWHM (axis angle ignored: circular-equivalent area)
    barea=1.1331*fa*fb/PS**2
    sl=np.sort(B.ravel()); neg=B[B<0].sum()/B.sum(); pos=B[B>0].sum()
    P(f'## {label}: tuned main lobe {fa:.2f}x{fb:.2f}" (target {maj}x{minr}); sum(dirty beam)/pixels in main-lobe-area units = {B.sum()/barea:.3f}; peak sidelobe {B[np.hypot(X,Y)>1.2*fa].max():.3f}, min {B.min():.3f}')
    out={}
    for fs in (0.5,1.0,1.5,2.0,3.0):
        # source FWHM fs (sky), circular gaussian at origin; visibility exp(-2 pi^2 s^2 q^2)
        q2=U**2+V**2; Vs=np.exp(-2*np.pi**2*(fs/2.355/206265.0)**2*q2)
        dirty=np.fft.fftshift(np.fft.ifft2(D*Vs).real)*1.0       # in units of source peak (D(0)=? normalise below)
        d0=D.sum()/N**2; dirty/=d0                                 # normalised so that dirty beam peak=1 -> Jy/dirty beam for unit-peak point source flux
        # ideal clean image: source (x) clean Gaussian beam (unit-peak), in Jy/clean beam
        sc=np.sqrt((fs/2.355)**2+(np.sqrt(fa*fb)/2.355)**2); Cimg=np.exp(-0.5*(X**2+Y**2)/sc**2)*(np.sqrt(fa*fb)/2.355)**2/sc**2*(fs/2.355*0+1)   # peak of source*beam: (b^2/(s^2+b^2)) * source peak(in Jy/beam units where point source flux=1)
        # in units: image in Jy/beam of unit-flux source: dirty beam peak=1 for point source; extended source of flux F: peak = F*b^2/(b^2+s^2)  (F=1 here)
        Cimg=np.exp(-0.5*(X**2+Y**2)/sc**2)*((np.sqrt(fa*fb)/2.355)**2/sc**2)
        # source total flux = 1 -> the Gaussian: peak surface brightness per clean beam = (b/ s_conv)^2
        for rad in (1.0,1.5,2.5):
            m=np.hypot(X,Y)<=rad
            Fd=dirty[m].sum()/barea; Fc=Cimg[m].sum()/barea
            Ft=(1-np.exp(-0.5*rad**2/(fs/2.355)**2))      # true flux in aperture for unit total flux
            out[f'{fs}|{rad}']=(Fd,Fc,Ft)
        r=out[f'{fs}|1.5']; P(f'   source FWHM {fs}": r=1.5" aperture flux (unit total flux): dirty-image {r[0]:.3f}, ideal clean image {r[1]:.3f}, true-in-aperture {r[2]:.3f}; epsilon_eff=dirty/clean {r[0]/r[1]:.3f}; dirty/true {r[0]/r[2]:.3f}; (r=1.0": {out[f"{fs}|1.0"][0]/out[f"{fs}|1.0"][1]:.3f}; r=2.5": {out[f"{fs}|2.5"][0]/out[f"{fs}|2.5"][1]:.3f})')
    return out
res={}
for model in ('G','P'):
    for umin in (13e3,26e3):
        res[f'Xb4d|{model}|{umin:.0f}']=run(f'Xb4d-like {model} hole {umin/1e3:.0f} klambda',model,0.93,0.74,-86,umin)
    res[f'Xb4b|{model}']=run(f'Xb4b-like {model}',model,0.25,0.23,-49,40e3)
json.dump({k:{kk:list(vv) for kk,vv in v.items()} for k,v in res.items()},open(OUT+'/s77b.json','w'),indent=1)
open(OUT+'/s77b.txt','w').write('\n'.join(txt))
