# Q7: (a) native-beam peak spectra at N and S (Xb4b, one beam); (b) compact-excess model: HCN = k*HCO+ + A*Gauss(s) (+c) at 0.30" beam about N.
import sys; sys.path.insert(0,'/workspace/scripts')
from q7_lib import *
from scipy.optimize import least_squares
H,O=Cube('Xb4b','HCN'),Cube('Xb4b','HCO+')
kb=lambda f,b: 1.222e6/(f**2*b[0]*b[1])   # K per Jy/beam (RJ)
F={'HCN':265.886434*(1-VSYS/C),'HCO+':267.557626*(1-VSYS/C)}
print('## (a) native-beam peak spectra (Xb4b, 0.25"): peak pixel of +-350 moment map within 0.6" of the nucleus')
for reg,(x0,y0) in {'N':(0.0,2.55),'S':(-0.15,-2.55)}.items():
    for c in (H,O):
        px,py=peak_near(c,c.m0(-350,350),x0,y0); ix=np.argmin(abs(c.X[0]-px)); iy=np.argmin(abs(c.Y[:,0]-py)); s=c.cube[:,iy,ix]
        pk=s[abs(c.v)<350].max(); vp=c.v[abs(c.v)<350][np.argmax(s[abs(c.v)<350])]+VSYS
        half=pk/2; m=(s>half)&(abs(c.v)<350); fw=m.sum()*c.dv   # crude FWHM (channels above half max)
        rms=c.ch_rms[iy,ix]; I=s[abs(c.v)<350].sum()*c.dv
        print(f'  {reg} {c.ln:5s} peak pos ({px:+.2f},{py:+.2f}) peak {pk*1e3:.1f} mJy/beam = {pk*kb(F[c.ln]*1,c.beam):.2f} K at v={vp:.0f} km/s FWHM~{fw:.0f} km/s  I(+-350)={I:.3f} Jy/b km/s  rms/chan(4.4 km/s)={rms*1e3:.2f} mJy/beam  S/N peak {pk/rms:.0f}')
    # ratio of peaks (K) and of integrals at native beams (beam areas differ <3%)
tb=(0.30,0.30,0.); mH,mO=smooth(H.m0(-350,350),H,tb),smooth(O.m0(-350,350),O,tb)
print('\n## (b) compact HCN-excess model at 0.30" beam: pixels within RMAX of the N peak; HCN=k*HCO+ + A*G(s conv beam)+c')
x0,y0=peak_near(O,mO,0.0,2.55); RMAX=2.0
m=np.hypot(O.X-x0,O.Y-y0)<RMAX; X,Y=O.X[m]-x0,O.Y[m]-y0; h=mH[m]; o=mO[m]
def G(s): 
    w=np.sqrt(0.30**2+s**2); g=np.exp(-4*np.log(2)*(X**2+Y**2)/w**2); return g/ (1.0)   # peak-normalised (Jy/beam of the beam-convolved source)
sig_pix=np.std(smooth(np.random.default_rng(1).standard_normal(O.m0(-350,350).shape)*0,O,tb))  # placeholder
res=[]
for s in [0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.6,1.0]:
    A_=np.vstack([o,G(s),np.ones_like(o)]).T; c,*_=np.linalg.lstsq(A_,h,rcond=None); r=h-A_@c; chi=np.sum(r**2)
    # integrated excess flux (Jy km/s): peak amplitude c[1] (Jy/b km/s) * (2D integral of normalised gaussian / beam area)
    w=np.sqrt(0.30**2+s**2); flux=c[1]*(np.pi*w**2/(4*np.log(2)))/(OM(tb)) 
    res.append((s,chi,c[0],c[1],c[2],flux))
c0=min(r[1] for r in res)
print('  s_FWHM[arcsec]  rss/rss_min   k     A[Jy/b km/s]  c   HCN_excess_flux[Jy km/s]')
for r in res: print(f'  {r[0]:.2f}  {r[1]/c0:6.3f}  {r[2]:.3f}  {r[3]:.3f}  {r[4]:+.3f}  {r[5]:.3f}')
# same with k fixed from 1-2" (no free k): and a pure linear k only
A_=np.vstack([o,np.ones_like(o)]).T; c,*_=np.linalg.lstsq(A_,h,rcond=None); print(f'  k-only (no core) rss/rss_min={np.sum((h-A_@c)**2)/c0:.2f}  k={c[0]:.3f}')
# annulus profile of excess: HCN - k*HCO+ with k=0.38 (S-like ratio) at 0.30"
for k in (0.38,0.42):
    print(f'  excess HCN-k*HCO+ (k={k}) cumulative, Jy km/s:',[round(annulus_flux(H,mH-k*mO,x0,y0,0,b,tb)[0],3) for b in (0.15,0.25,0.5,1,1.5,2,3)])
