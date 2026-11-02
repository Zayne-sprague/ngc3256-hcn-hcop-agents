# Q12 round 2, step f: forward-model test of the recovery/power-law/extrapolation recipe on SIMULATED skies passed through assumed uv coverage with NONLINEAR (Hogbom) cleaning.
# Sky: the s58d Xb4f HCO+ 3-2 component fit (N,S,Bridge cores/halos, E1,E2,W Gaussians, true Jy km/s) plus an optional extended Gaussian (FWHM 25/40/80", flux = f_ext x nuclear-complex total)
#   that no ACA/12-m data constrain (the unknown outer-scale structure).  Same sky at every J (kJ=1) unless stated.  Primary beam multiplies the sky (12 m 21.9", ACA 1.13 lambda/7m scaling).
# uv: baselines with radial pdf ~ L^-p (p=1), isotropic, between Lmin and Lmax (m); gridded natural weighting on a 512^2 / 0.4" grid; Hogbom clean (gain 0.1, 3000 it, stop 1.5% of peak);
#   restored with a Gaussian fitted to the dirty-beam core; PB-corrected; then smoothed to a 10" circular beam and compared to the TRUE sky smoothed to 10" at the six positions.
#   Arrays (ASSUMED, in metres): Xb4b-like 12-m: Lmin 60, Lmax 400; Xb4d-like 12-m: Lmin 29 (MRS ratio d/f = 0.30), Lmax 250; ACA (Xb4f/Xb47/Xb55): Lmin 8.7, Lmax 45 at every band (fixed in metres,
#   so MRS ~ lambda, beam ~ lambda). The 12-m C43 shortest baseline is really ~15 m; 29 m is an effective value that reproduces the nominal MRS ratio.
# Outputs: out/q12b/sim.txt/.json: absolute recovery R per array/position/scenario, power-law slopes d/f and b/d, predicted vs simulated R55/R4f and R47/R4f.
# usage: python3 -I s58f_simulate.py
import sys,json; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve
from regions import REGIONS
import s56_lib as L
CEN=L.CEN; COS=np.cos(np.deg2rad(CEN[1])); FW=2*np.sqrt(2*np.log(2)); NG=512; PX=0.4; log=[]
def P(s=''): print(s,flush=True); log.append(s)
CF=json.load(open('/workspace/out/q12b/clumpfit.json'))['HCO+|Xb4f']['V0']
x=(np.arange(NG)-NG//2)*PX; X,Y=np.meshgrid(x,x)
def sky_pos(p): return ((REGIONS[p][0]-CEN[0])*COS*3600,(REGIONS[p][1]-CEN[1])*3600)
def gauss(xc,yc,fw,flux): s=fw/FW; return flux*PX**2/(2*np.pi*s*s)*np.exp(-((X-xc)**2+(Y-yc)**2)/(2*s*s))    # Jy/pixel (flux = total, Jy km/s)
def true_sky(fext,wext):
    I=gauss(*sky_pos('N'),1.0,CF['N'])+gauss(*sky_pos('S'),1.0,CF['S'])+gauss(*sky_pos('Bridge'),2.0,CF['Bridge'])+gauss(0,0,5,CF.get('haloA',0))+gauss(0,0,14,CF['haloB'])
    for c in ('E1','E2','W'): I=I+gauss(*sky_pos(c),2.5,CF[c])
    if fext>0: I=I+gauss(0,0,wext,fext*CF['nuc_total'])
    return I
def pbmap(fw): return np.exp(-4*np.log(2)*(X**2+Y**2)/fw**2)
def uvmask(Lmin,Lmax,lam,p=1.0,M=400000,seed=0):
    rng=np.random.default_rng(seed); u=rng.uniform(size=M)
    Lr=(Lmin**(1-p)+u*(Lmax**(1-p)-Lmin**(1-p)))**(1/(1-p)) if p!=1 else Lmin*(Lmax/Lmin)**u
    th=rng.uniform(0,2*np.pi,M); q=Lr/lam/206265.0                       # cycles per arcsec
    du=1/(NG*PX); iu=np.rint(q*np.cos(th)/du).astype(int)%NG; iv=np.rint(q*np.sin(th)/du).astype(int)%NG
    W=np.zeros((NG,NG)); np.add.at(W,(iv,iu),1); W2=np.zeros((NG,NG)); np.add.at(W2,((-iv)%NG,(-iu)%NG),1); return (W+W2)/ (W+W2).sum()*M
def gfit(beam):
    c=NG//2; sub=beam[c-6:c+7,c-6:c+7]; yy,xx=np.mgrid[-6:7,-6:7]
    f=lambda xy,a,sx,sy:a*np.exp(-0.5*(xy[0]**2/sx**2+xy[1]**2/sy**2)); p,_=curve_fit(f,(xx.ravel(),yy.ravel()),sub.ravel(),p0=[1,2,2]); return abs(p[1])*PX,abs(p[2])*PX
def hogbom(dirty,beam,n=3000,gain=0.1,thr=0.015):
    res=dirty.copy(); comp=np.zeros_like(dirty); c=NG//2; pk0=np.abs(dirty).max()
    for _ in range(n):
        i=np.unravel_index(np.argmax(np.abs(res)),res.shape); v=res[i]
        if abs(v)<thr*pk0: break
        comp[i]+=gain*v; sh=np.roll(np.roll(beam,i[0]-c,0),i[1]-c,1); res-=gain*v*sh
    return comp,res
def smooth10(img,fwhm_in_sky=None): s=10/FW/PX; k=np.exp(-0.5*(np.arange(-40,41)/s)**2); k/=k.sum(); return fftconvolve(img,k[:,None]*k[None],mode='same')
def observe(sky,Lmin,Lmax,lam,pbfw,p=1.0,thr=0.015,n=3000):
    pb=pbmap(pbfw); W=uvmask(Lmin,Lmax,lam,p=p); ft=np.fft.fft2(np.fft.ifftshift(sky*pb))
    dirty=np.fft.fftshift(np.fft.ifft2(ft*W).real); beam=np.fft.fftshift(np.fft.ifft2(W).real); nb=beam.max(); dirty/=nb; beam/=nb
    comp,res=hogbom(dirty,beam,n=n,thr=thr); sa,sb=gfit(beam); 
    kern_x=np.exp(-0.5*(np.arange(-30,31)*PX/sa)**2); kern_y=np.exp(-0.5*(np.arange(-30,31)*PX/sb)**2); K=np.outer(kern_y,kern_x)
    img=(fftconvolve(comp,K,mode='same')+res)/(2*np.pi*(sa/PX)*(sb/PX))   # Jy/beam (unit-peak clean beam; comps in Jy) -> Jy per pixel (surface brightness)
    return img/np.where(pb>0.05,pb,np.nan),(sa,sb)
POS=['N','S','Bridge','E1','E2','W']
RR=np.hypot(X,Y)
def bowl_metrics(img,pbf,sky):
    """same recipe as s58b on the 10\" simulated pbcor image: bowl/peak, F_edge/F_max (r<25\" circle of the PB-uncorrected map), aperture loss r<8\" (obs/true)"""
    m=smooth10(np.nan_to_num(img)); u=m*pbmap(pbf); edges=np.arange(0,26,1.0)
    prof=np.array([m[(RR>=a)&(RR<b)].mean() for a,b in zip(edges[:-1],edges[1:])]); up=np.array([u[(RR>=a)&(RR<b)].mean() for a,b in zip(edges[:-1],edges[1:])])
    cog=np.array([u[RR<=r].sum() for r in edges[1:]]); pk=m[RR<4].max(); neg=np.where(up<0)[0]; i0=neg[0] if len(neg) else len(up)-1
    j=i0+np.argmin(prof[i0:]); tr=smooth10(sky)
    return dict(bowl_over_peak=float(prof[j]/pk),R_bowl=float(0.5*(edges[j]+edges[j+1])),Fedge_over_Fmax=float(cog[-1]/cog.max()),R0=float(edges[i0]),loss8=float(1-m[RR<8].sum()/tr[RR<8].sum()),loss25=float(1-m[RR<25].sum()/tr[RR<25].sum()))
def read10(img):
    m=smooth10(np.nan_to_num(img)); return {p:float(m[int(round(sky_pos(p)[1]/PX))+NG//2,int(round(sky_pos(p)[0]/PX))+NG//2]) for p in POS}
LAMS={'B6':2.998e8/267.56e9,'B5':2.998e8/178.375e9,'B7':2.998e8/356.73e9}
SCEN=[('no extended',0,40),('ext 25" f=0.3',0.3,25),('ext 40" f=0.5',0.5,40),('ext 80" f=1.0',1.0,80)]
OUT={}
for sn,fe,we in SCEN:
    sky=true_sky(fe,we); tr=read10(sky); R={}; BM={}
    arrs=[('Xb4b',60,400,'B6',21.9),('Xb4d',29,250,'B6',21.9),('Xb4f',8.7,45,'B6',37.6),('Xb55',8.7,45,'B5',56.4),('Xb47',8.7,45,'B7',28.2)]
    for a,Lmin,Lmax,b,pbf in arrs:
        img,bm=observe(sky,Lmin,Lmax,LAMS[b],pbf); o=read10(img); R[a]={p:o[p]/tr[p] for p in POS}; BM[a]=bowl_metrics(img,pbf,sky)
    OUT[sn]=R; OUT[sn+'|bowl']=BM
    P(f'\n### sky: {sn}   (true 10\" values {", ".join(f"{p} {tr[p]/PX**2*100:.2f}" for p in POS)} Jy km/s per 100 arcsec^2)')
    P('absolute recovery R = observed(10")/true(10"):'); 
    for a in R: P(f'  {a:5s} '+'  '.join(f'{p} {R[a][p]:.2f}' for p in POS))
    P('bowl diagnostics of the simulated 10" pbcor maps (same recipe as s58b) and the TRUE aperture loss (obs/true, r<8" / r<25"):')
    for a in BM: P(f'  {a:5s} bowl/peak {BM[a]["bowl_over_peak"]:+.3f} at r={BM[a]["R_bowl"]:.0f}"  zero crossing {BM[a]["R0"]:.0f}"  F_edge/F_max {BM[a]["Fedge_over_Fmax"]:.2f}   true loss r<8" {BM[a]["loss8"]:+.2f}  r<25" {BM[a]["loss25"]:+.2f}')
    P('slopes alpha=ln(R_hi/R_lo)/ln(MRS ratio) with nominal MRS ratios (b/d=Lmin_d/Lmin_b, d/f): ')
    mrs={'Xb4b':8.7/60,'Xb4d':8.7/29,'Xb4f':1.0}
    for p in POS:
        a_df=np.log(R['Xb4d'][p]/R['Xb4f'][p])/np.log(mrs['Xb4d']); a_bd=np.log(R['Xb4b'][p]/R['Xb4d'][p])/np.log(mrs['Xb4b']/mrs['Xb4d'])
        lam_ratio=1.5; r55=R['Xb55'][p]/R['Xb4f'][p]; r47=R['Xb47'][p]/R['Xb4f'][p]
        pred55=(1.5)**max(a_df,0); pred47=(0.75)**max(a_df,0)
        P(f'  {p:7s} alpha_df={a_df:5.2f} alpha_bd={a_bd:5.2f}   R55/R4f sim {r55:.2f} vs power-law(e=1,sat=1,kJ=1) {pred55:.2f};  R47/R4f sim {r47:.2f} vs {pred47:.2f}')

P('\n## sensitivity of the simulated relative recoveries to the clean depth (stop threshold, % of the dirty peak; <= 6000 it) and the uv radial density index p (pdf ~ L^-p); sky = ext 40" f=0.5')
sky=true_sky(0.5,40); tr=read10(sky); SENS={}
for p_,thr,nit in ((1.0,0.05,3000),(1.0,0.015,3000),(1.0,0.003,6000),(0.5,0.015,3000),(1.5,0.015,3000),(1.0,0.0,0)):
    R={}
    for a,Lmin,Lmax,b,pbf in arrs:
        img,bm=observe(sky,Lmin,Lmax,LAMS[b],pbf,p=p_,thr=thr,n=nit); o=read10(img); R[a]={p:o[p]/tr[p] for p in POS}
    SENS[f'p{p_}_thr{thr}']=R
    P(f'p={p_} thr={thr*100:.1f}%{" (no clean: dirty image)" if nit==0 else ""}: '+'; '.join(f'{pp}: R4f {R["Xb4f"][pp]:.2f} R55/R4f {R["Xb55"][pp]/R["Xb4f"][pp]:.2f} R47/R4f {R["Xb47"][pp]/R["Xb4f"][pp]:.2f} d/f {R["Xb4d"][pp]/R["Xb4f"][pp]:.2f} b/d {R["Xb4b"][pp]/R["Xb4d"][pp]:.2f}' for pp in ('N','E1')))
OUT['sens']=SENS
json.dump(OUT,open('/workspace/out/q12b/sim.json','w')); open('/workspace/out/q12b/sim.txt','w').write('\n'.join(log))
