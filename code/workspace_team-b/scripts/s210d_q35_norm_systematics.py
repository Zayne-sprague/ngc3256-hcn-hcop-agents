# Q35 fix 3+4: (a) alternative normalisers at the matched 1.7" beam (Harada+2018 Table 13): C18O 1-0/2-1, 13CO 1-0/2-1, CS 2-1/5-4 and the thin-tracer excitation N/S (C18O 2-1/1-0, 13CO 2-1/1-0);
# (b) free-free contribution to the 346 GHz continuum from H40alpha (99 GHz; RRL/continuum LTE relation, Rohlfs&Wilson 14.? : int T_L dv / T_C = 6.985e3 nu_GHz^1.1 Te^-1.15 /(1+y) km/s, Te 5000-10000 K, y=0.08, alpha_ff=-0.1);
# (c) velocity-resolved deblend of HCN 4-3, HCO+ 4-3 and CS 7-6 (Xb47, N/S/Bridge + constant, 4.31x2.79" beam, 100 km/s bins) -> line profile shapes / FWHM / dip at N and S and HCN/HCO+ vs velocity;
# (d) H13CN/H13CO+ 1-0 thin-isotopologue N/S with a proper posterior (S truncated at >=0), with 12C/13C and fractionation marginalised.
# usage: python3 -I s210d_q35_norm_systematics.py -> out/q35/s210d.txt/json
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, pandas as pd, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits; from astropy.wcs import WCS; import astropy.units as u
from astropy.convolution import convolve_fft; from radio_beam import Beam
from regions import REGIONS
from common import files
O='/workspace/out/q35/'; T=[]; R={}
def P(s=''): print(s,flush=True); T.append(s)
H=pd.read_csv('/workspace/lit/harada_table13_parsed.csv'); H['sp']=H.species_transition.str.replace(' ','')
def g(p,s): r=H[(H.pos==p)&(H.sp.str.startswith(s))].iloc[0]; return r.I_Kkms,r.eI_Kkms
rng=np.random.default_rng(5)
def draw(p,s,n=20000,sys=0.07): I,e=g(p,s); return rng.normal(I,np.hypot(e,sys*I),n)    # table rms + 7% flux scale (per line)
P('== (a) N/S double ratios against alternative normalisers (Harada 1.7" peak pixels, K km/s; 16-84% of MC; flux-scale 7% per line)')
P('  raw line N/S: '+'; '.join(f"{s} {np.median(draw('N',s)/draw('S',s)):.2f}" for s in ('HCN1','HCO+1','HCN3','HCO+3','CS2','CS5','C18O1','C18O2','13CO1','13CO2')))
pairs=[('HCN1','HCO+1','C18O1','1-0 / C18O 1-0'),('HCN1','HCO+1','13CO1','1-0 / 13CO 1-0'),('HCN1','HCO+1','CS2','1-0 / CS 2-1'),('HCN3','HCO+3','C18O2','3-2 / C18O 2-1'),('HCN3','HCO+3','13CO2','3-2 / 13CO 2-1'),('HCN3','HCO+3','CS5','3-2 / CS 5-4'),('HCN3','HCO+3','C18O1','3-2 / C18O 1-0 (J-mismatch)'),('HCN3','HCO+3','13CO1','3-2 / 13CO 1-0 (J-mismatch)')]
for a,b,c,lab in pairs:
    o=[]
    for nm,x in (('HCN',a),('HCO+',b)):
        v=(draw('N',x)/draw('N',c))/(draw('S',x)/draw('S',c)); q=np.percentile(v,[16,50,84]); o.append(f'{nm}/norm N/S = {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}]'); R[f'norm|{lab}|{nm}']=q.tolist()
    P(f'  {lab:32s}: '+'; '.join(o))
P('  thin-tracer excitation (K km/s ratio): '+'; '.join(f"{lab} N {g('N',a)[0]/g('N',b)[0]:.2f} S {g('S',a)[0]/g('S',b)[0]:.2f} -> N/S {np.median((draw('N',a)/draw('N',b))/(draw('S',a)/draw('S',b))):.2f} [{np.percentile((draw('N',a)/draw('N',b))/(draw('S',a)/draw('S',b)),16):.2f},{np.percentile((draw('N',a)/draw('N',b))/(draw('S',a)/draw('S',b)),84):.2f}]" for lab,a,b in (('C18O 2-1/1-0','C18O2','C18O1'),('13CO 2-1/1-0','13CO2','13CO1'),('CS 5-4/2-1','CS5','CS2'),('HCN 3-2/1-0','HCN3','HCN1'),('HCO+ 3-2/1-0','HCO+3','HCO+1'))))
P('  (C18O/13CO 2-1/1-0 are optically thinner, abundance-independent excitation tracers: both are higher at N, like CS 5-4/2-1; HCN and HCO+ r31 are not.)')
# (b) free-free
P('\n== (b) free-free at 346 GHz from H40alpha (99.02 GHz), 1.7" peak pixels')
for p in ('N','S','C'):
    I,e=g(p,'H40'); out=[]
    for Te in (5000,10000):
        k=6.985e3*99.02**1.1*Te**-1.15/1.08; Tc=I/k; S99=Tc*99.02**2*1.7**2/1222.   # mJy/beam (1.7" circular)
        S346=S99*(346./99.02)**-0.1; out.append(S346)
    R[f'ff|{p}']=out; P(f'  {p}: H40a {I:.2f}+-{e:.2f} K km/s -> S_ff(346 GHz, 1.7" beam, point-like) = {out[1]:.1f} (Te=1e4) - {out[0]:.1f} (Te=5e3) mJy/beam')
c={'N':39.3,'S':21.8,'Bridge':11.0}; ce=1.8
for lab,fN,fS in (('point-like ff (=1.7" peak flux, all in the 4.3x2.8" beam)',R['ff|N'][1],R['ff|S'][1]),('upper: Te=5e3 point-like',R['ff|N'][0],R['ff|S'][0]),('extended ff, 2x (flux in 4.3x2.8 beam twice the 1.7 peak; generous)',2*R['ff|N'][1],2*R['ff|S'][1])):
    cn=c['N']-fN; cs_=c['S']-fS; f=(c['N']/c['S'])
    P(f'  {lab}: dust-only continuum N {cn:.1f} S {cs_:.1f} mJy/beam; N/S = {cn/cs_:.2f} (uncorrected {f:.2f}) -> HCO+/cont N/S {0.67*f/(cn/cs_):.2f}, HCN/cont {1.09*f/(cn/cs_):.2f}, CS/cont {0.99*f/(cn/cs_):.2f}')
# (c) velocity resolved deblend
C=299792.458; tb=Beam(4.3122*u.arcsec,2.7884*u.arcsec,86.676*u.deg)
LINES={'CS 7-6':('spw22',342.88285),'HCN 4-3':('spw16',354.50547),'HCO+ 4-3':('spw18',356.73424)}
def bins(spw,f0,vc=2780.,edges=np.arange(-700,701,100)):
    h=fits.open(files(f'*Xb47*{spw}.cube.I.pbcor.fits')[0]); hd=h[0].header; d=np.nan_to_num(h[0].data[0]); nz=d.shape[0]; ps=abs(hd['CDELT2'])*3600
    nu=hd['CRVAL3']+(np.arange(nz)+1-hd['CRPIX3'])*hd['CDELT3']; v=C*(1-nu/(f0*1e9)); dv=abs(np.median(np.diff(v)))
    wc=WCS(hd).celestial; cx,cy=[float(t) for t in wc.world_to_pixel_values(156.96343,-43.90461)]; hw=int(14/ps); x0,x1,y0,y1=int(cx-hw),int(cx+hw),int(cy-hw),int(cy+hw)
    b=Beam(hd['BMAJ']*u.deg,hd['BMIN']*u.deg,hd['BPA']*u.deg); conv=(tb.major>b.major or tb.minor>b.minor)
    if conv: ker=tb.deconvolve(b).as_kernel(ps*u.arcsec,x_size=41,y_size=41).array; ker/=ker.sum(); fac=(tb.sr/b.sr).value
    sel=np.where(abs(v-vc)<1100)[0]; dd=np.array([(convolve_fft(d[i][y0:y1,x0:x1],ker,boundary='fill',allow_huge=True)*fac) if conv else d[i][y0:y1,x0:x1] for i in sel]); vs=v[sel]
    sb=abs(vs-vc)>500; X=np.vstack([np.ones(sb.sum()),(vs[sb]-vc)/600.]).T; Pm=np.linalg.pinv(X); coef=np.tensordot(Pm,dd[sb],axes=(1,0)); dd=dd-(coef[0][None]+coef[1][None]*((vs-vc)/600.)[:,None,None])
    out=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        m=(vs-vc>=lo)&(vs-vc<hi); out.append(dd[m].sum(0)*dv)
    return np.array(out),(wc,x0,y0,ps),vs,sb,dd
yy=None; amps={}
for ln,(spw,f0) in LINES.items():
    mb,(wc,x0,y0,ps),vs,sb,dd=bins(spw,f0); ny,nx=mb.shape[1:]; yy,xx=np.mgrid[:ny,:nx]
    pos={k:[float(t)-o for t,o in zip(wc.world_to_pixel_values(REGIONS[k][0],REGIONS[k][1]),(x0,y0))] for k in ('N','S','Bridge')}
    bm=tb.major.value/ps; bn=tb.minor.value/ps; th=np.deg2rad(tb.pa.value+90)
    def gb(px,py):
        dx=xx-px; dy=yy-py; xp=dx*np.cos(th)+dy*np.sin(th); yp=-dx*np.sin(th)+dy*np.cos(th); return np.exp(-0.5*((xp/(bm/2.355))**2+(yp/(bn/2.355))**2))
    mx,my=(pos['N'][0]+pos['S'][0])/2,(pos['N'][1]+pos['S'][1])/2; box=(abs(xx-mx)*ps<8)&(abs(yy-my)*ps<8)
    A=np.array([gb(*pos[c])[box] for c in ('N','S','Bridge')]+[np.ones(box.sum())]).T; Ai=np.linalg.pinv(A)
    # per-bin errors: scatter of the same fit over off-source baseline-only 100 km/s windows would need line-free bins; use the box residual std x sqrt(8)/sqrt(npix_per_beam-ish) as in s98 (errors x sqrt8)
    sol=np.array([Ai@m[box] for m in mb]); res=np.array([np.std(m[box]-A@s) for m,s in zip(mb,sol)]); cov=np.linalg.inv(A.T@A)[np.arange(3),np.arange(3)]
    err=np.sqrt(cov)[None,:]*res[:,None]*np.sqrt(8)
    amps[ln]=(sol[:,:3],err)
edges=np.arange(-700,701,100); cen=(edges[:-1]+edges[1:])/2
P('\n== (c) velocity-resolved deblended amplitudes (Jy/beam km/s per 100 km/s bin, 4.31x2.79" beam; v relative to 2780 km/s) and HCN/HCO+ 4-3 (K ratio)')
fk={'CS 7-6':339.7,'HCN 4-3':351.2,'HCO+ 4-3':353.4}
for ci_,comp in enumerate(('N','S')):
    P(f' {comp}: v    HCN        HCO+       CS7-6      HCN/HCO+(K)')
    for j,v0 in enumerate(cen):
        a,ea=amps['HCN 4-3'][0][j,ci_],amps['HCN 4-3'][1][j,ci_]; b,eb=amps['HCO+ 4-3'][0][j,ci_],amps['HCO+ 4-3'][1][j,ci_]; c_,ec=amps['CS 7-6'][0][j,ci_],amps['CS 7-6'][1][j,ci_]
        r=(a/b)*(fk['HCO+ 4-3']/fk['HCN 4-3'])**2 if b>2*eb and a>0 else float('nan')
        P(f'   {v0:+5.0f}  {a:6.2f}+-{ea:.2f} {b:6.2f}+-{eb:.2f} {c_:6.2f}+-{ec:.2f}  {r:5.2f}')
    for ln in LINES:
        s=amps[ln][0][:,ci_]; 
        if s.max()>0:
            hm=s.max()/2; idx=np.where(s>=hm)[0]; cc=(s*cen).sum()/s.sum()
            w=np.interp(hm,s[:idx[0]+1][::1],cen[:idx[0]+1]) if False else None
            R[f'prof|{comp}|{ln}']=dict(peak=float(s.max()),vpeak=float(cen[s.argmax()]),fwhm_bins=float(len(idx)*100),centroid=float(cc),sum=float(s.sum()))
            P(f'   {comp} {ln}: peak bin {s.max():.2f} at {cen[s.argmax()]:+.0f}; bins >= half max span {len(idx)*100} km/s; centroid {cc:+.0f}; ratio peak-bin/sum {s.max()/s.sum():.2f}')
P('  Peak/total of the HCO+ profile at N vs S (a self-absorption dip would show as a lower peak/total and a minimum between two bins).')
# (d) thin isotopologues: H13CN/H13CO+ 1-0 (Q29, 2.2", K km/s)  N 1.88+-0.31 / 2.50+-0.34; S 0.31+-0.32 / 1.66+-0.36; Bridge 1.80+-0.34/1.92+-0.35
n=200000
def cn(a,ea,b,eb,trunc=True):
    x=rng.normal(a,ea,n); y=rng.normal(b,eb,n); 
    if trunc: x=np.where(x<0,np.nan,x)
    return x/y
rN=cn(1.88,.31,2.50,.34); rS=cn(0.31,.32,1.66,.36); rB=cn(1.80,.34,1.92,.35)
# 12C/13C 30-100 common for HCN and HCO+ (cancels) ; independent fractionation/isotopic-chemistry factor per species: log-normal sigma 0.3 dex? -> use 0.15 dex per species (x1.4)
frac=lambda: 10**rng.normal(0,0.15,n)
f=2.9*10**rng.normal(0,0.1,n)   # RADEX conversion N(H13CN)/N(H13CO+) = f x intensity ratio (s188 table: 1.9-4.4 at n=1e4-1e5, Tk 25-100)
R['thin']={}
for nm,r in (('N',rN),('S',rS),('Bridge',rB)):
    q=np.nanpercentile(r,[16,50,84]); q3=np.nanpercentile(r,99.87); R['thin'][nm]=q.tolist()
    P(f'\n  H13CN/H13CO+ 1-0 intensity ratio {nm}: median {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}] (S truncated >=0: 3-sigma upper {q3:.2f})')
x=rN/rS; ok=np.isfinite(x); q=np.percentile(x[ok],[16,50,84]); R['thin']['N/S']=q.tolist()
P(f'  N/S of H13CN/H13CO+ (intensity) = {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}] (1-sigma range), P(N/S>1) = {(x[ok]>1).mean():.2f}, 5th percentile {np.percentile(x[ok],5):.2f}; NOT significant (S is a {0.31/0.32:.1f}-sigma H13CN value)')
xc=(rN*frac()/frac())/(rS*frac()/frac()); ok=np.isfinite(xc); q=np.percentile(xc[ok],[16,50,84]); P(f'  with independent 0.15 dex isotopic-fractionation factors per species and a common 12C/13C (30-100, cancels): N(HCN)/N(HCO+) N/S = {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}], P(>1) = {(xc[ok]>1).mean():.2f}')
for nm,r in (('N',rN),('S',rS),('Bridge',rB)):
    v=r*f; q=np.nanpercentile(v,[16,50,84]); P(f'  N(HCN)/N(HCO+) (thin, RADEX factor 2.9 x/÷1.25) {nm}: {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}]')
P('  blends in the H13CN 1-0 window: HCN ladders N/S H13CN/(HCN) not tested for SO/SiO/CH3OH here; X95d H13CN (86.34 GHz) is clear of SiO 2-1 (86.85) and the strong CH3OH lines, but H13CN 1-0 is a hyperfine triplet (F=1-1,2-1,0-1): the template integral used in Q29 includes all three.')
open(O+'s210d.txt','w').write('\n'.join(T)); json.dump(R,open(O+'s210d.json','w'),indent=1,default=float)
