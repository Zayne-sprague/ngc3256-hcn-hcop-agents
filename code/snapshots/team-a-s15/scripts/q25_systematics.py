# Q25 fixes: (1) zero-level (bowl) offset scenarios on the 8-12/4-8 factor; (2) integration-window test + annulus spectra.
import sys; sys.path.insert(0,'/workspace/scripts')
from q25_ratio10 import *
RB=(22,35)   # emission-free reference ring for bias (arcsec), PB>0.5 at 88 GHz extends to ~35"
def prep(path,f0):
    h=fits.open(path,memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C
    dv=abs(np.median(np.diff(v))); return h,hd,v,dv
def m0win(h,v,dv,lo,hi,sl=None):
    m=np.where((v>=lo)&(v<=hi))[0]; return np.nansum(np.array(h[0].data[0,m[0]:m[-1]+1],float),axis=0)*dv
def sm_(d,a): return sm(d,a)
def stats(d,s):
    r=np.hypot(d['dx']+0.075,d['dy']); omp=OM(tb)/d['pix']**2
    N=[((r>=a)&(r<b)).sum() for a,b in ANN]; fl=np.array([np.nansum(s[(r>=a)&(r<b)]) for a,b in ANN])/omp
    ring=(r>=RB[0])&(r<RB[1]); return fl,np.array(N)/omp,np.nanmean(s[ring]),np.nanstd(s[ring])
out={}
print('=== (1) bowl: mean of 10" moment-0 in r=22-35" ring (Jy/beam km/s) and offset scenarios')
for nm,p in P.items():
    for ln,f0 in REST.items():
        d=load10(p,f0); s=sm(d); fl,nb,mu,sd=stats(d,s)
        nl=[stats(d,sm(d,n))[2] for n in d['nulls']]
        print(f'{nm} {ln}: ring mean {mu:.3f} (std {sd:.3f}); line-free-window ring means {np.round(nl,3)}')
        out[(nm,ln)]=(fl,nb,mu)
for k,d,s_ in [(k,load(k,nullspec=False),None) for k in ('HCN 2-1','HCO+ 2-1')]:
    s=smooth(d,10.); fl,nb,mu,sd=stats(d,s); print(f'Xb55 {k}: ring mean {mu:.3f} (std {sd:.3f}) [PB>0.5 at 165-178 GHz ends ~r=17-18"! ring is beyond PB half-power]'); out[('Xb55',k)]=(fl,nb,mu)
def ratio_scen(nm,lh,lo_,frac):
    fh,nh,mh=out[(nm,lh)]; fo,no,mo=out[(nm,lo_)]
    bh=-min(mh,0)*frac; bo=-min(mo,0)*frac; F1=fh+bh*nh; F2=fo+bo*no; R=F1/F2; return R
print('\nfraction of the full offset (b = -ring mean, i.e. no real emission at 22-35") added back: 8-12/4-8 factor and R(0-4,4-8,8-12)')
for nm,a,b in (('Michiyama','HCN 1-0','HCO+ 1-0'),('Harada','HCN 1-0','HCO+ 1-0'),('Xb55','HCN 2-1','HCO+ 2-1')):
    for fr in (0,0.5,1.0):
        R=ratio_scen(nm,a,b,fr); print(f'{nm:9s} frac {fr}: R {R.round(3)} factor {R[3]/R[1]:.3f}')
# (2) windows
print('\n=== (2) integration window test (half width km/s), 8-12/4-8 factor and R')
for nm,p in P.items():
    H={}
    for ln,f0 in REST.items():
        h,hd,v,dv=prep(p,f0); d=load10(p,f0)
        for hw in (250,350,500):
            m0=m0win(h,v,dv,VSYS-hw,VSYS+hw); dd=dict(d); dd['m0']=m0; H[(ln,hw)]=stats(dd,sm(dd))[0]
    for hw in (250,350,500):
        R=H[('HCN 1-0',hw)]/H[('HCO+ 1-0',hw)]; print(f'{nm} +-{hw}: R {R.round(3)} factor {R[3]/R[1]:.3f}')
for hw in (250,350,500):
    H=[]
    for k in ('HCN 2-1','HCO+ 2-1'):
        d=load(k,nullspec=False,window=(VSYS-hw,VSYS+hw)); H.append(stats(d,smooth(d,10.))[0])
    R=H[0]/H[1]; print(f'Xb55 +-{hw}: R {R.round(3)} factor {R[3]/R[1]:.3f}')
# annulus spectra (10" beam), Michiyama
print('\nannulus-averaged spectra, Michiyama 1-0 (10" smoothed), mean Jy/beam: v, 4-8 HCN,HCO+, 8-12 HCN,HCO+, ring22-35 HCN,HCO+')
S={}
for ln,f0 in REST.items():
    h,hd,v,dv=prep(P['Michiyama'],f0); d=load10(P['Michiyama'],f0); r=np.hypot(d['dx']+0.075,d['dy'])
    k=kernel(d['hd'],tb,d['beam']); sel=np.where((v>VSYS-1000)&(v<VSYS+1000))[0]
    for i in sel:
        a=np.nan_to_num(np.array(h[0].data[0,i],float)); s=fftconvolve(a,k,mode='same')*OM(tb)/OM(d['beam'])
        S.setdefault(ln,[]).append((v[i],[np.nanmean(s[(r>=x)&(r<y)]) for x,y in ((4,8),(8,12),(22,35))]))
for (vh,a),(vo,b) in zip(S['HCN 1-0'],S['HCO+ 1-0']): print(f'{vh:7.0f} {a[0]*1e3:7.2f} {a[1]*1e3:7.2f} {a[2]*1e3:7.2f} | {b[0]*1e3:7.2f} {b[1]*1e3:7.2f} {b[2]*1e3:7.2f}')
