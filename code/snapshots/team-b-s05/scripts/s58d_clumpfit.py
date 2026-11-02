# Q12 round 2, step d: de-contaminated clump fluxes. Forward-model each array's native-beam moment-0 map (|v|<350 km/s, +-19" box about the N-S midpoint) as a sum of
# Gaussian sky components convolved with THAT array's (elliptical) beam: N core, S core, Bridge, two nuclear halos (FWHM 5" and 14", at the midpoint), the clumps E1, E2, W
# (Gaussian, FWHM W_c) and a free-sign constant (the bowl/zero level).  Linear in the amplitudes -> NNLS (const free sign).  The fitted TOTAL clump flux (Jy km/s) is
# beam-independent, so arrays and lines can be compared directly (no point value at a position).  Errors: (i) 16 sign-flip noise realisations refitted, (ii) spread over the
# model variants (clump FWHM 1.5/2.5/4", halo-B 10/14/20", no halo-A, no constant).  Sky-component positions are fixed (Q0 regions).
# usage: python3 -I s58d_clumpfit.py -> out/q12b/clumpfit.json, clumpfit.txt
import sys,json; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from scipy.optimize import nnls
from common import *
from regions import REGIONS
import s56_lib as L
CEN=L.CEN; COS=np.cos(np.deg2rad(CEN[1])); log=[]
def P(s=''): print(s,flush=True); log.append(s)
FW=2*np.sqrt(2*np.log(2))
def getmap(ous,line,J):
    """-> dx,dy (arcsec E,N of CEN) pixel grids, m0, noise stack, beam (maj,min,pa), flattened to the +-19" box"""
    if ous in ('Xb4b','Xb4d') or (ous=='Xb4f' and False):
        a=L.Arr(ous,line); m0,nz=a.m0,a.noise; n=m0.shape[0]; jj=np.arange(n); Y,X=np.mgrid[:n,:n]; xp=a.x0+X*a.k+(a.k-1)/2; yp=a.y0+Y*a.k+(a.k-1)/2; wc=a.wc; b=a.b0
    else:
        d=np.load(f'/workspace/out/q12b/full_{ous}_{line}.npz'); m0=d['m0']; nz=d['noise']; n=m0.shape[0]; Y,X=np.mgrid[:n,:n]; xp,yp=X,Y; b=tuple(d['beam'])
        spw=dict(Xb55={'HCN':'spw16','HCO+':'spw18'},Xb4f={'HCO+':'spw16','HCN':'spw18'},Xb47={'HCN':'spw16','HCO+':'spw18'})[ous][line]
        wc=WCS(fits.getheader(files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0])).celestial
    ra,dec=wc.pixel_to_world_values(xp,yp); dx=(ra-CEN[0])*COS*3600; dy=(dec-CEN[1])*3600
    ok=(abs(dx)<19)&(abs(dy)<19)
    return dx[ok],dy[ok],m0[ok],np.array([q[ok] for q in nz]),b
def beamcov(b):
    sa,sb,pa=b[0]/FW,b[1]/FW,np.deg2rad(b[2]); u=np.array([np.sin(pa),np.cos(pa)]); v=np.array([np.cos(pa),-np.sin(pa)])
    return sa**2*np.outer(u,u)+sb**2*np.outer(v,v),2*np.pi*sa*sb
def tmpl(dx,dy,xc,yc,fwhm,b):
    Sb,Om=beamcov(b); S=Sb+np.eye(2)*(fwhm/FW)**2; Si=np.linalg.inv(S); x=dx-xc; y=dy-yc
    q=Si[0,0]*x*x+2*Si[0,1]*x*y+Si[1,1]*y*y; return Om/(2*np.pi*np.sqrt(np.linalg.det(S)))*np.exp(-0.5*q)
def sky(p): return ((REGIONS[p][0]-CEN[0])*COS*3600,(REGIONS[p][1]-CEN[1])*3600)
VARS={'V0':dict(wc=2.5,hb=14,ha=5,const=True),'Wc1.5':dict(wc=1.5,hb=14,ha=5,const=True),'Wc4':dict(wc=4,hb=14,ha=5,const=True),'hB10':dict(wc=2.5,hb=10,ha=5,const=True),
      'hB20':dict(wc=2.5,hb=20,ha=5,const=True),'noA':dict(wc=2.5,hb=14,ha=None,const=True),'noconst':dict(wc=2.5,hb=14,ha=5,const=False)}
NAMES=['N','S','Bridge','haloA','haloB','E1','E2','W']
def design(dx,dy,b,v):
    cols={}
    cols['N']=tmpl(dx,dy,*sky('N'),1.0,b); cols['S']=tmpl(dx,dy,*sky('S'),1.0,b); cols['Bridge']=tmpl(dx,dy,*sky('Bridge'),2.0,b)
    if v['ha']: cols['haloA']=tmpl(dx,dy,0,0,v['ha'],b)
    cols['haloB']=tmpl(dx,dy,0,0,v['hb'],b)
    for c in ('E1','E2','W'): cols[c]=tmpl(dx,dy,*sky(c),v['wc'],b)
    names=list(cols); A=np.array([cols[k].ravel() for k in names]).T
    if v['const']: A=np.hstack([A,np.ones((A.shape[0],1)),-np.ones((A.shape[0],1))])
    return names,A
def fit(A,y):
    x,_=nnls(A,y.ravel()); k=A.shape[1]
    return x
if __name__=='__main__':
    RES={}
    JOBS=[(o,l,J) for l in ('HCO+','HCN') for o,J in (('Xb4b',3),('Xb4d',3),('Xb4f',3),('Xb55',2),('Xb47',4))]
    for ous,line,J in JOBS:
        dx,dy,m0,nz,b=getmap(ous,line,J); pk=m0.max()
        out={}
        for vn,v in VARS.items():
            names,A=design(dx,dy,b,v); x=fit(A,m0); mod=A@x; res=m0.ravel()-mod
            # fit amplitudes are total fluxes (Jy km/s); const gives bowl offset
            d={n:float(x[i]) for i,n in enumerate(names)}; d['const']=float(x[len(names)]-x[len(names)+1]) if v['const'] else 0.0
            d['rms_over_peak']=float(res.std()/pk); d['nuc_total']=float(sum(d[n] for n in ('N','S','Bridge','haloA','haloB') if n in d))
            if vn=='V0':
                ens=[fit(A,m0+q) for q in nz-nz.mean(0)] if False else [fit(A,q) for q in nz]
                d['err']={n:float(np.std([e[i] for e in ens])) for i,n in enumerate(names)}
            out[vn]=d
        RES[f'{line}|{ous}']=out
        s='  '.join(f'{c} {out["V0"][c]:.2f}+-{out["V0"]["err"][c]:.2f}(sys {np.std([out[v][c] for v in out]):.2f})' for c in ('E1','E2','W'))
        P(f'{line:5s}{ous} beam {b[0]:.2f}x{b[1]:.2f}  clump total flux [Jy km/s] {s}; N {out["V0"]["N"]:.1f} S {out["V0"]["S"]:.1f} nuc_total {out["V0"]["nuc_total"]:.1f} const {out["V0"]["const"]:.2f}; resid rms/peak {out["V0"]["rms_over_peak"]:.3f}')
    json.dump(RES,open('/workspace/out/q12b/clumpfit.json','w'),indent=1); open('/workspace/out/q12b/clumpfit.txt','w').write('\n'.join(log))
