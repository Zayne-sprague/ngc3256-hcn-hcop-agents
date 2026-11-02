# Reality tests of the S-nucleus blue wing and the 'Sblue' blob (Xb4d 0.94", Xb4b 0.25", HCO+ and HCN 3-2)
#  T1 empirical noise: ap-averaged blue-wing map vs pseudo-wing maps (line-free windows) at all pixels / annulus around S
#  T2 sidelobe test: Sblue spectrum vs a*S+b*N (the dirty-beam response is channel independent -> same spectral shape)
#  T3 beam-area scaling of an extended wing (Jy/beam ratio expected = beam-area ratio) in K km/s
#  T4 Xb4b convolved to the Xb4d beam
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
from scipy import ndimage, signal
import json
SETS={'HCO+':('spw25',267.55763),'HCN':('spw27',265.88618)}
FOV=12.0   # arcsec half box about midpoint of N,S
def load(ous,spw,f0,vmax=950,fov=FOV):
    fn=files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]
    h,hd,w,nu=cube_info(fn); v=vrad(nu,f0*1e9)-VSYS; ps=abs(hd['CDELT1'])*3600
    wc=WCS(hd).celestial; cx,cy=wc.world_to_pixel_values(*[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)])
    r=int(fov/ps); x0,x1=int(cx-r),int(cx+r); y0,y1=int(cy-r),int(cy+r)
    ok=np.where(abs(v)<vmax)[0]; i0,i1=ok[0],ok[-1]+1
    cube=np.array(h[0].data[0,i0:i1,y0:y1,x0:x1],dtype=np.float32)
    return cube,v[i0:i1],hd,wc,(x0,y0),ps
def disk(r):
    n=int(np.ceil(r))*2+1; yy,xx=np.mgrid[:n,:n]-n//2; k=(np.hypot(xx,yy)<=r).astype(float); return k/k.sum()
def gk(bmaj,bmin,pa,ps,n=None):   # unit-area elliptical Gaussian kernel, bmaj/bmin FWHM arcsec, pa deg E of N
    sx,sy=bmin/2.3548/ps,bmaj/2.3548/ps; n=n or int(6*bmaj/2.3548/ps)|1
    yy,xx=np.mgrid[:n,:n]-n//2; th=np.deg2rad(pa)   # major axis along PA (from N toward E; x axis points W when RA flipped -> use sky coords xE=-x)
    xe=-xx; u=xe*np.sin(th)+yy*np.cos(th); t=-xe*np.cos(th)+yy*np.sin(th)   # u along major axis, t along minor
    k=np.exp(-0.5*((u/(bmaj/2.3548/ps))**2+(t/(bmin/2.3548/ps))**2)); return k/k.sum()
def winmap(cube,v,a,b,dv):
    m=(v>=a)&(v<=b); return np.nansum(cube[m],axis=0)*dv,m.sum()
out={}; log=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); log.append(s)
for ous in ('Xb4d','Xb4b'):
  for line,(spw,f0) in SETS.items():
    cube,v,hd,wc,(x0,y0),ps=load(ous,spw,f0); dv=abs(np.median(np.diff(v)))
    bmaj,bmin=hd['BMAJ']*3600,hd['BMIN']*3600; apr=0.5*np.sqrt(bmaj*bmin)/ps
    free=(abs(v)>500)
    base=np.nanmedian(cube[free],axis=0); cube=cube-base          # per-pixel baseline from line-free channels
    K=disk(apr)
    def apmap(m): return ndimage.convolve(np.nan_to_num(m),K,mode='constant')
    wins={'blue':(-450,-150),'red':(150,450),'pseudoR':(450,750),'pseudoB':(-900,-600)}
    M={}
    for k,(a,b) in wins.items():
        m,n=winmap(cube,v,a,b,dv)
        if n>=0.8*(b-a)/dv: M[k]=apmap(m)
    pos={k:wc.world_to_pixel_values(r[0],r[1]) for k,r in REGIONS.items()}
    pos={k:(p[0]-x0,p[1]-y0) for k,p in pos.items()}
    yy,xx=np.mgrid[:cube.shape[1],:cube.shape[2]]
    def rad(k): return np.hypot(xx-pos[k][0],yy-pos[k][1])*ps
    # exclude |offset|<2" from a nucleus? keep: sidelobe environment is what matters. Use r<FOV*0.85 from centre
    cen=np.hypot(xx-cube.shape[2]/2,yy-cube.shape[1]/2)*ps<0.85*FOV
    # refined Sblue peak of the blue map within 1.2" of nominal
    d=rad('Sblue')<1.2; iy,ix=np.unravel_index(np.argmax(np.where(d,M['blue'],-1e9)),M['blue'].shape)
    P(f'## {ous} {line}: beam {bmaj:.2f}x{bmin:.2f}, ap diameter {2*apr*ps:.2f}", dv={dv:.1f}, v range {v.min():.0f}..{v.max():.0f}')
    for tag,(cx_,cy_) in {'Sblue_nominal':pos['Sblue'],'Sblue_peak':(ix,iy),'S':pos['S']}.items():
        i,j=int(round(cy_)),int(round(cx_))
        vals={k:M[k][i,j] for k in M}
        P(f'  {tag:14s} '+'  '.join(f'{k}={vals[k]:+.4f}' for k in vals))
    for k in ('pseudoR','pseudoB'):
        if k not in M: continue
        sel=cen&(rad('S')>1.5)&(rad('N')>1.5)
        x=M[k][sel]; s=1.4826*np.median(abs(x-np.median(x)))
        ann=cen&(rad('S')>1.5)&(rad('S')<6)
        xa=M[k][ann]; sa=1.4826*np.median(abs(xa-np.median(xa)))
        P(f'  {k}: field rms(MAD) {s:.4f} (n={sel.sum()}, 99.9% |x|={np.percentile(abs(x),99.9):.4f}, max {x.max():.4f}, min {x.min():.4f}); annulus 1.5-6" of S: sigma {sa:.4f}, max {xa.max():.4f}, min {xa.min():.4f}')
        out[f'{ous}_{line}_{k}']=dict(sig=s,p999=float(np.percentile(abs(x),99.9)),sig_ann=sa,max_ann=float(xa.max()))
    # percentile of Sblue blue-wing value among pseudo-wing values / blue map elsewhere
    bs=M['blue'][iy,ix]
    allp=np.concatenate([M[k][cen&(rad('S')>1.5)&(rad('N')>1.5)].ravel() for k in ('pseudoR','pseudoB') if k in M])
    P(f'  Sblue blue-wing(ap, peak pix)={bs:.4f} Jy/beam km/s; fraction of pseudo-wing pixels >= this: {(allp>=bs).mean():.2e} (n={allp.size}); N-sigma (MAD) = {bs/(1.4826*np.median(abs(allp-np.median(allp)))):.1f}')
    # --- T2 sidelobe/spectral shape test (ap spectra)
    z=np.load(f'/workspace/out/regspec/{ous}.{spw}.cube.pbcor.npz'); vv=vrad(z['nu'],f0*1e9)-VSYS
    def sp(k): 
        y=z[k+'_ap']; fl=(abs(vv)>500)&np.isfinite(y); return y-np.nanmedian(y[fl])
    sb,sS,sN=sp('Sblue'),sp('S'),sp('N')
    wing=(abs(vv)<=450)&(abs(vv)>=150); core=abs(vv)<150; fl=(abs(vv)>500)&np.isfinite(sb)
    sig=np.std(sb[fl])
    # (i) scale S by the Sblue/S ratio of the core, predict wing; (ii) fit a*S+b*N on core channels
    A=np.vstack([sS[core],sN[core]]).T; coef,*_=np.linalg.lstsq(A,sb[core],rcond=None)
    pred=coef[0]*sS+coef[1]*sN
    for nm,(a,b) in {'blue':(-450,-150),'red':(150,450)}.items():
        m=(vv>=a)&(vv<=b)
        P(f'  T2 {nm} wing: Sblue obs={sb[m].sum()*dv:+.4f}+-{sig*np.sqrt(m.sum())*dv:.4f}, sidelobe-model (a*S+b*N from core fit a={coef[0]:.3f},b={coef[1]:.3f}) predicts {pred[m].sum()*dv:+.4f}; ratio wing/core: Sblue {sb[m].sum()/sb[core].sum():.3f} vs S {sS[m].sum()/sS[core].sum():.3f} vs N {sN[m].sum()/sN[core].sum():.3f}')
    # max sidelobe level allowed: wing excess if a, b were as large as the largest plausible dirty-beam sidelobe 10%
    m=(vv>=-450)&(vv<=-150); P(f'     a 10% sidelobe-level copy of S would give a blue-wing flux {0.1*sS[m].sum()*dv:.4f}, of N {0.1*sN[m].sum()*dv:.4f} Jy/beam km/s')
    out[f'{ous}_{line}']=dict(blue_peak=float(bs),coef=[float(c) for c in coef])
    # --- T3 K km/s
    # Tb = 1.222e3 * S[Jy/beam] /(nu^2 [GHz] bmaj bmin ["]) K
    nu_obs=f0*(1-VSYS/C)
    fac=1.222e3*1e3/(nu_obs**2*bmaj*bmin)   # K per Jy/beam (1.222e3 is per mJy/beam)
    for k,(a,b) in {'blue':(-450,-150)}.items():
        P(f'  T3 {ous} {line} Sblue blue-wing = {bs*fac:.2f} K km/s ; S-core(-150..150, ap) = {sS[core].sum()*dv*fac:.2f} K km/s')
    out[f'{ous}_{line}']['blue_Kkms']=float(bs*fac)
    if ous=='Xb4b':   # T4 convolve to Xb4d beam
        tgt=dict(bmaj=0.93,bmin=0.74,pa=-86.0)
        from scipy.signal import fftconvolve
        # deconvolve beams (covariance subtraction in sky coords)
        def cov(a,b,pa):
            th=np.deg2rad(pa); R=np.array([[np.sin(th),np.cos(th)],[np.cos(th),-np.sin(th)]])  # (E,N) comps of major and minor axes
            Sg=np.diag([(a/2.3548)**2,(b/2.3548)**2]); 
            e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
        C2=cov(tgt['bmaj'],tgt['bmin'],tgt['pa'])-cov(bmaj,bmin,hd['BPA'])
        w_,U=np.linalg.eigh(C2); w_=np.clip(w_,1e-6,None)
        # kernel in pixel coords (xE = -x, y=N)
        n=int(8*np.sqrt(w_.max())/ps)|1; ky,kx=np.mgrid[:n,:n]-n//2; E=-kx*ps; Nn=ky*ps
        q=np.stack([E.ravel(),Nn.ravel()]); Ci=np.linalg.inv(C2); ker=np.exp(-0.5*np.einsum('ij,ij->j',q,Ci@q)).reshape(n,n); ker/=ker.sum()
        omega_ratio=(tgt['bmaj']*tgt['bmin'])/(bmaj*bmin)
        # per-channel convolution, then windows
        small=np.zeros_like(cube)
        for i in range(cube.shape[0]): small[i]=fftconvolve(np.nan_to_num(cube[i]),ker,mode='same')*omega_ratio
        Ks=disk(0.5*np.sqrt(0.93*0.74)/ps)
        Ms={}
        for k,(a,b) in wins.items():
            m,nn=winmap(small,v,a,b,dv)
            if nn>=0.8*(b-a)/dv: Ms[k]=ndimage.convolve(m,Ks,mode='constant')
        cen_=cen
        P(f'  T4 Xb4b convolved to 0.93x0.74 (beam-area factor {omega_ratio:.1f}); Sblue-peak window values (Jy/beam km/s):')
        for k in Ms: 
            sel=cen&(rad('S')>1.5)&(rad('N')>1.5); xs=Ms[k][sel]
            P(f'     {k}: at Sblue peak {Ms[k][iy,ix]:+.4f}; at S {Ms[k][int(round(pos["S"][1])),int(round(pos["S"][0]))]:+.4f}; field MAD-sigma {1.4826*np.median(abs(xs-np.median(xs))):.4f}')
        out[f'Xb4b_{line}_smoothed']={k:float(Ms[k][iy,ix]) for k in Ms}
        np.savez(f'/workspace/out/s11_Xb4b_{line}_smoothed_maps.npz',**{k:Ms[k] for k in Ms},x0=x0,y0=y0)
    np.savez(f'/workspace/out/s11_{ous}_{line}_maps.npz',**M,x0=x0,y0=y0,ix=ix,iy=iy)
open('/workspace/out/s11_wing_tests.txt','w').write('\n'.join(log)); json.dump(out,open('/workspace/out/s11_wing_tests.json','w'),indent=1)
