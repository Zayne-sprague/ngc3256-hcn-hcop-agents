# Cross-OUS flux consistency: Xb4b (0.25") convolved to the Xb4d (0.94x0.74") beam vs Xb4d, same cube definition (HCO+ & HCN 3-2 core moment-0, continuum) at S, N, E1, E2
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
from scipy.signal import fftconvolve
import json
def cov(a,b,pa):
    th=np.deg2rad(pa); e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
def kernel(ps,b0,tgt):
    C2=cov(*tgt)-cov(*b0); w_,U=np.linalg.eigh(C2); C2=U@np.diag(np.clip(w_,1e-6,None))@U.T
    n=int(8*np.sqrt(w_.max())/ps)|1; ky,kx=np.mgrid[:n,:n]-n//2; q=np.stack([(-kx*ps).ravel(),(ky*ps).ravel()])
    k=np.exp(-0.5*np.einsum('ij,ij->j',q,np.linalg.inv(C2)@q)).reshape(n,n); return k/k.sum()
tgt=(0.93,0.74,-86.0); out={}
def val(d,wc,x0,y0,ps,reg,r=0.5):
    px,py=wc.world_to_pixel_values(REGIONS[reg][0],REGIONS[reg][1]); px-=x0; py-=y0
    Y,X=np.mgrid[:d.shape[0],:d.shape[1]]; m=np.hypot(X-px,Y-py)*ps<=r; return float(np.max(np.where(m,d,-1e9)))
# continuum
rows=[]
dd={}
for ous,pat in (('Xb4b','*Xb4b*25_27_29_31.cont.I.pbcor.fits'),('Xb4d','*Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')):
    d,h=load2d(files(pat)[0]); dd[ous]=(d,h)
d,h=dd['Xb4b']; ps=abs(h['CDELT1'])*3600; k=kernel(ps,(h['BMAJ']*3600,h['BMIN']*3600,h['BPA']),tgt)
c=fftconvolve(np.nan_to_num(d),k,mode='same')*(tgt[0]*tgt[1])/(h['BMAJ']*3600*h['BMIN']*3600)
wcb=WCS(h).celestial; wcd=WCS(dd['Xb4d'][1]).celestial; psd=0.15
print('continuum peak (mJy/beam): region  Xb4b-conv  Xb4d  ratio')
for reg in ('N','S','E1','E2'):
    a=val(c,wcb,0,0,ps,reg)*1e3; b=val(dd['Xb4d'][0],wcd,0,0,psd,reg)*1e3; print(f'  {reg} {a:.3f} {b:.3f} {a/b:.3f}'); out[f'cont_{reg}']=a/b
# lines (core |v-VSYS|<350 moment-0 with per-pixel line-free baseline), crop +-12"
for line,spw,f0 in (('HCO+','spw25',267.55763),('HCN','spw27',265.88618)):
    M={}
    for ous in ('Xb4b','Xb4d'):
        fn=files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]; h_,hd,w,nu=cube_info(fn); v=vrad(nu,f0*1e9)-VSYS; ps_=abs(hd['CDELT1'])*3600
        wc=WCS(hd).celestial; cx,cy=wc.world_to_pixel_values(*[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)]); r=int(12/ps_)
        x0,y0=int(cx-r),int(cy-r); ok=np.where(abs(v)<950)[0]; i0,i1=ok[0],ok[-1]+1
        cube=np.array(h_[0].data[0,i0:i1,y0:y0+2*r,x0:x0+2*r],dtype=np.float32); v=v[i0:i1]; dv=abs(np.median(np.diff(v)))
        cube-=np.nanmedian(cube[abs(v)>500],axis=0)
        if ous=='Xb4b':
            b0=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); k=kernel(ps_,b0,tgt); fac=(tgt[0]*tgt[1])/(b0[0]*b0[1])
            cube=np.array([fftconvolve(np.nan_to_num(cube[i]),k,mode='same')*fac for i in range(cube.shape[0])])
        m0=np.nansum(cube[abs(v)<350],axis=0)*dv; M[ous]=(m0,wc,x0,y0,ps_)
    print(f'{line} 3-2 core moment-0 (Jy/beam km/s, Xb4d beam): region Xb4b-conv Xb4d ratio')
    for reg in ('N','S','Bridge','E1','E2'):
        a=val(*M['Xb4b'][:1],M['Xb4b'][1],M['Xb4b'][2],M['Xb4b'][3],M['Xb4b'][4],reg); b=val(M['Xb4d'][0],M['Xb4d'][1],M['Xb4d'][2],M['Xb4d'][3],M['Xb4d'][4],reg)
        print(f'  {reg} {a:.3f} {b:.3f} {a/b:.3f}'); out[f'{line}_{reg}']=a/b
json.dump(out,open('/workspace/out/s15_crossous.json','w'),indent=1)
