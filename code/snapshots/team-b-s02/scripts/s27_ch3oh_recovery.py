# Q8/Q9: does smoothing Xb4d to the Xb4f beam recover the Xb4f CH3OH 254 GHz flux (and HCO+ 3-2 in the same windows)?
# Moment-0 (CH3OH frame v0=2783 +-150 km/s, per-pixel linear baseline, as s21) of Xb4d, convolved to the Xb4f beam (Fourier Gaussian, Jy/beam rescaled),
# vs Xb4f native; Xb4f native beam of spw20 (CH3OH) / spw16 (HCO+).
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.wcs import WCS
from regions import REGIONS
import importlib.util
exec(open('/workspace/scripts/s21_ch3oh254_maps.py').read().split('def beam(H)')[0])   # reuse load(), mom0(), constants
def cov(a,b,pa):
    th=np.deg2rad(pa); e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
def smooth(m,ps,b0,tgt):
    P=2*max(m.shape); C2=cov(*tgt)-cov(*b0); w_,U=np.linalg.eigh(C2); C2=U@np.diag(np.clip(w_,1e-8,None))@U.T
    FX,FY=np.meshgrid(np.fft.fftfreq(P,ps),np.fft.fftfreq(P,ps)); q=np.stack([-FX.ravel(),FY.ravel()])
    T=np.exp(-2*np.pi**2*np.einsum('ij,ij->j',q,C2@q)).reshape(P,P)
    a=np.zeros((P,P)); a[:m.shape[0],:m.shape[1]]=np.nan_to_num(m)
    return np.fft.ifft2(np.fft.fft2(a)*T).real[:m.shape[0],:m.shape[1]]*tgt[0]*tgt[1]/(b0[0]*b0[1])
out={}
M={}
for name,(spC,spH) in {'CH3OH':('29','25'),'HCO+':('29','25')}.items(): pass
for ous,spf in (('Xb4d',{'CH3OH':'29','HCO+':'25'}),('Xb4f',{'CH3OH':'20','HCO+':'16'})):
    for tr in ('CH3OH','HCO+'):
        d,H,nu=load(ous,spf[tr])
        rest,v0=(REST,V0) if tr=='CH3OH' else (HCO,V0-9)
        m,sig,dv=mom0(d,H,nu,rest,v0,HW); M[(ous,tr)]=(m,H,WCS(H).celestial)
tgt={tr:(M[('Xb4f',tr)][1]['BMAJ']*3600,M[('Xb4f',tr)][1]['BMIN']*3600,M[('Xb4f',tr)][1]['BPA']) for tr in ('CH3OH','HCO+')}
for tr in ('CH3OH','HCO+'):
    mD,HD,wD=M[('Xb4d',tr)]; ps=abs(HD['CDELT1'])*3600; b0=(HD['BMAJ']*3600,HD['BMIN']*3600,HD['BPA'])
    sm=smooth(mD,ps,b0,tgt[tr]); mF,HF,wF=M[('Xb4f',tr)]
    # noise: MAD of moment-0 over an annulus 7-14" from centre (Xb4f) / smoothed Xb4d
    def mad(m,w,psx,r0,r1):
        ny,nx=m.shape; yy,xx=np.mgrid[:ny,:nx]; cx,cy=w.all_world2pix(156.96346,-43.90461,0); r=np.hypot(xx-cx,yy-cy)*psx; s=(r>r0)&(r<r1)&(m!=0)&np.isfinite(m)
        return 1.4826*np.nanmedian(np.abs(m[s]-np.nanmedian(m[s])))
    sD=mad(sm,wD,ps,10,15); sF=mad(mF,wF,0.7,10,15)
    d_={'beam_target':tgt[tr],'mad_smoothedXb4d':sD,'mad_Xb4f':sF}
    for k in ('N','S','Bridge','E1','E2','W','Sblue'):
        res=[]
        for m,w,p in ((sm,wD,ps),(mF,wF,0.7)):
            x,y=w.all_world2pix(REGIONS[k][0],REGIONS[k][1],0); x=int(round(float(x))); y=int(round(float(y))); r=int(round(1.0/p))
            sub=m[max(y-r,0):y+r+1,max(x-r,0):x+r+1]; res.append(float(np.nanmax(sub)) if r>0 else float(m[y,x]))
        d_[k]={'Xb4d_smoothed':res[0],'Xb4f':res[1],'ratio':res[0]/res[1] if res[1]!=0 else None}
    out[tr]=d_
    print(tr,'target beam',[round(x,2) for x in tgt[tr]],f'MAD smoothed Xb4d {sD:.3f}, Xb4f {sF:.3f}')
    for k in ('N','S','Bridge','E1','E2','W','Sblue'): print(f"  {k:7s} Xb4d_sm {d_[k]['Xb4d_smoothed']:7.3f}  Xb4f {d_[k]['Xb4f']:7.3f}  ratio {d_[k]['ratio']:.2f}")
json.dump(out,open('/workspace/out/s27_ch3oh_recovery.json','w'),indent=1)
