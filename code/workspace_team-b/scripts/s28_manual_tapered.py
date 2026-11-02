# Q8: hand-made QA2 products vs standard cubes. Moment-0 (raw channel sum x channel width, over the channel-edge range of the coarse cube; no baseline fit
# because the coarse cubes are short) at region peaks (max within 0.5" of region), after convolving to a common beam:
#  - Xb4d tapered (0.96x0.83") vs Xb4d standard 3.9-MHz cube convolved to the tapered beam (HCO+ spw25, HCN spw27)
#  - Xb4b manual (0.22x0.20") convolved to the standard Xb4b beam (0.25x0.23) vs Xb4b standard (HCO+ spw25); manual non-pbcor vs pbcor.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
from scipy.ndimage import gaussian_filter
import json
exec(open('/workspace/scripts/s27_ch3oh_recovery.py').read().split('out={}')[0].split('exec(open')[0])  # imports
def cov(a,b,pa):
    th=np.deg2rad(pa); e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
def smooth(m,ps,b0,tgt):
    n=m.shape[-1]; P=2*n; C2=cov(*tgt)-cov(*b0); w_,U=np.linalg.eigh(C2); C2=U@np.diag(np.clip(w_,1e-8,None))@U.T
    FX,FY=np.meshgrid(np.fft.fftfreq(P,ps),np.fft.fftfreq(P,ps)); q=np.stack([-FX.ravel(),FY.ravel()])
    T=np.exp(-2*np.pi**2*np.einsum('ij,ij->j',q,C2@q)).reshape(P,P); a=np.zeros((P,P)); a[:n,:n]=np.nan_to_num(m)
    return np.fft.ifft2(np.fft.fft2(a)*T).real[:n,:n]*tgt[0]*tgt[1]/(b0[0]*b0[1])
def crop_m0(fn,vlo,vhi,f0,half,beamfix=None):
    h_,hd,w,nu=cube_info(fn); v=vrad(nu,f0)-VSYS; ps=abs(hd['CDELT1'])*3600; dv=abs(np.median(np.diff(v)))
    sel=np.where((v>vlo)&(v<vhi))[0]; wc=WCS(hd).celestial; cx,cy=wc.world_to_pixel_values(*[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)])
    r=int(half/ps); x0,y0=int(round(cx-r)),int(round(cy-r))
    c=np.nan_to_num(np.array(h_[0].data[0,sel[0]:sel[-1]+1,y0:y0+2*r,x0:x0+2*r],dtype=np.float32))
    return c.sum(0)*dv,(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),ps,wc,x0,y0,v,dv
def peak(m,wc,x0,y0,ps,reg,rad=0.5):
    px,py=wc.world_to_pixel_values(REGIONS[reg][0],REGIONS[reg][1]); px-=x0; py-=y0; Y,X=np.mgrid[:m.shape[0],:m.shape[1]]
    return float(np.max(np.where(np.hypot(X-px,Y-py)*ps<=rad,m,-1e9)))
RG=('N','S','Bridge','E1','E2','Sblue'); out={}; lines=[]
for line,f0,spw in (('HCO+',267.55763e9,'spw25'),('HCN',265.88618e9,'spw27')):
    ft=files(f'*Xb4d_ngc_3256_sci.{spw}.cube.I.tapered.image.pbcor.fits')[0]
    ht=cube_info(ft); vt=vrad(ht[3],f0)-VSYS; dvt=abs(np.median(np.diff(vt))); lo,hi=vt.min()-dvt/2,vt.max()+dvt/2
    mt,bt,pst,wct,xt,yt,_,_=crop_m0(ft,lo-1,hi+1,f0,12.0)
    fs=files(f'*Xb4d.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]
    ms,bs,pss,wcs_,xs,ys,vs,dvs=crop_m0(fs,lo,hi,f0,12.0)
    # edge-weighted: standard channels centred within [lo,hi]; width mismatch (<=1 channel=3.9 km/s over range %.0f km/s) is negligible
    # (tapered cube pixel 0.2", standard 0.15"): evaluate on native grids after smoothing standard to the tapered beam
    ms=smooth(ms,pss,bs,bt)
    out[f'Xb4d_tapered_{line}']={'range_kms':[float(lo),float(hi)],'beam_tapered':bt,'beam_std':bs}
    lines.append(f'## Xb4d tapered vs standard (conv. to tapered beam {bt[0]:.2f}x{bt[1]:.2f}), {line} 3-2, v in [{lo:.0f},{hi:.0f}] km/s, Jy/beam km/s')
    for r in RG:
        a=peak(mt,wct,xt,yt,pst,r); b=peak(ms,wcs_,xs,ys,pss,r); out[f'Xb4d_tapered_{line}'][r]=[a,b,a/b]; lines.append(f'  {r:7s} tapered {a:8.3f}  standard {b:8.3f}  ratio {a/b:5.3f}')
# Xb4b manual
f0=267.55763e9
fm=files('*Xb4b*spw25.cube.I.manual.image.pbcor.fits')[0]; fmn=files('*Xb4b*spw25.cube.I.manual.image.fits')[0]
hm=cube_info(fm); vm=vrad(hm[3],f0)-VSYS; dvm=abs(np.median(np.diff(vm))); lo,hi=-540.,320.
mm,bm,psm,wcm,xm,ym,_,_=crop_m0(fm,lo,hi,f0,6.0); mmn,_,_,_,_,_,_,_=crop_m0(fmn,lo,hi,f0,6.0)
fb=files('*Xb4b.ngc_3256_sci.spw25.cube.I.pbcor.fits')[0]; mb,bb,psb,wcb,xb,yb,_,_=crop_m0(fb,lo,hi,f0,6.0)
mmc=smooth(mm,psm,bm,bb); mmnc=smooth(mmn,psm,bm,bb)
# pixel scales differ (0.038 vs 0.047"): compare region peaks only
lines.append(f'## Xb4b manual (smoothed from {bm[0]:.2f}x{bm[1]:.2f} to standard beam {bb[0]:.3f}x{bb[1]:.3f}) vs standard, HCO+ 3-2, v in [{lo:.0f},{hi:.0f}], Jy/beam km/s; also non-pbcor manual (flat-noise) at same places')
out['Xb4b_manual_HCO+']={}
for r in ('N','S','Bridge','E1','Sblue'):
    a=peak(mmc,wcm,xm,ym,psm,r); a2=peak(mmnc,wcm,xm,ym,psm,r); b=peak(mb,wcb,xb,yb,psb,r); out['Xb4b_manual_HCO+'][r]=[a,a2,b,a/b]; lines.append(f'  {r:7s} manual_pbcor {a:8.3f} manual_nonpb {a2:8.3f} standard {b:8.3f}  ratio {a/b:5.3f}')
json.dump(out,open('/workspace/out/s28_manual_tapered.json','w'),indent=1); open('/workspace/out/s28_manual_tapered.txt','w').write('\n'.join(lines)); print('\n'.join(lines))
