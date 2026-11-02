# Q8: cross-array line flux consistency. Xb4b (0.25") and Xb4d (0.94x0.74") HCO+/HCN 3-2 cubes are convolved (analytic Fourier Gaussian,
# Jy/beam rescaled by beam-area ratio) to the Xb4f beam (header beam), per-pixel line-free baseline removed, and integrated over |v-VSYS|<350 km/s
# (core) and the full-window |v|<600 (all) in circular apertures (centred on the N-S midpoint) and at region peaks. Output: out/s26_recovery.json/.txt, s26_spectra_*.npz
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
import json
LINES=(('HCO+',267.55763e9,{'Xb4f':'spw16','Xb4d':'spw25','Xb4b':'spw25'}),('HCN',265.88618e9,{'Xb4f':'spw18','Xb4d':'spw27','Xb4b':'spw27'}))
HALF=16.0   # crop half-width ["]
CEN=[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)]
def cov(a,b,pa):
    th=np.deg2rad(pa); e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
def transfer(shape,ps,b0,tgt):
    # Fourier transfer function of the Gaussian that takes beam b0 to tgt (arcsec, arcsec, PA deg E of N); x is RA (east = -x), y is Dec
    C2=cov(*tgt)-cov(*b0); w_,U=np.linalg.eigh(C2); C2=U@np.diag(np.clip(w_,1e-8,None))@U.T
    fy=np.fft.fftfreq(shape[0],ps); fx=np.fft.fftfreq(shape[1],ps); FX,FY=np.meshgrid(fx,fy)
    # covariance in (east, north); image x increases to west -> east comp = -x
    q=np.stack([-FX.ravel(),FY.ravel()])
    return np.exp(-2*np.pi**2*np.einsum('ij,ij->j',q,C2@q)).reshape(shape)
def load(ous,spw,f0,pad):
    fn=files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]; h_,hd,w,nu=cube_info(fn)
    v=vrad(nu,f0)-VSYS; ps=abs(hd['CDELT1'])*3600; wc=WCS(hd).celestial
    cx,cy=wc.world_to_pixel_values(*CEN); r=int(HALF/ps); x0,y0=int(round(cx-r)),int(round(cy-r))
    n=2*r; ok=np.where((v>-1400)&(v<600))[0]; i0,i1=ok[0],ok[-1]+1
    cube=np.array(h_[0].data[0,i0:i1,y0:y0+n,x0:x0+n],dtype=np.float32); cube=np.nan_to_num(cube)
    beam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
    return cube,v[i0:i1],ps,wc,x0,y0,beam
def conv(cube,ps,b0,tgt):
    if b0[0]*b0[1]>0.5*tgt[0]*tgt[1]*0+1e9: return cube
    n=cube.shape[1]; P=2*n; T=transfer((P,P),ps,b0,tgt); fac=tgt[0]*tgt[1]/(b0[0]*b0[1]); out=np.empty_like(cube)
    for i in range(cube.shape[0]):
        a=np.zeros((P,P),np.float32); a[:n,:n]=cube[i]; out[i]=np.fft.irfft2(np.fft.rfft2(a)*T[:, :P//2+1],s=(P,P))[:n,:n]*fac
    return out
res={}; txt=[]
for line,f0,spws in LINES:
    cf,vf,psf,wcf,xf,yf,tgt=load('Xb4f',spws['Xb4f'],f0,0)
    D={'Xb4f':(cf,vf,psf,wcf,xf,yf)}
    for ous in ('Xb4d','Xb4b'):
        c,v,ps,wc,x0,y0,b0=load(ous,spws[ous],f0,0)
        # block-average Xb4b to ~0.19" before convolving to save time (beam 0.25" -> pixel sampling still fine for 5" target)
        if ous=='Xb4b':
            k=4; n=c.shape[1]//k*k; c=c[:,:n,:n].reshape(c.shape[0],n//k,k,n//k,k).mean(axis=(2,4)); ps*=k
            # pixel-average preserves Jy/beam; x0,y0 offsets handled via a shifted WCS below
            wc=WCS(wc.to_header()); 
        # recentre: use sky position of crop centre rather than index
        cc=conv(c,ps,b0,tgt)
        D[ous]=(cc,v,ps,wc,x0,y0)
        print(line,ous,'beam',b0,'->',tgt,'cube',cc.shape,flush=True)
    # per-OUS: baseline-subtract, integrate
    ap={}; spec={}
    for ous,(c,v,ps,wc,x0,y0) in D.items():
        ref=(v<-600)|((v>450)&(v<550)); core=abs(v)<350; allv=(v>-600)&(v<550); dv=abs(np.median(np.diff(v)))
        base=np.median(c[ref],axis=0); c=c-base
        # noise per channel in aperture-integrated spectrum later
        n=c.shape[1]
        # pixel coords of region centres in the crop frame (Xb4b binned: pixel size ps, origin at crop corner)
        wcs_full=wc
        Y,X=np.mgrid[:n,:n]
        def pxy(ra,dec):
            if ous=='Xb4b': # binned: index (full-crop) -> centre-of-block convention
                px,py=wc.world_to_pixel_values(ra,dec); return (px-x0-1.5)/4,(py-y0-1.5)/4
            px,py=wc.world_to_pixel_values(ra,dec); return px-x0,py-y0
        beamarea=1.1331*tgt[0]*tgt[1]/ps**2
        m0c=c[core].sum(0)*dv; m0a=c[allv].sum(0)*dv
        d={}
        for R in (4,8,12):
            px,py=pxy(*CEN); m=np.hypot(X-px,Y-py)*ps<=R
            d[f'ap{R}_core']=float(m0c[m].sum()/beamarea); d[f'ap{R}_all']=float(m0a[m].sum()/beamarea)
            spec[f'{ous}_ap{R}']=c[:,m].sum(1)/beamarea
            # noise: scatter of aperture sum in line-free channels (blue baseline)
            lf=c[(v<-700)][:,m].sum(1)/beamarea; d[f'ap{R}_chan_rms_Jy']=float(np.std(lf)); d[f'ap{R}_core_err']=float(np.std(lf)*dv*np.sqrt(core.sum()))
        for reg in ('N','S','Bridge','E1','E2','W','Sblue'):
            px,py=pxy(REGIONS[reg][0],REGIONS[reg][1]); m=np.hypot(X-px,Y-py)*ps<=1.0
            idx=np.unravel_index(np.argmax(np.where(m,m0c,-1e9)),m0c.shape)
            d[f'pk_{reg}_core']=float(m0c[idx]); d[f'pk_{reg}_all']=float(m0a[idx]); d[f'pk_{reg}_peakchan_mJy']=float(c[:,idx[0],idx[1]].max()*1e3)
            spec[f'{ous}_{reg}']=c[:,idx[0],idx[1]]
        ap[ous]=d
    np.savez(f'/workspace/out/s26_spectra_{line}.npz',v=D['Xb4f'][1],**{k:v for k,v in spec.items()},vb=D['Xb4b'][1],vd=D['Xb4d'][1])
    res[line]=ap
    txt.append(f'## {line} 3-2 at the Xb4f beam {tgt[0]:.2f}x{tgt[1]:.2f}" PA{tgt[2]:.0f}; Jy km/s (aperture sum / beam area), core |v|<350, all -600..+550')
    for key in [k for k in ap['Xb4f'] if not k.endswith('mJy') and 'rms' not in k and 'err' not in k]:
        a=ap['Xb4f'][key]; txt.append(f'  {key:18s} Xb4f {a:8.3f}  Xb4d {ap["Xb4d"][key]:8.3f} ({ap["Xb4d"][key]/a:5.2f})  Xb4b {ap["Xb4b"][key]:8.3f} ({ap["Xb4b"][key]/a:5.2f})')
    for R in (4,8,12): txt.append(f'  ap{R} core 1-sigma err: Xb4f {ap["Xb4f"][f"ap{R}_core_err"]:.3f} Xb4d {ap["Xb4d"][f"ap{R}_core_err"]:.3f} Xb4b {ap["Xb4b"][f"ap{R}_core_err"]:.3f}')
json.dump(res,open('/workspace/out/s26_recovery.json','w'),indent=1)
open('/workspace/out/s26_recovery.txt','w').write('\n'.join(txt)); print('\n'.join(txt))
