# Q12 step a: native-resolution moment-0 maps (|v-VSYS|<350 km/s, per-pixel median baseline from line-free channels v<-600 | 450<v<600 km/s, pbcor, Jy/beam km/s)
# for HCO+ and HCN 3-2 in Xb4b (binned x4 to 0.19"), Xb4d and Xb4f, on a +-HALF" crop centred on the N-S midpoint. Also NREAL sign-flipped line-free
# pseudo-moment maps (noise realisations with the same spatial correlation, scaled to the core channel count). Output: out/q12/maps_<ous>_<line>.npz
# usage: python3 -I s56a_maps.py <ous> <line> [half-width arcsec, default 19]   (run for all 6 combinations by run_q12.sh)
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
HALF=float(sys.argv[3]) if len(sys.argv)>3 else 19.0; NREAL=16; SUF='' if HALF==19.0 else f'_h{HALF:g}'
CEN=[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)]
LINES={'HCO+':(267.55763e9,{'Xb4f':'spw16','Xb4d':'spw25','Xb4b':'spw25'}),'HCN':(265.88618e9,{'Xb4f':'spw18','Xb4d':'spw27','Xb4b':'spw27'})}
def run(ous,line):
    f0,spws=LINES[line]; fn=files(f'*{ous}.ngc_3256_sci.{spws[ous]}.cube.I.pbcor.fits')[0]; h_,hd,w,nu=cube_info(fn)
    v=vrad(nu,f0)-VSYS; ps=abs(hd['CDELT1'])*3600; wc=WCS(hd).celestial
    cx,cy=wc.world_to_pixel_values(*CEN); r=int(HALF/ps); x0,y0=int(round(cx-r)),int(round(cy-r)); n=2*r
    ok=np.where((v>-1100)&(v<700))[0]; i0,i1=ok[0],ok[-1]+1
    cube=np.nan_to_num(np.array(h_[0].data[0,i0:i1,y0:y0+n,x0:x0+n],dtype=np.float32)); v=v[i0:i1]
    k=4 if ous=='Xb4b' else 1
    if k>1:
        n2=n//k*k; cube=cube[:,:n2,:n2].reshape(cube.shape[0],n2//k,k,n2//k,k).mean(axis=(2,4))
    ref=(v<-600)|((v>450)&(v<600)); core=abs(v)<350; dv=abs(np.median(np.diff(v)))
    c=cube-np.median(cube[ref],axis=0)
    m0=c[core].sum(0)*dv
    rng=np.random.default_rng(1); lf=c[ref]; sc=np.sqrt(core.sum()/ref.sum())
    noise=np.stack([(lf*rng.choice([-1,1],lf.shape[0])[:,None,None]).sum(0)*dv*sc for _ in range(NREAL)])
    # pixel centre of block j (binned) in full-crop pixel units: (j*k+(k-1)/2); sky->index: (px-x0-(k-1)/2)/k
    hdr=dict(ps=ps*k,x0=x0,y0=y0,k=k,beam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),dv=dv,nlf=int(ref.sum()),ncore=int(core.sum()),
             crpix=(hd['CRPIX1'],hd['CRPIX2']))
    # sky positions of every binned pixel centre via WCS for later region lookup
    jj=np.arange(m0.shape[0]); X=x0+jj*k+(k-1)/2; ra,dec=wc.pixel_to_world_values(X,np.full_like(X,y0+(m0.shape[0]//2)*k+(k-1)/2)); 
    np.savez(f'/workspace/out/q12/maps_{ous}_{line}{SUF}.npz',m0=m0,noise=noise,**{kk:np.array(vv) for kk,vv in hdr.items()})
    print(ous,line,'beam',hdr['beam'],'m0 peak',m0.max(),'rms_lf(median pix)',np.median(noise.std(0)),flush=True)
if __name__=='__main__': run(sys.argv[1],sys.argv[2])
