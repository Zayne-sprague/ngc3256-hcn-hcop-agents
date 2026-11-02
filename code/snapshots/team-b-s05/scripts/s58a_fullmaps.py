# Q12 round 2, step a: FULL-FIELD native-beam moment-0 maps (|v-VSYS|<350 km/s, per-pixel median baseline from line-free channels v<-600 | 450<v<600 km/s, pbcor,
# Jy/beam km/s) for the three ACA sets Xb55 (J=2-1), Xb4f (3-2), Xb47 (4-3), HCO+ and HCN; 16 sign-flipped line-free noise realisations; fitted Gaussian PB (from the
# noise-map ~1/PB) so that maps can be un-corrected (map*PB) for the zero-spacing / bowl analysis. Also the line-free continuum (cont.I.pbcor) maps.
# usage: python3 -I s58a_fullmaps.py   -> out/q12b/full_<ous>_<line>.npz, out/q12b/cont_<ous>.npz
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy.optimize import least_squares
NREAL=16
JOBS={'Xb55':{'HCN':('spw16',177.26111e9),'HCO+':('spw18',178.37506e9)},'Xb4f':{'HCO+':('spw16',267.55763e9),'HCN':('spw18',265.88618e9)},
      'Xb47':{'HCN':('spw16',354.50547e9),'HCO+':('spw18',356.73424e9)}}
def fit_pb(e,ps):
    n=e.shape[0]; Y,X=np.mgrid[:n,:n]; ok=np.isfinite(e)&(e>0)
    f=lambda p:(p[0]+p[3]*((X-p[1])**2+(Y-p[2])**2)*ps**2)[ok]-np.log(e[ok])
    s=least_squares(f,[np.log(np.median(e[ok])),n/2,n/2,1e-3]); a,x0,y0,b=s.x
    return np.sqrt(4*np.log(2)/b),x0,y0
def run(ous,line):
    spw,f0=JOBS[ous][line]; fn=files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]; h_,hd,w,nu=cube_info(fn)
    v=vrad(nu,f0)-VSYS; ps=abs(hd['CDELT1'])*3600
    ok=np.where((v>-1100)&(v<700))[0]; i0,i1=ok[0],ok[-1]+1
    cube=np.nan_to_num(np.array(h_[0].data[0,i0:i1],dtype=np.float32)); v=v[i0:i1]
    ref=(v<-600)|((v>450)&(v<600)); core=abs(v)<350; dv=abs(np.median(np.diff(v)))
    c=cube-np.median(cube[ref],axis=0); m0=c[core].sum(0)*dv
    rng=np.random.default_rng(1); lf=c[ref]; sc=np.sqrt(core.sum()/ref.sum())
    noise=np.stack([(lf*rng.choice([-1,1],lf.shape[0])[:,None,None]).sum(0)*dv*sc for _ in range(NREAL)])
    fw,px,py=fit_pb(noise.std(0),ps)
    wc=WCS(hd).celestial; ra,dec=wc.pixel_to_world_values(px,py)
    np.savez(f'/workspace/out/q12b/full_{ous}_{line}.npz',m0=m0,noise=noise,ps=ps,beam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),pbfwhm=fw,pbx=px,pby=py,pbra=ra,pbdec=dec,nref=int(ref.sum()),ncore=int(core.sum()))
    print(ous,line,'shape',m0.shape,'ps',ps,'beam',hd['BMAJ']*3600,hd['BMIN']*3600,'PB fwhm %.1f"'%fw,'centre',ra,dec,'nref',ref.sum(),'m0 max',m0.max(),flush=True)
def cont(ous):
    fn=files(f'*{ous}.ngc_3256_sci.spw*_*.cont.I.pbcor.fits')[0]; d,hd=load2d(fn); ps=abs(hd['CDELT1'])*3600
    np.savez(f'/workspace/out/q12b/cont_{ous}.npz',img=np.nan_to_num(d),ps=ps,beam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),nu=hd.get('CRVAL3',0)); print('cont',ous,d.shape,ps)
if __name__=='__main__':
    from multiprocessing import Pool
    jobs=[(o,l) for o in JOBS for l in JOBS[o]]
    with Pool(6) as p: p.starmap(run,jobs)
    for o in ('Xb55','Xb4f','Xb47','Xb4d','Xb4b'):
        try: cont(o)
        except Exception as ex: print('cont fail',o,ex)
