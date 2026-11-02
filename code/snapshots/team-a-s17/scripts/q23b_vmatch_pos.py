# Q23 round 2: (1) CO T_b velocity-matched to the HCN/HCO+ lines, (2) position error budget (offset grid + Monte Carlo),
# using the 0.18" cube convolved to circular beams with proper elliptical kernels (radio_beam). Output out/q23b/vmatch_pos.json
import numpy as np, json, warnings, pandas as pd; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve_fft
from astropy import units as u
from radio_beam import Beam
from scipy.ndimage import map_coordinates
C=299792.458; H=6.62607e-27; KB=1.380649e-16; CL=2.99792458e10
P='/archive/files/2015.1.00714.S/member.uid___A001_X2d6_Xae.ari_l.NGC_3256_sci.spw1_229288MHz.12m.cube.I.pbcor.fits'
RA0,DE0=156.96346,-43.90461
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
h=fits.open(P,memmap=True); hd=h[0].header; cube=h[0].data[0]
nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/230.538e9)*C
nat=Beam(hd['BMAJ']*u.deg,hd['BMIN']*u.deg,hd['BPA']*u.deg); pix=abs(hd['CDELT2'])*3600
w=WCS(hd).celestial; HW=int(round(6/pix))
def TRJ(S,nu0,a,b): return 1.222e3*S*1e3/((nu0/1e9)**2*a*b)
def TPl(S,nu0,a,b):
    I=S/(1.1331*a*b*(np.pi/180/3600)**2)*1e-23
    return H*nu0/KB/np.log1p(2*H*nu0**3/CL**2/I)
wd=pd.read_csv('/workspace/out/q5/q5e_width.csv')
def dense(line,src,ap=0.0):
    r=wd[(wd.line==line)&(wd.src==src)&(np.isclose(wd.aperture_r,ap))].iloc[0]; return float(r.v0),float(r.fwhm)
ks=np.where((v>2250)&(v<2845))[0]; koff=np.where((v>1700)&(v<2250))[0]
chan=np.concatenate([koff,ks]); isl=np.r_[np.zeros(len(koff),bool),np.ones(len(ks),bool)]
vs=v[chan]; nus=nu[chan]
BEAMS=[0.26,0.38,0.57,1.0,1.7]; SIGPOS=[0.05,0.07,0.10]
rng=np.random.default_rng(1)
def sample(S,x,y):  # bilinear at fractional pix (x,y) arrays -> (n,nch)
    x0=np.floor(x).astype(int); y0=np.floor(y).astype(int); fx=x-x0; fy=y-y0
    return (S[:,y0,x0]*(1-fx)*(1-fy)+S[:,y0,x0+1]*fx*(1-fy)+S[:,y0+1,x0]*(1-fx)*fy+S[:,y0+1,x0+1]*fx*fy).T
res={'pix':pix,'native_beam':[nat.major.to(u.arcsec).value,nat.minor.to(u.arcsec).value],'cores':{}}
for nm,(x0,y0) in CORES.items():
    ra=RA0+x0/3600/np.cos(np.deg2rad(DE0)); de=DE0+y0/3600
    fx,fy=[float(a) for a in w.world_to_pixel_values(ra,de)]; ix,iy=int(round(fx)),int(round(fy))
    sub=np.array(cube[chan][:,iy-HW:iy+HW+1,ix-HW:ix+HW+1],float); sub=np.nan_to_num(sub)
    cx,cy=HW+(fx-ix),HW+(fy-iy)   # fractional centre in subcube
    lines={}
    for ln in ['HCN','HCO+']:
        v0,fw=dense(ln,nm); lines[ln]=(v0,fw)
    cr={'dense':{k:list(x) for k,x in lines.items()},'beams':{}}
    for b in BEAMS:
        tgt=Beam(b*u.arcsec,b*u.arcsec,0*u.deg)
        ker=tgt.deconvolve(nat).as_kernel(pix*u.arcsec)
        fac=b*b/(nat.major.to(u.arcsec).value*nat.minor.to(u.arcsec).value)
        S=np.array([convolve_fft(s,ker,boundary='fill',fill_value=0,allow_huge=True,normalize_kernel=True) for s in sub])*fac
        rms=np.nanstd(S[~isl][:,HW,HW])
        T=TRJ(S,nus[:,None,None],b,b)  # RJ K, shape (nch,ny,nx)
        e_T=TRJ(rms,230e9/1e0*0+nus[isl].mean(),b,b)
        dv=abs(np.median(np.diff(v)))
        def stats(sp):  # sp: RJ K spectrum on chan grid
            s=sp[isl]; vv=vs[isl]; k=np.argmax(s); out={'peak_RJ':float(s[k]),'v_peak':float(vv[k])}
            for ln,(v0,fw) in lines.items():
                m=np.abs(vv-v0)<=fw/2; cov=m.sum()*dv/fw
                g=np.exp(-0.5*((vv-v0)/(fw/2.3548))**2); gw=float((s*g).sum()/g.sum()); gcov=float(g.sum()/(fw/2.3548*np.sqrt(2*np.pi)/dv))
                out[ln]={'T_at_v0':float(np.mean(s[np.abs(vv-v0)<=dv*1.01])),'T_mean_FWHM':float(s[m].mean()) if m.any() else np.nan,'cov_FWHM':float(cov),
                         'T_gauss_wt':gw,'cov_gauss':gcov,'T_max_in_FWHM':float(s[m].max()) if m.any() else np.nan}
            return out
        sp0=T[:,int(round(cy)),int(round(cx))]
        rec={'noise_RJ_K':float(e_T),'at_core_pixel':stats(T[:,int(round(cy)),int(round(cx))]),'at_core_interp':stats(sample(T,np.array([cx]),np.array([cy]))[0])}
        # offset grid: 8 azimuths
        og={}
        for off in [0,0.03,0.06,0.10,0.15]:
            ang=np.arange(8)*np.pi/4 if off>0 else np.array([0.])
            xs=cx+off*np.cos(ang)/pix; ys=cy+off*np.sin(ang)/pix
            sp=sample(T,xs,ys)
            pk=[stats(s) for s in sp]
            og[str(off)]={'peak_mean':float(np.mean([p['peak_RJ'] for p in pk])),'peak_min':float(np.min([p['peak_RJ'] for p in pk])),'peak_max':float(np.max([p['peak_RJ'] for p in pk])),
                'HCN_mean_FWHM_mean':float(np.nanmean([p['HCN']['T_mean_FWHM'] for p in pk])),'HCN_mean_FWHM_min':float(np.nanmin([p['HCN']['T_mean_FWHM'] for p in pk])),'HCN_mean_FWHM_max':float(np.nanmax([p['HCN']['T_mean_FWHM'] for p in pk]))}
        rec['offset_grid']=og
        # Monte Carlo of core centre
        mc={}
        for sg in SIGPOS:
            n=1500; xs=cx+rng.normal(0,sg,n)/pix; ys=cy+rng.normal(0,sg,n)/pix
            sp=sample(T,xs,ys); pk=[stats(s) for s in sp]
            def q(f): a=np.array([f(p) for p in pk]); return [float(x) for x in np.nanpercentile(a,[16,50,84])]
            mc[str(sg)]={'peak_RJ':q(lambda p:p['peak_RJ']),'HCN_T_at_v0':q(lambda p:p['HCN']['T_at_v0']),'HCN_mean_FWHM':q(lambda p:p['HCN']['T_mean_FWHM']),
                'HCN_gauss_wt':q(lambda p:p['HCN']['T_gauss_wt']),'HCO+_mean_FWHM':q(lambda p:p['HCO+']['T_mean_FWHM']),'HCO+_gauss_wt':q(lambda p:p['HCO+']['T_gauss_wt'])}
        rec['mc']=mc
        # CO peak (max over region within 0.3" and within 1 beam of core, any channel in window)
        yy,xx=np.mgrid[:T.shape[1],:T.shape[2]]; rr=np.hypot(xx-cx,yy-cy)*pix
        for rad,key in [(0.3,'max_r0.3'),(max(b,0.3),'max_1beam')]:
            m=rr<rad; sub_T=T[isl][:,m]; k,j=np.unravel_index(np.argmax(sub_T),sub_T.shape)
            rec[key]={'T_RJ':float(sub_T[k,j]),'T_Pl':float(TPl(sub_T[k,j]/TRJ(1,nus[isl][k],b,b),nus[isl][k],b,b)),'v':float(vs[isl][k]),'dx':float((xx[m][j]-cx)*pix),'dy':float((yy[m][j]-cy)*pix)}
        cr['beams'][str(b)]=rec
        print(nm,b,'peak',round(rec['at_core_interp']['peak_RJ'],1),'@',round(rec['at_core_interp']['v_peak']),
              'HCN meanFWHM',round(rec['at_core_interp']['HCN']['T_mean_FWHM'],1),'cov',round(rec['at_core_interp']['HCN']['cov_FWHM'],2),'gw',round(rec['at_core_interp']['HCN']['T_gauss_wt'],1),
              'HCO+ meanFWHM',round(rec['at_core_interp']['HCO+']['T_mean_FWHM'],1),'MC.07 peak',[round(x,1) for x in mc['0.07']['peak_RJ']],'HCNmean',[round(x,1) for x in mc['0.07']['HCN_mean_FWHM']],
              'max1beam',round(rec['max_1beam']['T_RJ'],1),flush=True)
    cr['v']=vs[isl].tolist()
    res['cores'][nm]=cr
json.dump(res,open('/workspace/out/q23b/vmatch_pos.json','w'),indent=1)
