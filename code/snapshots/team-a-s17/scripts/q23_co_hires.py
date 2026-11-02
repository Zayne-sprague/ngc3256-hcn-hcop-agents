# Q23: CO(2-1) T_b at N and S cores from the 0.18" cube (2015.1.00714.S 12m spw1), native and convolved to circular beams
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel, convolve_fft
C=299792.458; H=6.62607e-27; KB=1.380649e-16; CL=2.99792458e10
P='/archive/files/2015.1.00714.S/member.uid___A001_X2d6_Xae.ari_l.NGC_3256_sci.spw1_229288MHz.12m.cube.I.pbcor.fits'
RA0,DE0=156.96346,-43.90461
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
h=fits.open(P,memmap=True); hd=h[0].header; cube=h[0].data[0]
F0=230.538e9; nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/F0)*C
bmaj,bmin,pa=hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']; pix=abs(hd['CDELT2'])*3600
print('beam',bmaj,bmin,pa,'pix',pix,'v',v.min(),v.max(),'dv',np.median(np.diff(v)))
w=WCS(hd).celestial
HW=int(round(6/pix))  # subcube +-6"
def cen(x0,y0):
    ra=RA0+x0/3600/np.cos(np.deg2rad(DE0)); de=DE0+y0/3600
    x,y=w.world_to_pixel_values(ra,de); return int(round(float(x))),int(round(float(y)))
sel=(v>2300)&(v<2844); ks=np.where(sel)[0]; off=np.where(v<2300)[0]  # line-free: v<2300 (+/-480 of vsys)
res={'beam':[bmaj,bmin,pa],'pix':pix,'dv':float(abs(np.median(np.diff(v)))),'cores':{}}
def TRJ(S,nu0,a,b): return 1.222e3*S*1e3/((nu0/1e9)**2*a*b)
def TPl(S,nu0,a,b):
    I=S/(1.1331*a*b*(np.pi/180/3600)**2)*1e-23
    return H*nu0/KB/np.log1p(2*H*nu0**3/CL**2/I)
for nm,(x0,y0) in CORES.items():
    ix,iy=cen(x0,y0)
    sub=np.array(cube[:,iy-HW:iy+HW+1,ix-HW:ix+HW+1],float)  # (nch,ny,nx)
    c0=HW
    d={'center_pix':[ix,iy],'beams':{}}
    for lab,fwhm in [('native',None),('0.26',0.26),('0.38',0.38),('0.57',0.57),('0.75',0.75),('1.0',1.0),('1.7',1.7)]:
        if fwhm is None:
            S=sub; a,b=bmaj,bmin
            # keep native; sample at the centre and best in half-beam
        else:
            # convolve to circular beam fwhm (approx: deconvolve geometric-mean native beam; elliptical beam residual ignored)
            g=np.sqrt(max(fwhm**2-bmaj*bmin,1e-6)); 
            if fwhm<bmaj: continue
            sig=g/2.3548/pix
            k=Gaussian2DKernel(sig); 
            # brightness conserved in Jy/beam: scale by (new beam area / old) via kernel normalisation x area ratio
            fac=fwhm**2/(bmaj*bmin)
            S=np.array([convolve_fft(sub[i],k,boundary='fill',fill_value=0,allow_huge=True) for i in range(sub.shape[0])])*fac
            a=b=fwhm
        sp=S[:,c0,c0]; rms=np.nanstd(sp[off]) if True else 0
        ksel=ks; kpk=ks[np.nanargmax(sp[ks])]
        m=np.hypot(*np.mgrid[-HW:HW+1,-HW:HW+1])*pix<0.5*min(a,b) if fwhm is None else None
        rec={'beam':[a,b],'peak_Jy':float(sp[kpk]),'v_peak':float(v[kpk]),'rms_Jy':float(rms),
             'Tb_RJ':float(TRJ(sp[kpk],nu[kpk],a,b)),'Tb_Planck':float(TPl(sp[kpk],nu[kpk],a,b)),'e_Tb_RJ':float(TRJ(rms,nu[kpk],a,b))}
        if m is not None:
            ss=S[ks][:,m]; kk,jj=np.unravel_index(np.nanargmax(ss),ss.shape)
            rec['Tb_RJ_halfbeam']=float(TRJ(ss[kk,jj],nu[ks[kk]],a,b)); rec['Tb_Pl_halfbeam']=float(TPl(ss[kk,jj],nu[ks[kk]],a,b))
        # max T over +-1 beam box too (for CO peaks offset from dense-gas centre)
        r=np.hypot(*np.mgrid[-HW:HW+1,-HW:HW+1])*pix<max(a,b)
        rec['Tb_RJ_max_1beam']=float(TRJ(np.nanmax(S[ks][:,r]),nu[kpk],a,b))
        d['beams'][lab]=rec
        print(nm,lab,{k:round(x,3) if isinstance(x,float) else x for k,x in rec.items()})
        if lab in('native','0.38'): d.setdefault('spec',{})[lab]=sp[ks].tolist()
    d['v']=v[ks].tolist()
    res['cores'][nm]=d
json.dump(res,open('/workspace/out/q23/co_hires.json','w'),indent=1)
