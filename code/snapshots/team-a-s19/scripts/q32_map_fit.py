# Q32b: 2D elliptical Gaussian + constant fit to native 0.18" CO(2-1) maps at S and N cores; deconvolve beam.
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
C=299792.458
P='/archive/files/2015.1.00714.S/member.uid___A001_X2d6_Xae.ari_l.NGC_3256_sci.spw1_229288MHz.12m.cube.I.pbcor.fits'
RA0,DE0=156.96346,-43.90461
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
WIN={'S':(2711,100),'N':(2806,155)}   # HCN v0,FWHM (Q5, 0.25" aperture)
h=fits.open(P,memmap=True); hd=h[0].header; cube=h[0].data[0]
nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/230.538e9)*C
bmaj,bmin,pa=hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']; pix=abs(hd['CDELT2'])*3600
w=WCS(hd).celestial; HW=int(round(1.5/pix))
K=1.222e3/((nu/1e9)**2*bmaj*bmin)*1e3   # K per Jy/beam, RJ
out={}
def g2(par,X,Y,thm=None):
    A,x0,y0,smaj,smin,phi,c=par
    t=np.deg2rad(phi); xr=(X-x0)*np.cos(t)+(Y-y0)*np.sin(t); yr=-(X-x0)*np.sin(t)+(Y-y0)*np.cos(t)
    return A*np.exp(-4*np.log(2)*((xr/smaj)**2+(yr/smin)**2))+c
for nm,(x0,y0) in CORES.items():
    ra=RA0+x0/3600/np.cos(np.deg2rad(DE0)); de=DE0+y0/3600
    ix,iy=[int(round(float(a))) for a in w.world_to_pixel_values(ra,de)]
    sub=np.array(cube[:,iy-HW:iy+HW+1,ix-HW:ix+HW+1],float)*K[:,None,None]
    v0,fw=WIN[nm]; ks=np.where((v>v0-fw/2)&(v<v0+fw/2)&(v<2844))[0]
    # peak channel at the core pixel within +-0.3"
    kpk=np.where((v>2300)&(v<2844))[0]; kpk=kpk[np.argmax(sub[kpk][:,HW,HW])]
    maps={'peak_ch':sub[kpk],'HCNwin_mean':sub[ks].mean(0)}
    Y,X=(np.mgrid[-HW:HW+1,-HW:HW+1])*pix   # X: +offset in pixel x (=-RA direction), only relative
    for lab,m in maps.items():
        # fit observed Gaussian+const; beam PA is on sky; pixel x axis is east-negative so PA sign flip handled by fitting phi free
        best=None
        for s0 in (0.25,0.4,0.7):
            for ph in (0,60,120):
                try: r=least_squares(lambda p:(g2(p,X,Y)-m).ravel(),[m.max()-m.min(),0,0,s0,s0*0.8,ph,m.min()],
                  bounds=([0,-0.4,-0.4,0.1,0.1,-90,-5],[200,0.4,0.4,2,2,270,60]))
                except Exception: continue
                if best is None or r.cost<best.cost: best=r
        A,xc,yc,smaj,smin,phi,c=best.x
        # beam in the same frame: pixel x axis is -RA; fit phi measured from x axis. Beam PA (N->E) -> phi_b = 90-(-PA)... compute both, report deconvolved sizes for circular-equivalent
        sgeo=np.sqrt(smaj*smin); bgeo=np.sqrt(bmaj*bmin)
        dec_geo=np.sqrt(max(sgeo**2-bgeo**2,0))
        # error from residual rms (white-noise approx; correlated by beam so inflate by beam pixels)
        resid=(g2(best.x,X,Y)-m); npb=1.133*bmaj*bmin/pix**2
        J=best.jac; cov=np.linalg.inv(J.T@J)*resid.var()*npb; e=np.sqrt(np.diag(cov))
        # peak at A+c; T0 intrinsic = A*(smaj*smin)/(dec_maj*dec_min) approx geo
        T0=A*sgeo**2/max(dec_geo**2,1e-4) if dec_geo>0.02 else np.nan
        rec={'A_K':A,'const_K':c,'center_offset_arcsec':[xc,yc],'obs_fwhm':[smaj,smin],'phi':phi,'deconv_geo_fwhm':dec_geo,'e_obs_fwhm':[e[3],e[4]],'T0_intrinsic_RJ':T0,'beam_geo':bgeo,'cost':best.cost,'peak_v':float(v[kpk]) if lab=='peak_ch' else [float(v[ks[0]]),float(v[ks[-1]])]}
        # ff at 0.26" for a Gaussian of deconvolved size s: s^2/(s^2+0.26^2)
        rec['ff_0.26']=dec_geo**2/(dec_geo**2+0.26**2); rec['ff_native']=dec_geo**2/(dec_geo**2+bgeo**2)
        out[f'{nm}_{lab}']=rec; print(nm,lab,json.dumps({k:(np.round(x,3).tolist() if not isinstance(x,str) else x) for k,x in rec.items()}))
json.dump(out,open('out/q32/map_fit.json','w'),indent=1,default=float)
