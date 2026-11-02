# Q23: CO(2-1) peak brightness temperature at the N and S dense-gas core positions, archive cube 2015.1.00902.S (0.73")
import numpy as np, json, sys, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
C=299792.458; H=6.62607e-27; KB=1.380649e-16; CL=2.99792458e10
P='/archive/files/2015.1.00902.S/member.uid___A001_X2fb_X881.ari_l.ngc_3256_sci.spw2_228401MHz.12m.cube.I.pbcor.fits'
RA0,DE0=156.96346,-43.90461
NOFF=float(sys.argv[1]) if len(sys.argv)>1 else 550.
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}   # arcsec offsets from phase centre, Q5 0.24" HCN/HCO+ Gaussian-fit centres (out/q5/q5g_clumps.csv, q5d_fits.csv)
h=fits.open(P,memmap=True); hd=h[0].header; cube=h[0].data[0]
F0=230.538e9; nu=hd['CRVAL3']+(np.arange(hd['NAXIS3']))*hd['CDELT3']
v=(1-nu/F0)*C   # radio-def LSRK velocities rel. to rest 230.538 (not sys-shifted)
vsys=2780.; fl=(1-vsys/C)
dv=abs(np.median(np.diff(v)))
bmaj,bmin=hd['BMAJ']*3600,hd['BMIN']*3600
w=WCS(hd).celestial; pix=abs(hd['CDELT2'])*3600
def TRJ(Sjy,nu0): # Jy/beam -> K  (Rayleigh-Jeans), beam in arcsec^2: T=1.222e3*S[mJy]/(nu[GHz]^2 th^2) K
    return 1.222e3*Sjy*1e3/((nu0/1e9)**2*bmaj*bmin)
def TPl(Sjy,nu0):
    I=Sjy/(1.1331*bmaj*bmin*(np.pi/180/3600)**2)*1e-23
    return H*nu0/KB/np.log1p(2*H*nu0**3/CL**2/I)
print('beam',bmaj,bmin,hd['BPA'],'pix',pix,'dv',dv,'nch',len(v), 'RJ conv 1 Jy/b=',TRJ(1,228.4e9))
# noise: channels free of line. line spans vsys +-~500 km/s -> restfreq velocity centre = vsys
vc=v-(0)  # rest-freq velocity; the line sits near vsys (2780)
off=np.abs(vc-vsys)>NOFF
print('line-free channels',off.sum(), 'v range',v.min(),v.max())
res={'beam':[bmaj,bmin,hd['BPA']],'dv':dv,'pix':pix,'cores':{}}
ny,nx=cube.shape[1:]
Y,X=np.mgrid[:ny,:nx]; ra,de=w.pixel_to_world_values(X,Y)
dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600
rms_ch=np.nanstd(cube[off],axis=0)   # per-pixel noise Jy/beam
for nm,(x0,y0) in CORES.items():
    r=np.hypot(dx-x0,dy-y0); iy,ix=np.unravel_index(np.argmin(r),r.shape)
    sp=np.array(cube[:,iy,ix],float); sig=np.nanstd(sp[off])
    k=np.nanargmax(sp)
    # local max within one beam (0.5 beam radius) too
    m=r<0.5*bmaj; sub=np.nan_to_num(cube[:, m],nan=-9)
    kk,jj=np.unravel_index(np.argmax(sub),sub.shape)
    bj=np.where(m)
    d={'pix':(int(iy),int(ix)),'offset_arcsec':float(r[iy,ix]),'peak_Jy':float(sp[k]),'v_peak':float(v[k]),'noise_Jy':float(sig),
       'Tb_RJ':float(TRJ(sp[k],nu[k])),'Tb_Planck':float(TPl(sp[k],nu[k])),'e_Tb_RJ':float(TRJ(sig,nu[k])),
       'peak_in_halfbeam_Jy':float(sub[kk,jj]),'Tb_RJ_halfbeam':float(TRJ(sub[kk,jj],nu[kk])),'v_halfbeam':float(v[kk]),
       'M0_Jykms':float(np.nansum(sp[~off])*dv),'sp':sp.tolist()}
    # mom0-based peak over velocity in 3-ch smoothed (to compare with Brunetti channel peak)
    res['cores'][nm]=d
    print(nm,{a:b for a,b in d.items() if a!='sp'})
# also measure the max Tb over the whole map & the disk percentile for context
pk=np.nanmax(cube,axis=0); T=TRJ(pk,228.4e9)
res['map_peak_K']=float(np.nanmax(T))
for q in [16,50,84]: res[f'p{q}_K_pixels_gt3sig']=float(np.nanpercentile(T[pk>3*np.nanmedian(rms_ch)],q))
iy,ix=np.unravel_index(np.nanargmax(pk),pk.shape); res['map_peak_pos_arcsec']=[float(dx[iy,ix]),float(dy[iy,ix])]
print({k:res[k] for k in res if k.startswith('p') or k.startswith('map')})
res['noff']=int(off.sum())
json.dump(res,open('/workspace/out/q23/co_tb.json','w'),indent=1)
np.savez('/workspace/out/q23/co_peakmap.npz',T=T,dx=dx,dy=dy)
