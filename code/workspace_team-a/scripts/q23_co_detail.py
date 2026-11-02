import numpy as np, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load('/workspace/out/q23/co_peakmap.npz'); T,dx,dy=d['T'],d['dx'],d['dy']
r=json.load(open('/workspace/out/q23/co_tb.json'))
from astropy.io import fits
P='/archive/files/2015.1.00902.S/member.uid___A001_X2fb_X881.ari_l.ngc_3256_sci.spw2_228401MHz.12m.cube.I.pbcor.fits'
h=fits.open(P,memmap=True); hd=h[0].header; cube=h[0].data[0]
nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/230.538e9)*299792.458
conv=r['cores']['N']['Tb_RJ']/r['cores']['N']['peak_Jy']
CORES={'N':(-0.06,2.65),'S':(-0.15,-2.59)}
fig,ax=plt.subplots(1,3,figsize=(15,4.5))
im=ax[0].imshow(T,origin='lower',extent=[dx.min(),dx.max(),dy.min(),dy.max()] if False else None,vmax=20,cmap='viridis')
iy,ix=np.indices(T.shape)
# ring statistics: Tb_peak vs offset from core position
out={}
for nm,(x0,y0) in CORES.items():
    rr=np.hypot(dx-x0,dy-y0)
    sp=cube[:,rr<0.1].mean(axis=1) if False else None
    print(nm,'peak-T (K, RJ) vs radius from core:')
    for a,b in [(0,.1),(.1,.2),(.2,.4),(.4,.6),(.6,.8),(.8,1.2),(1.2,2.0)]:
        m=(rr>=a)&(rr<b); print(f'  r={a}-{b}\" n={m.sum()} median={np.nanmedian(T[m]):.1f} max={np.nanmax(T[m]):.1f} min={np.nanmin(T[m]):.1f}')
        out[f'{nm}_{a}_{b}']=[float(np.nanmedian(T[m])),float(np.nanmax(T[m]))]
    ax[0].plot(*np.array([[np.argmin(abs(dx[0]-x0))],[np.argmin(abs(dy[:,0]-y0))]]),'rx')
    ax[1 if nm=='N' else 2].plot(v,r['cores'][nm]['sp'] and np.array(r['cores'][nm]['sp'])*conv,'k-',label=f'{nm} core pixel')
    # best-in-halfbeam spectrum
    m=np.hypot(dx-x0,dy-y0)<0.41; jj=np.unravel_index(np.nanargmax(np.nanmax(cube[:,m] if False else np.where(m,cube,-9),axis=0)),m.shape)
    a=ax[1 if nm=='N' else 2]; a.plot(v,cube[:,jj[0],jj[1]]*conv,'r-',label=f'max pixel within 0.4" ({dx[jj]:.2f},{dy[jj]:.2f})')
    a.set_xlabel('v_rest-radio (km/s)'); a.set_ylabel('T_b RJ (K)'); a.legend(); a.set_title(nm)
fig.savefig('/workspace/out/q23/co_core_spectra.png',dpi=90)
json.dump(out,open('/workspace/out/q23/co_ring.json','w'))
