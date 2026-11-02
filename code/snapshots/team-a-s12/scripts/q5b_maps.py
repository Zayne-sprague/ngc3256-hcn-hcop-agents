# Q5b: Xb4b native moment 0/1/2 maps (HCN, HCO+ 3-2), local maxima (clump census), figure.
# python3 -I scripts/q5b_maps.py -> out/q5/q5b.log, maps_Xb4b.npz, q5b_maps.png
import numpy as np, sys
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load('/workspace/out/q5/cube_native_Xb4b.npz'); O='/workspace/out/q5/'
dx,dy,pix=d['dx'],d['dy'],float(d['pix']); log=open(O+'q5b.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
V0=2780.; W=(2780-400,2780+400)   # km/s window for moments
res={}
for k,nm in (('H','HCN'),('O','HCO+')):
    a=d[k].astype(np.float32); v=d['v'+k]; dv=float(d['dv'+k]); rms=d['rms'+k]; beam=d['beam'+k]
    m=(v>=W[0])&(v<=W[1]); aw=a[m]; vw=v[m]
    # channel mask: > 2.5 sigma after 3-channel boxcar smoothing (dilated by 1), 3x3 spatial no smoothing (keep native)
    sm=ndi.uniform_filter1d(aw,3,axis=0)*np.sqrt(3)
    msk=sm>2.5*rms[None]; msk=ndi.binary_dilation(msk,structure=np.ones((3,1,1)))
    m0_all=aw.sum(0)*dv                      # full +-400 window, no mask (unbiased, flux conserving)
    am=np.where(msk,aw,0)
    m0=am.sum(0)*dv; w0=am.sum(0)
    with np.errstate(all='ignore'):
        m1=(am*vw[:,None,None]).sum(0)/w0
        m2=np.sqrt(np.clip((am*(vw[:,None,None]-m1)**2).sum(0)/w0,0,None))
    n_free=m.sum(); sig0=rms*dv*np.sqrt(n_free)   # formal noise of unmasked window sum, per pixel (uncorrelated channels)
    res[k]=dict(m0_all=m0_all,m0=m0,m1=m1,m2=m2,sig0=sig0,peak=aw.max(0),vpk=vw[aw.argmax(0)])
    pr(nm,'beam',beam,'dv',dv,'channels in window',n_free,'median rms mJy',np.median(rms)*1e3,'median sig0 (Jy/b km/s)',np.median(sig0))
    pr(nm,'peak of m0_all',m0_all.max(),'at',dx[np.unravel_index(m0_all.argmax(),m0_all.shape)],dy[np.unravel_index(m0_all.argmax(),m0_all.shape)])
np.savez(O+'maps_Xb4b.npz',dx=dx,dy=dy,pix=pix,**{f'{k}_{n}':res[k][n] for k in res for n in res[k]})
# clump census: local maxima of m0_all smoothed by 0.1" (suppress noise spikes) with S/N>sn_min and spacing>=0.25"
def census(m0,sig0,sn_min,sep_pix,box=None):
    s=ndi.gaussian_filter(m0,2.0)           # 0.094" sigma
    mx=(s==ndi.maximum_filter(s,size=2*sep_pix+1))&(m0/np.median(sig0)>sn_min)
    iy,ix=np.nonzero(mx); return [(dx[y,x],dy[y,x],m0[y,x],m0[y,x]/np.median(sig0)) for y,x in zip(iy,ix)]
for k,nm in (('H','HCN'),('O','HCO+')):
    sg=np.median(res[k]['sig0'])
    pr(f'\n# {nm}: local maxima of Xb4b native m0_all (+-400 km/s unmasked), noise per pixel {sg:.4f} Jy/b km/s; S/N>10, spacing >=0.25"')
    for x,y,f,sn in sorted(census(res[k]['m0_all'],res[k]['sig0'],10,5),key=lambda t:-t[2]):
        if np.hypot(x,y)<5.5: pr(f'  ({x:+.2f},{y:+.2f})"  m0={f:.3f}  S/N={sn:.0f}')
fig,ax=plt.subplots(2,3,figsize=(15,9.5))
ext=[dx.max()+pix/2,dx.min()-pix/2,dy.min()-pix/2,dy.max()+pix/2]
for i,(k,nm) in enumerate((('H','HCN'),('O','HCO+'))):
    r=res[k]; ok=r['m0']/np.median(r['sig0'])>5
    ax[i,0].imshow(np.arcsinh(r['m0_all']/0.3),origin='lower',extent=ext,cmap='magma'); ax[i,0].set_title(nm+' m0 (arcsinh)')
    im=ax[i,1].imshow(np.where(ok,r['m1'],np.nan),origin='lower',extent=ext,cmap='RdBu_r',vmin=2650,vmax=2900); plt.colorbar(im,ax=ax[i,1]); ax[i,1].set_title(nm+' m1')
    im=ax[i,2].imshow(np.where(ok,r['m2'],np.nan),origin='lower',extent=ext,cmap='viridis',vmin=0,vmax=120); plt.colorbar(im,ax=ax[i,2]); ax[i,2].set_title(nm+' m2')
for a in ax.ravel(): a.set_xlim(3,-3); a.set_ylim(-5,5)
plt.savefig(O+'q5b_maps.png',dpi=90)
