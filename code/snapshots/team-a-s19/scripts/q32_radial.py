# Q32c: azimuthally averaged radial profile (native 0.18" cube) about the CO peak near each core vs the CLEAN beam; plots map
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('scripts/q32_map_fit.py').read().split("out={}")[0].split("h=fits.open")[0])
h=fits.open(P,memmap=True); hd=h[0].header; cube=h[0].data[0]
nu=hd['CRVAL3']+np.arange(hd['NAXIS3'])*hd['CDELT3']; v=(1-nu/230.538e9)*C
bmaj,bmin,pa=hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']; pix=abs(hd['CDELT2'])*3600
w=WCS(hd).celestial; HW=int(round(2.0/pix)); K=1.222e3/((nu/1e9)**2*bmaj*bmin)*1e3
res={}; fig,ax=plt.subplots(2,3,figsize=(14,8))
for r_,(nm,(x0,y0)) in enumerate(CORES.items()):
    ra=RA0+x0/3600/np.cos(np.deg2rad(DE0)); de=DE0+y0/3600
    ix,iy=[int(round(float(a))) for a in w.world_to_pixel_values(ra,de)]
    sub=np.array(cube[:,iy-HW:iy+HW+1,ix-HW:ix+HW+1],float)*K[:,None,None]
    v0,fw=WIN[nm]; ks=np.where((v>v0-fw/2)&(v<v0+fw/2)&(v<2844))[0]
    kk=np.where((v>2300)&(v<2844))[0]
    for lab,m in [('peak_ch',None),('win',None)]:
        if lab=='peak_ch':
            k=kk[np.argmax(sub[kk][:,HW-4:HW+5,HW-4:HW+5].max((1,2)))]; M=sub[k]
        else: M=sub[ks].mean(0)
        # centre: max within 0.15" of nominal core
        Y,X=np.mgrid[-HW:HW+1,-HW:HW+1]; rr=np.hypot(X,Y)*pix
        mm=np.where(rr<0.15,M,-1e9); jy,jx=np.unravel_index(mm.argmax(),M.shape)
        R=np.hypot(X-jx+HW,Y-jy+HW)*pix; edges=np.array([0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.7,0.9,1.2,1.6])
        prof=[M[(R>=a)&(R<b)].mean() for a,b in zip(edges[:-1],edges[1:])]; rc=0.5*(edges[1:]+edges[:-1])
        bg=np.sqrt(bmaj*bmin); beam=M[jy,jx]*np.exp(-4*np.log(2)*rc**2/bg**2)
        # half-power radius of observed profile vs beam
        p=np.array(prof); half=np.interp(-M[jy,jx]/2,-p,rc)
        res[f'{nm}_{lab}']={'peak_K':float(M[jy,jx]),'peak_offset':[float((jx-HW)*pix),float((jy-HW)*pix)],'r':rc.tolist(),'prof':p.tolist(),'beam_prof':beam.tolist(),'r_half_obs':float(half),'r_half_beam':float(bg/2)}
        print(nm,lab,'peak',round(M[jy,jx],1),'K  r_half obs',round(half,3),'beam',round(bg/2,3)); print('  prof',np.round(p,1)); print('  beam',np.round(beam,1))
        a=ax[r_,0 if lab=='peak_ch' else 1]; a.plot(rc,p,'o-',label='obs'); a.plot(rc,beam,'k--',label='beam'); a.set_title(f'{nm} {lab}'); a.legend()
    ax[r_,2].imshow(sub[k],origin='lower',extent=[-HW*pix,HW*pix,-HW*pix,HW*pix]); 
plt.savefig('out/q32/radial.png',dpi=90); json.dump(res,open('out/q32/radial.json','w'),indent=1)
