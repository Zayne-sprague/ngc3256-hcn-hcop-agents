# Geometry of the Xb4d blue-wing (-450..-150) emission around S: centroid, extent, PA of the S->blob vector; figure of blue, pseudo-wing, Xb4b convolved
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import ndimage
res=[]
fig,ax=plt.subplots(2,3,figsize=(15,9))
for r,line in enumerate(('HCO+','HCN')):
    d=np.load(f'/workspace/out/s11_Xb4d_{line}_maps.npz'); s=np.load(f'/workspace/out/s11_Xb4b_{line}_smoothed_maps.npz')
    hd=cube_info(files('*Xb4d.ngc_3256_sci.spw25.cube.I.pbcor.fits')[0])[1]; wc=WCS(hd).celestial
    ps=0.15; x0,y0=int(d['x0']),int(d['y0']); B=d['blue']; sig=1.4826*np.median(abs(d['pseudoB']-np.median(d['pseudoB'])))
    yy,xx=np.mgrid[:B.shape[0],:B.shape[1]]
    ra,dec=wc.pixel_to_world_values(xx+x0,yy+y0)
    S=REGIONS['S']; dE=(ra-S[0])*np.cos(np.deg2rad(S[1]))*3600; dN=(dec-S[1])*3600
    # blob: connected region >4 sigma containing the max within 4" S of S nucleus
    msk=(B>4*sig)&(np.hypot(dE,dN)<7)
    lab,n=ndimage.label(msk)
    for l in range(1,n+1):
        m=lab==l
        if m.sum()<8: continue
        w=B[m]; cE=(dE[m]*w).sum()/w.sum(); cN=(dN[m]*w).sum()/w.sum()
        cov=np.cov(np.vstack([dE[m],dN[m]]),aweights=w); ev,evec=np.linalg.eigh(cov); pa=np.degrees(np.arctan2(evec[0,1],evec[1,1]))%180
        res.append((line,l,m.sum()*ps**2,w.sum()*ps**2/ (0.93*0.74*1.1331) ,cE,cN,np.hypot(cE,cN),np.degrees(np.arctan2(cE,cN))%360,pa,np.sqrt(ev[1])*2.355,np.sqrt(ev[0])*2.355))
    ext=[dE.max(),dE.min(),dN.min(),dN.max()]
    for c,(m_,t) in enumerate([(B,'Xb4d blue wing'),(d['pseudoB'],'Xb4d pseudo-wing -900..-600'),(s['blue'],'Xb4b conv. blue wing')]):
        a=ax[r,c]; im=a.imshow(m_,origin='lower',extent=ext,cmap='magma',vmin=-0.1,vmax=0.3); a.set_title(f'{line} {t} [Jy/b km/s]'); a.set_xlim(8,-8); a.set_ylim(-9,5)
        a.contour(dE,dN,B/sig,[3,5,8],colors='c',linewidths=.5)
        a.plot(0,0,'w+'); plt.colorbar(im,ax=a)
plt.tight_layout(); plt.savefig('/workspace/out/s11_wing_maps.png',dpi=60)
print('line label area_arcsec2 flux_beams(Jy km/s sum/beam) cE cN r PA_centroid(deg E of N) PA_major FWHM_major FWHM_minor [arcsec relative to S]')
for x in res: print(x[0],x[1],'%.2f %.2f %.2f %.2f %.2f %.0f %.0f %.2f %.2f'%x[2:])
