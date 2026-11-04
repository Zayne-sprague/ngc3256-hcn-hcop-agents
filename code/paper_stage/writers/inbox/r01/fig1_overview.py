"""R01 / Figure 1: data overview at the common 0.26" (50 pc) beam, Xb4b (2017.1.00379.S).
Pipeline = s170_q26_ratios.py: pbcor cube cutout -> RJ T_b per channel -> flux-conserving convolution of T_b to 0.26" circular (native 0.246x0.232" PA -49)
-> baseline = per-pixel median of line-free channels (v_rel<-450, >600 km/s) -> HCN interpolated onto the HCO+ velocity grid.
Panels: (a,b) peak T_b = max over |v-2780|<350 km/s of the 3-channel-smoothed T_b (K, RJ), masked <3 sigma (per-pixel line-free rms of the smoothed spectrum);
(c) ratio of moment-0 (K km/s, window +-350 km/s about 2780), masked where either moment-0 < 3 sigma. Equal beams => T_b ratio.
usage: python3 -I fig1_overview.py"""
import sys, json, numpy as np, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
C=299792.458; V0=2780.; RA0,DE0=156.96346,-43.90461; TH=0.26; HALF=7.0
LINES={'HCO+':('spw25',267.5576),'HCN':('spw27',265.8864)}
def tbk(nu,a,b): return 1.222e3/(nu**2*a*b)
def cov(maj,mn,pa):
    s1,s2=maj/2.3548,mn/2.3548; pa=np.deg2rad(pa); u=np.array([-np.sin(pa),np.cos(pa)]); w=np.array([np.cos(pa),np.sin(pa)])
    return s1**2*np.outer(u,u)+s2**2*np.outer(w,w)
def kernel(Ct,Cn,pix):
    Ck=(Ct-Cn)/pix**2; ev=np.linalg.eigvalsh(Ck); assert ev.min()>-1e-9
    sg=np.sqrt(np.sort(ev)[::-1]); R=int(np.ceil(5*sg[0]))+1; y,x=np.mgrid[-R:R+1,-R:R+1]
    Ci=np.linalg.inv(Ck+1e-9*np.eye(2)); q=Ci[0,0]*x*x+2*Ci[0,1]*x*y+Ci[1,1]*y*y; k=np.exp(-0.5*q); return k/k.sum()
def load(spw,f0):
    h=fits.open('/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.%s.cube.I.pbcor.fits'%spw,memmap=True); hd=h[0].header; d=h[0].data
    while d.ndim>3: d=d[0]
    nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']; v=C*(1-nu/(f0*1e9))-V0
    pix=abs(hd['CDELT2'])*3600; w=WCS(hd).celestial; cx,cy=[float(t) for t in w.world_to_pixel_values(RA0,DE0)]; n=int(HALF/pix)
    ix,iy=int(round(cx)),int(round(cy)); y0=iy-n; x0=ix-n
    cut=np.array(d[:,y0:iy+n+1,x0:ix+n+1],dtype=np.float32)*1e3
    bm=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); nuo=f0*(1-V0/C)
    return cut*np.float32(tbk(nuo,bm[0],bm[1])),v,pix,bm,(cx-x0,cy-y0)
sm={}
for ln,(spw,f0) in LINES.items():
    Tb,v,pix,bm,cen=load(spw,f0); k=kernel(cov(TH,TH,0),cov(*bm),pix)
    c=np.empty_like(Tb)
    for i in range(Tb.shape[0]): c[i]=fftconvolve(np.nan_to_num(Tb[i]),k,'same')
    lf=(v<-450)|(v>600); c-=np.median(c[lf],axis=0); sm[ln]=(c,v)
    if ln=='HCO+': vh=v; pix_,cen_=pix,cen
keep=vh>-1400; cH,vH=sm['HCN']; fi=np.interp(vh[keep],vH[::-1],np.arange(len(vH))[::-1]); i0=np.clip(np.floor(fi).astype(int),0,len(vH)-2); fr=(fi-i0).astype(np.float32)[:,None,None]
cube={'HCN':cH[i0]*(1-fr)+cH[i0+1]*fr,'HCO+':sm['HCO+'][0][keep]}; v=vh[keep]; dv=abs(v[1]-v[0]); lf=(v<-450)|(v>600); pix=pix_
ny,nx=cube['HCO+'].shape[1:]
# peak T_b maps
sel=np.abs(v)<350; Pk={}; Pe={}
for ln,c in cube.items():
    s3=np.apply_along_axis(lambda a:np.convolve(a,np.ones(3)/3,'same'),0,c)
    Pe[ln]=np.std(s3[lf],axis=0); Pk[ln]=s3[sel].max(axis=0)
# moment-0
kw=np.clip(np.minimum(v+dv/2,350)-np.maximum(v-dv/2,-350),0,dv)/dv; sk=np.sqrt((kw**2).sum())*dv
M={ln:np.tensordot(kw,c,axes=(0,0))*dv for ln,c in cube.items()}; S={ln:np.std(c[lf],axis=0)*sk for ln,c in cube.items()}
ok=(M['HCN']>3*S['HCN'])&(M['HCO+']>3*S['HCO+']); Rm=np.where(ok,M['HCN']/M['HCO+'],np.nan)
okp={ln:Pk[ln]>3*Pe[ln] for ln in Pk}
np.savez_compressed('/workspace/paper_answers/r01/fig1_maps.npz',pkHCN=Pk['HCN'],pkHCO=Pk['HCO+'],R=Rm,pix=pix,cen=cen_)
# axes: offsets (E positive, arcsec) from RA0,DE0
cx,cy=cen_; xo=-(np.arange(nx)-cx)*pix; yo=(np.arange(ny)-cy)*pix   # east = left in image => +RA offset = -x ... use x offset positive east
ext=[xo[0]+pix/2,xo[-1]-pix/2,yo[0]-pix/2,yo[-1]+pix/2]
XL=(3.6,-3.6); YL=(-5.4,5.4)
plt.rcParams.update({'font.size':10})
fig,axs=plt.subplots(1,3,figsize=(12.5,6.4),sharey=True)
panels=[(np.where(okp['HCO+'],Pk['HCO+'],np.nan),'inferno',0,3.0,r'HCO$^+$ J=3–2 peak $T_b$ (K)'),(np.where(okp['HCN'],Pk['HCN'],np.nan),'inferno',0,3.0,r'HCN J=3–2 peak $T_b$ (K)'),(Rm,'viridis',0,1.3,r'HCN/HCO$^+$ $T_b$ ratio')]
for i,(ax,(im,cm,lo,hi,lab)) in enumerate(zip(axs,panels)):
    cmap=plt.get_cmap(cm).copy(); cmap.set_bad('0.85' if i==2 else 'k')
    h=ax.imshow(im,origin='lower',extent=ext,cmap=cmap,vmin=lo,vmax=hi,interpolation='nearest'); ax.set_xlim(*XL); ax.set_ylim(*YL)
    cb=fig.colorbar(h,ax=ax,orientation='horizontal',pad=0.14,fraction=0.045); cb.set_label(lab)
    col='w' if i<2 else 'k'
    for nm,(dx,dy) in {'N':(-0.15,2.55),'S':(-0.15,-2.55)}.items():
        ax.plot([-dx],[dy],'+',color=col if i<2 else 'w',ms=11,mew=1.6) if False else ax.plot([dx],[dy],'+',color='#56B4E9' if i<2 else 'w',ms=12,mew=1.8)
        ax.text(dx+0.35,dy+0.2 if nm=='N' else dy-0.45,nm,color='#56B4E9' if i<2 else 'w',weight='bold',fontsize=11,ha='right')
    ax.plot([0.0],[0.0],'x',color='#56B4E9' if i<2 else 'w',ms=8,mew=1.4); ax.text(-0.35,0.0,'Bridge',color='#56B4E9' if i<2 else 'k',fontsize=9,va='center')
    ax.add_patch(Ellipse((3.0,-4.7),TH,TH,fc='w',ec='k',lw=1,zorder=5)); ax.text(3.4,-4.15,'beam 0.26" (50 pc)',ha='left',fontsize=8,color='w' if i<2 else 'k')
    ax.set_xlabel(r'$\Delta$RA (arcsec, east to the left)'); ax.set_title(lab if i<2 else 'ratio (>3σ both lines)',fontsize=10)
    ax.xaxis.set_ticks([3,2,1,0,-1,-2,-3])
    s=ax.secondary_xaxis('top',functions=(lambda x:x*194,lambda x:x/194)); s.set_xlabel('pc (1" = 194 pc)')
axs[0].set_ylabel(r'$\Delta$Dec (arcsec)'); fig.text(0.5,0.005,'offsets from RA 156.96346°, Dec −43.90461° (J2000)',ha='center',fontsize=8)
sy=axs[2].secondary_yaxis('right',functions=(lambda x:x*194,lambda x:x/194)); sy.set_ylabel('pc')
fig.tight_layout(); fig.savefig('/workspace/paper_answers/r01/fig1_data_overview.png',dpi=180)
# numbers for caption
def at(m,dx,dy,r=0.13):
    X,Y=np.meshgrid(xo,yo); s=np.hypot(X-dx,Y-dy)<=r; return np.nanmedian(m[s])
for nm,(dx,dy) in {'N':(-0.15,2.55),'S':(-0.15,-2.55),'Bridge':(0,0)}.items():
    print(nm,'peak HCO+ %.2f HCN %.2f ratio(m0) %.2f'%(at(Pk['HCO+'],dx,dy),at(Pk['HCN'],dx,dy),at(Rm,dx,dy)))
j=np.unravel_index(np.nanargmax(np.where(yo[:,None]>0,M['HCO+'],0)),M['HCO+'].shape); print('HCO+ m0 peak N offset (E,N) %.2f %.2f'%(xo[j[1]],yo[j[0]]))
j=np.unravel_index(np.nanargmax(np.where(yo[:,None]<-1,M['HCO+'],0)),M['HCO+'].shape); print('HCO+ m0 peak S offset (E,N) %.2f %.2f'%(xo[j[1]],yo[j[0]]))
print('typical sigma_peak HCO+ %.3f HCN %.3f K; sigma m0 %.2f %.2f'%(np.median(Pe['HCO+']),np.median(Pe['HCN']),np.median(S['HCO+']),np.median(S['HCN'])))
