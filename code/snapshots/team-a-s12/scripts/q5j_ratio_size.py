# Q5j: model-independent N/S core size comparison HCN vs HCO+ at matched resolution (both convolved to a common circular 0.30" beam),
# radial profiles about the HCO+ peak, enclosed-flux half-light radii, HCN/HCO+ ratio vs radius, null errors from line-free sums.
# python3 -I scripts/q5j_ratio_size.py -> out/q5h/q5j.log, q5j_profiles.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q5lib import *
import pandas as pd
from scipy import ndimage as ndi
d=load(); dx=d['dx']; dy=d['dy']; pix=float(d['pix']); log=open('/workspace/out/q5h/q5j.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
def conv_to(img,bm,fw):
    """convolve a map with native beam bm=(maj,min,pa) to circular fwhm fw (Gaussian kernel, flux-per-beam rescaled)"""
    Bc=beamcov(*bm); K=np.eye(2)*(fw/2.3548)**2-Bc
    w,v=np.linalg.eigh(K); assert w.min()>0
    hk=int(4*fw/pix); yy,xx=np.mgrid[-hk:hk+1,-hk:hk+1]*pix; Ki=np.linalg.inv(K); ker=np.exp(-0.5*(Ki[0,0]*xx**2+2*Ki[0,1]*xx*yy+Ki[1,1]*yy**2)); ker/=ker.sum()
    from scipy.signal import fftconvolve
    sc=(fw**2)/(bm[0]*bm[1]); return fftconvolve(img,ker,mode='same')*sc
FW=0.30
mH,_=win_map(d,'H',*VW); mO,_=win_map(d,'O',*VW)
cH=conv_to(mH,d['beamH'],FW); cO=conv_to(mO,d['beamO'],FW)
nH=[conv_to(m,d['beamH'],FW) for m in null_maps_scaled(d,'H',int((VW[1]-VW[0])/float(d['dvH'])))]
nO=[conv_to(m,d['beamO'],FW) for m in null_maps_scaled(d,'O',int((VW[1]-VW[0])/float(d['dvO'])))]
sH=np.mean([m[np.hypot(dx,dy)<4].std() for m in nH]); sO=np.mean([m[np.hypot(dx,dy)<4].std() for m in nO])
pr(f'common beam {FW}" ; null rms (Jy/b km/s) HCN {sH:.4f} HCO+ {sO:.4f}')
edges=np.array([0,0.1,0.2,0.3,0.45,0.6,0.8,1.0,1.5,2.0]); rows=[]
for src,(x,y) in (('N',(-0.05,2.63)),('S',(-0.15,-2.60))):
    m=np.hypot(dx-x,dy-y)<1.0; z=np.where(m,cO,-9); iy,ix=np.unravel_index(z.argmax(),z.shape); xc,yc=dx[iy,ix],dy[iy,ix]
    zh=np.where(m,cH,-9); jy,jx=np.unravel_index(zh.argmax(),zh.shape)
    pr(f'\n{src}: HCO+ peak at ({xc:.3f},{yc:.3f}) ; HCN peak at ({dx[jy,jx]:.3f},{dy[jy,jx]:.3f}) at {FW}" ; peak HCN/HCO+ = {cH[jy,jx]/cO[iy,ix]:.3f}')
    r=np.hypot(dx-xc,dy-yc); cumH=[];cumO=[]
    for a,b in zip(edges[:-1],edges[1:]):
        s=(r>=a)&(r<b); n=s.sum(); IH=cH[s].mean(); IO=cO[s].mean(); nb=n*pix**2/(1.13*FW**2)   # independent beams
        eH=sH/np.sqrt(max(nb,1)); eO=sO/np.sqrt(max(nb,1)); R=IH/IO; eR=R*np.hypot(eH/IH,eO/IO)
        rows.append(dict(src=src,r_in=a,r_out=b,I_HCN=IH,I_HCO=IO,ratio=R,e_ratio=eR)); pr(f'  r {a:.2f}-{b:.2f}: <I> HCN {IH:.3f} HCO+ {IO:.3f}  ratio {R:.3f} +- {eR:.3f}')
    # enclosed flux half-light radius within 1.0" (background-subtracted with the 1.5-2.0" annulus level)
    for nm,im in (('HCN',cH),('HCO+',cO)):
        bg=im[(r>1.5)&(r<2.0)].mean(); rr=np.sort(np.unique(np.round(np.arange(0.05,1.01,0.05),3))); enc=np.array([(im[r<q]-bg).sum() for q in rr]); tot=enc[-1]
        hl=np.interp(0.5*tot,enc,rr); fr=[np.interp(f*tot,enc,rr) for f in (0.25,0.75)]
        pr(f'  {nm}: enclosed (r<1", bg-subtracted {bg:.3f}) flux {tot*pix**2/(1.13*FW**2):.2f} Jy km/s ; r_25,r_50,r_75 = {fr[0]:.2f},{hl:.2f},{fr[1]:.2f} arcsec')
pd.DataFrame(rows).to_csv('/workspace/out/q5h/q5j_profiles.csv',index=False)
# ---- r50 / r75 noise + background systematics: add randomly rolled line-free noise maps, vary the background annulus
rng=np.random.default_rng(3); pr('\n# half-light radius difference r50(HCN)-r50(HCO+) [arcsec]: noise MC (200) and background choice')
for src,(x,y) in (('N',(-0.05,2.63)),('S',(-0.15,-2.60))):
    m=np.hypot(dx-x,dy-y)<1.0; z=np.where(m,cO,-9); iy,ix=np.unravel_index(z.argmax(),z.shape); xc,yc=dx[iy,ix],dy[iy,ix]; r=np.hypot(dx-xc,dy-yc)
    rr=np.arange(0.05,1.01,0.05)
    def r50(im,bg_sel):
        bg=0. if bg_sel is None else im[bg_sel].mean(); enc=np.array([(im[r<q]-bg).sum() for q in rr]); return np.interp(0.5*enc[-1],enc,rr),np.interp(0.75*enc[-1],enc,rr)
    for lab,sel in (('bg 1.5-2.0"',(r>1.5)&(r<2.0)),('bg 1.0-1.5"',(r>1.0)&(r<1.5)),('no bg',None)):
        d50=[];d75=[]
        for t in range(200):
            sh=lambda m:np.roll(m,(rng.integers(-60,60),rng.integers(-60,60)),(0,1))
            a=r50(cH+sh(nH[rng.integers(len(nH))]),sel); b=r50(cO+sh(nO[rng.integers(len(nO))]),sel); d50.append(a[0]-b[0]); d75.append(a[1]-b[1])
        a0=r50(cH,sel); b0=r50(cO,sel)
        pr(f'{src} {lab}: r50 HCN {a0[0]:.3f} HCO+ {b0[0]:.3f} diff {a0[0]-b0[0]:+.3f} (noise sd {np.std(d50):.3f}); r75 diff {a0[1]-b0[1]:+.3f} (sd {np.std(d75):.3f})')
