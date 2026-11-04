"""r12 (ii): statistical errors on the peak T_b (K) of HCN and HCO+ 3-2 at N and S, Xb4b, 50/55/120/250 pc beams (0.26/0.2835/0.6186/1.2887", 194 pc/"),
reproducing the 'Tchan_peak' definition of scripts/s170_q26_ratios.py (HCO+ +-350 km/s m0 peak pixel within 0.4" of nominal; per-pixel spectrum smoothed with a 3-channel
boxcar; maximum over |v|<350 km/s) but with the noise of THAT pixel's smoothed spectrum (std of line-free channels v<-450 or >600), instead of the aperture-mean noise s170 used.
Also: noise-bias of taking a maximum (Monte Carlo: add noise-only realisations drawn by circularly shifting the line-free baseline? no -- simple Gaussian MC on the observed spectrum,
10000 draws, error = std of max, bias = mean(max)-max), and N/S contrasts.  Flux scale (10%/line) cancels in N/S of the same line.
Usage: python3 -I r12_tb_err.py"""
import sys, json, numpy as np, warnings
sys.argv=['x','Xb4b']
sys.path.insert(0,'/workspace/scripts'); warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
from regions import REGIONS
C=299792.458; V0=2780.
LINES={'HCO+':('spw25',267.5576),'HCN':('spw27',265.8864)}
SC={'50pc(0.26")':0.26,'55pc(0.2835")':0.2835,'120pc(0.6186")':0.6186,'250pc(1.2887")':1.2887}; HALF=9.0; CEN=(156.9645,-43.9046)
LF=lambda v:(v<-450)|(v>600)
def P(spw): return '/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.%s.cube.I.pbcor.fits'%spw
def tbk(nu,a,b): return 1.222e3/(nu**2*a*b)
def cov(maj,mn,pa):
    s1,s2=maj/2.3548,mn/2.3548; pa=np.deg2rad(pa); u=np.array([-np.sin(pa),np.cos(pa)]); w=np.array([np.cos(pa),np.sin(pa)])
    return s1**2*np.outer(u,u)+s2**2*np.outer(w,w)
def kernel(Ct,Cn,pix):
    Ck=(Ct-Cn)/pix**2; ev=np.linalg.eigvalsh(Ck); assert ev.min()>-1e-9,ev
    sg=np.sqrt(np.sort(ev)[::-1]); R=int(np.ceil(5*sg[0]))+1; y,x=np.mgrid[-R:R+1,-R:R+1]
    Ci=np.linalg.inv(Ck+1e-9*np.eye(2)); q=Ci[0,0]*x*x+2*Ci[0,1]*x*y+Ci[1,1]*y*y; k=np.exp(-0.5*q); return k/k.sum()
def load(spw,f0):
    h=fits.open(P(spw),memmap=True); hd=h[0].header; d=h[0].data
    while d.ndim>3: d=d[0]
    nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']; v=C*(1-nu/(f0*1e9))-V0
    pix=abs(hd['CDELT2'])*3600; w=WCS(hd).celestial; cx,cy=[float(t) for t in w.world_to_pixel_values(*CEN)]; n=int(HALF/pix)
    ix,iy=int(round(cx)),int(round(cy)); y0=max(iy-n,0); x0=max(ix-n,0)
    cut=np.array(d[:,y0:iy+n+1,x0:ix+n+1],dtype=np.float32)*1e3
    bm=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); nuo=f0*(1-V0/C)
    def pxy(pos):
        a,b=w.world_to_pixel_values(*pos); return float(a)-x0,float(b)-y0
    return cut*np.float32(tbk(nuo,bm[0],bm[1])),v,pix,bm,pxy
def kern_w(v,dv,lo,hi): return np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv
DAT={ln:load(*LINES[ln]) for ln in LINES}; rng=np.random.default_rng(1); res={}; L=[]
for sn,theta in SC.items():
    sm={}
    for ln,(Tb,v,pix,bm,pxy) in DAT.items():
        k=kernel(cov(theta,theta,0),cov(*bm),pix); c=np.empty_like(Tb)
        for i in range(Tb.shape[0]): c[i]=fftconvolve(np.nan_to_num(Tb[i]),k,'same')
        c-=np.median(c[LF(v)],axis=0); sm[ln]=(c,v)
    vh=DAT['HCO+'][1]; keep=vh>-1400; cH,vH=sm['HCN']
    fi=np.interp(vh[keep],vH[::-1],np.arange(len(vH))[::-1]); i0=np.clip(np.floor(fi).astype(int),0,len(vH)-2); fr=(fi-i0).astype(np.float32)[:,None,None]
    sm['HCN']=cH[i0]*(1-fr)+cH[i0+1]*fr; sm['HCO+']=sm['HCO+'][0][keep]; v=vh[keep]; dv=abs(v[1]-v[0]); pix=DAT['HCO+'][2]; lf=LF(v)
    ny,nx=sm['HCO+'].shape[1:]; Y,X=np.mgrid[:ny,:nx]; kw=kern_w(v,dv,-350,350); m0=np.tensordot(kw,sm['HCO+'],axes=(0,0))*dv
    res[sn]={}
    for rn in ('N','S'):
        px,py=DAT['HCO+'][4](REGIONS[rn][:2]); mm=np.where(np.hypot(X-px,Y-py)*pix<0.4,m0,-np.inf); j=np.argmax(mm); iy,ix=j//nx,j%nx
        res[sn][rn]={}
        for ln in ('HCN','HCO+'):
            s3=np.convolve(sm[ln][:,iy,ix],np.ones(3)/3,'same'); sel=np.abs(v)<350
            pk=float(s3[sel].max()); sig=float(np.std(s3[lf]))
            # MC: perturb spectrum by gaussian noise with the smoothed-spectrum sigma (independent per smoothed channel is conservative-neutral at the peak)
            mc=(s3[sel][None,:]+rng.normal(0,sig,(10000,sel.sum()))).max(axis=1)
            res[sn][rn][ln]=dict(pk=pk,sig=sig,mc_std=float(mc.std()),bias=float(mc.mean()-pk),snr=pk/sig,pix=(int(ix),int(iy)))
    for ln in ('HCN','HCO+'):
        n,s=res[sn]['N'][ln],res[sn]['S'][ln]; r=n['pk']/s['pk']; e=r*np.hypot(n['sig']/n['pk'],s['sig']/s['pk']); res[sn][f'NS_{ln}']=[r,e]
    for ln in ('HCN','HCO+'):
        L.append(f'{sn:16s} {ln:5s} N {res[sn]["N"][ln]["pk"]:.3f}+-{res[sn]["N"][ln]["sig"]:.3f} K (S/N {res[sn]["N"][ln]["snr"]:.0f}, MC bias {res[sn]["N"][ln]["bias"]:+.3f})  S {res[sn]["S"][ln]["pk"]:.3f}+-{res[sn]["S"][ln]["sig"]:.3f} K (S/N {res[sn]["S"][ln]["snr"]:.0f}, bias {res[sn]["S"][ln]["bias"]:+.3f})  N/S {res[sn]["NS_"+ln][0]:.3f}+-{res[sn]["NS_"+ln][1]:.3f}')
open('r12_tb_err.txt','w').write('\n'.join(L)); print('\n'.join(L)); json.dump(res,open('r12_tb_err.json','w'),indent=1)
