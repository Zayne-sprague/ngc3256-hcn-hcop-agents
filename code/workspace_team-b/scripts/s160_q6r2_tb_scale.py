"""Q6 r2 fix 1: T_b(scale). Xb4b cubes -> T_b (RJ) maps per channel, convolved (flux-conserving, unit-area kernel on the T_b map =
beam-averaged brightness at the new resolution) to circular 55/120/250 pc beams (194 pc/") and to the Xb4d elliptical beam;
compare with native Xb4d. Core+pedestal beam-convolved model fit. python3 -I s160_q6r2_tb_scale.py"""
import sys, json, numpy as np, warnings
sys.path.insert(0,'/workspace/scripts'); warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
from regions import REGIONS
C=299792.458; V0=2780.; HW=350.; PC=194.
LINES={'HCO+3-2':('spw25',267.5576),'HCN3-2':('spw27',265.8864)}
REG={k:REGIONS[k] for k in ('N','S','Bridge','E1','E2')}
def P(obs,spw): return '/data/member.uid___A001_X1273_%s.ngc_3256_sci.%s.cube.I.pbcor.fits'%(obs,spw)
def tbk(nu,bmaj,bmin): return 1.222e3/(nu**2*bmaj*bmin)   # K per mJy/beam
def cov(maj,minr,pa_deg):  # pixel-space (x,y) covariance, FWHM->sigma ; x increases westward, PA from N through E
    s1,s2=maj/2.3548,minr/2.3548; pa=np.deg2rad(pa_deg); u=np.array([-np.sin(pa),np.cos(pa)]); w=np.array([np.cos(pa),np.sin(pa)])
    return s1**2*np.outer(u,u)+s2**2*np.outer(w,w)
def kernel(Ct,Cn,pix,n=None):
    Ck=(Ct-Cn)/pix**2; ev=np.linalg.eigvalsh(Ck); assert ev.min()>-1e-9,ev
    sg=np.sqrt(np.sort(ev)[::-1]); R=int(np.ceil(5*sg[0]))+1; y,x=np.mgrid[-R:R+1,-R:R+1]
    Ci=np.linalg.inv(Ck+1e-9*np.eye(2)); q=Ci[0,0]*x*x+2*Ci[0,1]*x*y+Ci[1,1]*y*y; k=np.exp(-0.5*q); return k/k.sum()
def boxpeak(s,v):
    sm=np.convolve(s,np.ones(3)/3,'same'); win=np.abs(v-V0)<HW; j=np.where(win)[0][np.argmax(sm[win])]; return sm[j],j
def load(obs,spw,f0,ra,de,half):
    h=fits.open(P(obs,spw),memmap=True); hd=h[0].header; d=h[0].data
    while d.ndim>3: d=d[0]
    nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']; v=C*(1-nu/(f0*1e9))
    pix=abs(hd['CDELT2'])*3600; w=WCS(hd).celestial; cx,cy=[float(t) for t in w.world_to_pixel_values(ra,de)]; n=int(half/pix)
    ix,iy=int(round(cx)),int(round(cy)); sel=np.abs(v-V0)<HW; lf=np.abs(v-V0)>600
    cut=np.array(d[:,iy-n:iy+n+1,ix-n:ix+n+1],dtype=np.float64)*1e3  # mJy/beam
    cut=cut-np.nanmedian(cut[lf],axis=0); sig=np.nanmedian(np.nanstd(cut[lf],axis=0))
    return cut,v,sel,pix,(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),(cx-(ix-n),cy-(iy-n)),sig,f0*(1-V0/C)
def peak_in(cube,v,sel,pix,c,rad):
    m0=np.nansum(cube[sel],axis=0); y,x=np.indices(m0.shape); r=np.hypot(x-c[0],y-c[1])*pix
    k=np.argmax(np.where(r<rad,m0,-np.inf)); py,px=np.unravel_index(k,m0.shape); return py,px
def main():
    out={}
    for ln,(spw,f0) in LINES.items():
        out[ln]={}
        for rn,(ra,de,_) in REG.items():
            cb,v,sel,pb,bb,cc,sigb,nuo=load('Xb4b',spw,f0,ra,de,7.0)
            Tb=cb*tbk(nuo,bb[0],bb[1])   # K per channel map
            Cn=cov(bb[0],bb[1],bb[2]); res={}
            # native
            py,px=peak_in(Tb,v,sel,pb,cc,0.4); pk,j=boxpeak(Tb[:,py,px],v)
            sg=np.std(Tb[~sel&(np.abs(v-V0)>600)][:,py,px])
            res['native 0.24"']=dict(theta=float(np.sqrt(bb[0]*bb[1])),Tb=float(pk),sig=float(sg/np.sqrt(3)))
            # coarse peak positions are the peaks of smoothed m0 within 0.4" of the nominal position
            targets={'55pc(0.28")':(0.2835,0.2835,0),'120pc(0.62")':(0.6186,0.6186,0),'250pc(1.29")':(1.2887,1.2887,0)}
            cd,vd,seld,pd_,bd,ccd,sigd,nud=load('Xb4d',spw,f0,ra,de,7.0)
            targets['Xb4d beam 0.93x0.74"']=(bd[0],bd[1],bd[2])
            for tn,(a,b,pa) in targets.items():
                k=kernel(cov(a,b,pa),Cn,pb)
                sm=np.array([fftconvolve(np.nan_to_num(Tb[i]),k,'same') if sel[i] or abs(v[i]-V0)>600 else np.zeros(Tb.shape[1:]) for i in range(len(v))])
                py,px=peak_in(sm,v,sel,pb,cc,0.4); pk,j=boxpeak(sm[:,py,px],v)
                lfm=np.abs(v-V0)>600; sg=np.std(sm[lfm][:,py,px]); res[tn]=dict(theta=float(np.sqrt(a*b)),Tb=float(pk),sig=float(sg/np.sqrt(3)),v=float(v[j]))
            # native Xb4d
            Td=cd*tbk(nud,bd[0],bd[1]); py,px=peak_in(Td,vd,seld,pd_,ccd,0.5); pk,j=boxpeak(Td[:,py,px],vd)
            sg=np.std(Td[np.abs(vd-V0)>600][:,py,px]); res['Xb4d native']=dict(theta=float(np.sqrt(bd[0]*bd[1])),Tb=float(pk),sig=float(sg/np.sqrt(3)),v=float(vd[j]))
            out[ln][rn]=res
            print(ln,rn,' | '.join('%s %.2f+-%.2f'%(k,x['Tb'],x['sig']) for k,x in res.items()),flush=True)
    json.dump(out,open('/workspace/out/q6r2/s160.json','w'),indent=1)
main()
