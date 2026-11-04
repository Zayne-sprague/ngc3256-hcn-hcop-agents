"""R07 (paper stage): copy of s170_q26_ratios.py with extra beams 1.0/1.5/2.2" (N,S,Bridge only; no maps) for the same-pipeline beam dependence of the N/S contrast.
Original doc: Q26: HCN/HCO+ 3-2 brightness-temperature ratio at native-matched 0.26" (50 pc) and 55/120/250 pc beams (0.2835/0.6186/1.2887", 194 pc/")
for N, S, Bridge, E1, E2, W, Sblue, S-lobe (blue), red ridge, with moment-0 and velocity-resolved (wing/core) values, plus Xb4d/Xb4f consistency.
Method: pbcor cube cutout -> RJ T_b per channel (native beam) -> flux-conserving convolution (unit-area kernel on T_b, so T_b is the beam-averaged
brightness at the new resolution; kernel = target - native covariance) -> baseline = per-pixel median of line-free channels (v_rel<-450 or >600)
-> moment maps. Ratio = T_b(HCN)/T_b(HCO+) = [I_HCN/I_HCO+]*(nu_HCO+/nu_HCN)^2 at equal beams.
Noise: per-pixel / aperture-mean std of line-free channels x dv x sqrt(sum k^2) (includes spatial correlation). Flux scale 10% per line on the ratio
(quadrature, 14%), reported separately; it cancels in region-to-region contrasts.
usage: python3 -I s170_q26_ratios.py [Xb4b|Xb4d|Xb4f]  -> /workspace/out/q26/s170_<obs>.json (+ npy maps for Xb4b)"""
import os, sys, json, numpy as np, warnings
sys.path.insert(0,'/workspace/scripts'); warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
from regions import REGIONS
C=299792.458; V0=2780.; PC=194.
LINES={'HCO+':('spw25',267.5576),'HCN':('spw27',265.8864)}
NUR=(267.5576/265.8864)**2
OBS=sys.argv[1] if len(sys.argv)>1 else 'Xb4b'
if OBS=='Xb4f': LINES={'HCO+':('spw16',267.5576),'HCN':('spw18',265.8864)}
OUT='/workspace/out/q26/'
TAG=os.environ.get('Q26TAG','')
LF=(lambda v:(v<-450)) if TAG=='_blueonly' else (lambda v:(v<-450)|(v>600))
def off(pos,dE,dN): return (pos[0]+dE/3600/np.cos(np.deg2rad(pos[1])),pos[1]+dN/3600)
S0=REGIONS['S'][:2]
REG={k:REGIONS[k][:2] for k in ('N','S','Bridge','E1','E2','W','Sblue')}
REG['Slobe']=off(S0,0.15,-1.85); REG['Ridge']=off(S0,0.15,3.3)
SPECIAL={'Slobe':(-320,-180),'Ridge':(180,320)}   # windows rel. to 2780 (as Q17)
FULL_FIXED=('N','S','Bridge')
if OBS=='Xb4b':
    SC={'0.26"(50pc)':0.26,'55pc(0.2835")':0.2835,'120pc(0.6186")':0.6186,'250pc(1.2887")':1.2887}; HALF=9.0; CEN=(156.9645,-43.9046); RG=list(REG)
elif OBS=='Xb4d':
    SC={'250pc(1.2887")':1.2887,'6.5"':6.5}; HALF=22.0; CEN=(156.9638,-43.9046); RG=['N','S','Bridge','E1','E2','Slobe','Ridge']
else:
    SC={'6.5"':6.5}; HALF=22.0; CEN=(156.9638,-43.9046); RG=['N','S','Bridge']
SC={'1.0"':1.0,'1.5"':1.5,'2.2"':2.2}; RG=['N','S','Bridge']
if OBS=='Xb4d': SC={'1.0"':1.0,'2.2"':2.2}
def P(spw): return '/data/member.uid___A001_X1273_%s.ngc_3256_sci.%s.cube.I.pbcor.fits'%(OBS,spw)
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
def main():
    DAT={ln:load(*LINES[ln]) for ln in LINES}
    import os
    res={}; vc=json.load(open(OUT+'vc_%s.json'%OBS)) if os.path.exists(OUT+'vc_%s.json'%OBS) else {}
    order=sorted(SC.items(),key=lambda kv:-kv[1]) if OBS=='Xb4b' else list(SC.items())
    order=[o for o in order if o[1]>0.5]+[o for o in order if o[1]<=0.5]
    for sn,theta in order:
        sm={}
        for ln,(Tb,v,pix,bm,pxy) in DAT.items():
            k=kernel(cov(theta,theta,0),cov(*bm),pix) if theta>max(bm[:2])*1.001 else None
            if k is None: raise SystemExit('target beam smaller than native')
            c=np.empty_like(Tb)
            for i in range(Tb.shape[0]): c[i]=fftconvolve(np.nan_to_num(Tb[i]),k,'same')
            lf=LF(v); c-=np.median(c[lf],axis=0)
            sm[ln]=(c,v)
        # the two spws have different velocity grids (HCN grid offset by ~117 km/s, channel 4.40 vs 4.38 km/s): put HCN on the HCO+ grid (linear interp.)
        vh=DAT['HCO+'][1]; keep=vh>-1400
        cH,vH=sm['HCN']; fi=np.interp(vh[keep],vH[::-1],np.arange(len(vH))[::-1]); i0=np.clip(np.floor(fi).astype(int),0,len(vH)-2); fr=(fi-i0).astype(np.float32)[:,None,None]
        sm['HCN']=cH[i0]*(1-fr)+cH[i0+1]*fr; sm['HCO+']=sm['HCO+'][0][keep]; del cH
        v=vh[keep]; dv=abs(v[1]-v[0]); pix=DAT['HCO+'][2]; lf=LF(v)
        ny,nx=sm['HCO+'].shape[1:]; Y,X=np.mgrid[:ny,:nx]
        res[sn]={}
        for rn in RG:
            pos=REG[rn]; px,py=DAT['HCO+'][4](pos)
            if rn in SPECIAL: lo,hi=SPECIAL[rn]; rad=max(theta,0.6)
            else: rad=theta
            def win(lo,hi): k=kern_w(v,dv,lo,hi); return k
            if rn in SPECIAL: kw=win(lo,hi); mpos=(px,py)
            else:
                kw=win(-350,350) if rn in FULL_FIXED or rn not in vc else win(vc[rn]-250,vc[rn]+250)
                m0=np.tensordot(kw,sm['HCO+'],axes=(0,0))*dv
                r0=np.hypot(X-px,Y-py)*pix; mm=np.where(r0<0.4,m0,-np.inf); j=np.argmax(mm); mpos=(float(X.flat[j]),float(Y.flat[j]))
            m0={ln:np.tensordot(kw,sm[ln],axes=(0,0))*dv for ln in sm}
            sn_k=np.sqrt((kw**2).sum())*dv
            sgp={ln:np.std(sm[ln][lf],axis=0)*sn_k for ln in sm}   # per-pixel sigma of m0 (K km/s)
            ipx=(int(round(mpos[1])),int(round(mpos[0])))
            ap=(np.hypot(X-mpos[0],Y-mpos[1])*pix<=rad)
            # own peaks (peak/peak) within 0.4" of nominal
            r0=np.hypot(X-px,Y-py)*pix
            own={}
            for ln in sm:
                mm=np.where(r0<0.4,m0[ln],-np.inf); j=np.argmax(mm); own[ln]=(j//nx,j%nx)
            out=dict(pos_pix=[float(mpos[0]),float(mpos[1])],rad=rad,npix_ap=int(ap.sum()),window=[float(v[kw>0][0]),float(v[kw>0][-1])])
            spec={}; sap={}
            for ln in sm:
                s=sm[ln][:,ap].mean(axis=1); spec[ln]=s; sap[ln]=float(np.std(s[lf]))
                m=m0[ln][ap].mean(); e=sap[ln]*sn_k
                out[ln]=dict(m0_pix=float(m0[ln][ipx]),m0_pix_e=float(sgp[ln][ipx]),m0_ap=float(m),m0_ap_e=float(e),
                             Tpk_own=float(m0[ln][own[ln]]/ (kw.sum()*dv)),
                             Tchan_peak=float(np.max(np.convolve(sm[ln][:,ipx[0],ipx[1]],np.ones(3)/3,'same')[np.abs(v)<350])),sig_ch_ap=sap[ln])
            a,b=out['HCN'],out['HCO+']
            def rat(x,ex,y,ey):
                r=x/y; return (float(r),float(abs(r)*np.hypot(ex/x,ey/y))) if x>0 and y>0 else (float(r),float('nan'))
            out['R_pix']=rat(a['m0_pix'],a['m0_pix_e'],b['m0_pix'],b['m0_pix_e']); out['snr_pix']=[a['m0_pix']/a['m0_pix_e'],b['m0_pix']/b['m0_pix_e']]
            out['R_ap']=rat(a['m0_ap'],a['m0_ap_e'],b['m0_ap'],b['m0_ap_e']); out['snr_ap']=[a['m0_ap']/a['m0_ap_e'],b['m0_ap']/b['m0_ap_e']]
            out['R_peakpeak']=rat(a['Tpk_own'],a['m0_pix_e']/(kw.sum()*dv),b['Tpk_own'],b['m0_pix_e']/(kw.sum()*dv))
            out['R_chanpeak']=rat(a['Tchan_peak'],a['sig_ch_ap']/np.sqrt(3),b['Tchan_peak'],b['sig_ch_ap']/np.sqrt(3))
            # velocity centroid (aperture spectrum, HCO+, |v|<350) -> vc
            s=spec['HCO+']; sel=np.abs(v)<350; w=np.clip(s[sel],0,None)
            if rn not in vc and rn not in SPECIAL and w.sum()>0: vc[rn]=float((v[sel]*w).sum()/w.sum())
            vcr=vc.get(rn,0.0) if rn not in SPECIAL else 0.0
            out['vc']=vcr
            # velocity-resolved: 50 km/s bins rel. to vc and core/wing
            bins=[]
            for lo2 in range(-350,350,50):
                k=kern_w(v,dv,vcr+lo2,vcr+lo2+50); sk=np.sqrt((k**2).sum())*dv
                A=float((k*spec['HCN']).sum()*dv); Bm=float((k*spec['HCO+']).sum()*dv); eA=sap['HCN']*sk; eB=sap['HCO+']*sk
                bins.append([lo2+25,A,eA,Bm,eB])
            out['bins']=bins
            cw={}
            for nm,(l2,h2) in {'blue':(-300,-100),'core':(-100,100),'red':(100,300)}.items():
                k=kern_w(v,dv,vcr+l2,vcr+h2); sk=np.sqrt((k**2).sum())*dv
                A=float((k*spec['HCN']).sum()*dv); Bm=float((k*spec['HCO+']).sum()*dv); eA=sap['HCN']*sk; eB=sap['HCO+']*sk
                cw[nm]=dict(HCN=A,eHCN=eA,HCO=Bm,eHCO=eB,R=A/Bm if Bm>0 and A>0 else None,eR=(A/Bm)*np.hypot(eA/A,eB/Bm) if Bm>0 and A>0 else None)
            out['cw']=cw
            res[sn][rn]=out
            print(sn,rn,'R_pix %.2f+-%.2f (snr %.0f,%.0f) R_ap %.2f+-%.2f R_pk %.2f R_chan %.2f | Tb_pk HCN %.2f HCO+ %.2f | vc %.0f'%(*out['R_pix'],*out['snr_pix'],*out['R_ap'],out['R_peakpeak'][0],out['R_chanpeak'][0],a['Tchan_peak'],b['Tchan_peak'],vcr),flush=True)
        # maps (Xb4b only, window +-350 about systemic)
        if False:
            kw=win(-350,350); sk=np.sqrt((kw**2).sum())*dv
            M={ln:np.tensordot(kw,sm[ln],axes=(0,0))*dv for ln in sm}; S={ln:np.std(sm[ln][lf],axis=0)*sk for ln in sm}
            ok=(M['HCN']>5*S['HCN'])&(M['HCO+']>5*S['HCO+']); Rm=np.where(ok,M['HCN']/M['HCO+'],np.nan)
            tag=sn.replace('"','').replace('(','_').replace(')','')
            np.savez_compressed(OUT+'maps_%s.npz'%tag,R=Rm.astype(np.float32),HCN=M['HCN'],HCO=M['HCO+'],sHCN=S['HCN'],sHCO=S['HCO+'],pix=pix)
            for nm,(l2,h2) in {'blue':(-350,-100),'red':(100,350)}.items():
                k=kern_w(v,dv,l2,h2); sk=np.sqrt((k**2).sum())*dv
                Mw={ln:np.tensordot(k,sm[ln],axes=(0,0))*dv for ln in sm}; Sw={ln:np.std(sm[ln][lf],axis=0)*sk for ln in sm}
                okw=(Mw['HCN']>3*Sw['HCN'])&(Mw['HCO+']>3*Sw['HCO+'])
                np.savez_compressed(OUT+'maps_%s_%s.npz'%(tag,nm),R=np.where(okw,Mw['HCN']/Mw['HCO+'],np.nan).astype(np.float32))
        del sm
    json.dump(res,open(OUT+'s400_%s%s.json'%(OBS,TAG),'w'),indent=1); pass
main()
