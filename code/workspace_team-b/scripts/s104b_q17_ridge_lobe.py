# Q17 step b: red ridge vs blue lobe in MIRRORED windows (v about VSYS=2780 radio): flux, HCN/HCO+, wing/core, <v> (MC errors), 2D centroids and geometry, narrow-field velocity range.
# Data: Xb4d 3-2 HCN spw27 / HCO+ spw25 at 1.05" circular (flux rescaled), pbcor, per-pixel baseline median(-900..-450). Noise: random r=1" apertures in field (s77 field_pixels), per window.
# Positions: ridge (+0.15" E, +3.3" N of S [s77 RED] and the fitted centroid below); lobe (+0.15" E, -1.85" S of S [s77 LOBE]).
# usage: python3 -I s104b_q17_ridge_lobe.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.signal import fftconvolve
from scipy.optimize import least_squares
OUT='/workspace/out/q17'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); rng=np.random.default_rng(17)
def wflux(v,dv,sp,lo,hi):
    ov=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return float((sp*ov).sum()*dv)
res={}; DAT={}; VM={}
cosd=np.cos(np.deg2rad(S0[1]))
def rel(pos): return ((pos[0]-S0[0])*3600*cosd,(pos[1]-S0[1])*3600)     # (dE,dN) arcsec from S
def pa_deg(dE,dN): return float(np.rad2deg(np.arctan2(dE,dN))%360)
# ---- literature geometry (Michiyama+2018 Table 1)
def hms(h,m,s): return (h+m/60+s/3600)*15
M18={'N':(hms(10,27,51.23),-(43+54/60+14.00/3600)),'S':(hms(10,27,51.22),-(43+54/60+19.20/3600)),'OF-Red':(hms(10,27,51.20),-(43+54/60+15.64/3600)),
     'North-OF-Blue':(hms(10,27,51.02),-(43+54/60+13.12/3600)),'South-OF-Blue':(hms(10,27,51.17),-(43+54/60+21.43/3600))}
def drel(a,b): return ((a[0]-b[0])*3600*cosd,(a[1]-b[1])*3600)
for nm,(a,b) in {'OF-Red from M18 S':('OF-Red','S'),'South-OF-Blue from M18 S':('South-OF-Blue','S'),'OF-Red from M18 N':('OF-Red','N'),'North-OF-Blue from M18 N':('North-OF-Blue','N')}.items():
    dE,dN=drel(M18[a],M18[b]); P(f'M18 {nm}: dE={dE:+.2f}" dN={dN:+.2f}" r={np.hypot(dE,dN):.2f}" PA={pa_deg(dE,dN):.0f} deg')
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,1.05,24.0); v=R['v']; dv=R['dv']; ps=R['ps']; n=R['n']
    VM[mol]=v; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0); lf=((v>-900)&(v<-450))|((v>450)&(v<650))
    def Wm(lo,hi):
        k=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return np.tensordot(k,c,axes=(0,0))*dv
    hat1=(np.hypot(*np.mgrid[-int(1.2/ps):int(1.2/ps)+1,-int(1.2/ps):int(1.2/ps)+1])*ps<=1.0).astype(float)/R['ba']
    fm=field_pixels(R,rmin=3.5,rmax=22.0)
    def sigfield(M,rad=1.0):
        h=(np.hypot(*np.mgrid[-int(rad/ps+1):int(rad/ps+1)+1,-int(rad/ps+1):int(rad/ps+1)+1])*ps<=rad).astype(float)/R['ba']
        A=fftconvolve(M,h,mode='same'); return float(np.std(A[fm]))
    def apf(M,pos,rad=1.0):
        m=apmask(n,ps,R['pxy'](*pos),rad); return float(M[m].sum()/R['ba'])
    RED=offset_pos(S0,0.15,3.3); LOBE=offset_pos(S0,0.15,-1.85); SPOS=S0; NPOS=N0
    WIN={'180-320':(180,320),'150-350':(150,350),'240-360':(240,360)}
    for wn,(a,b) in WIN.items():
        Mr=Wm(a,b); Mb=Wm(-b,-a); sr=sigfield(Mr); sb=sigfield(Mb)
        fr=apf(Mr,RED); fb=apf(Mb,LOBE); frs=apf(Mr,offset_pos(S0,0.15,1.85)); fbN=apf(Mb,offset_pos(N0,0.15,1.85)); 
        P(f'{mol} window |v|={wn}: ridge(+{a}..{b}) at +3.3" N of S {fr:.3f}+-{sr:.3f} Jy km/s | lobe({-b}..{-a}) at -1.85" {fb:.3f}+-{sb:.3f} | red mirror position +1.85" N of S {frs:.3f} | blue at mirror-of-ridge-about-N (1.85" N of N) {fbN:.3f} | ridge/lobe = {fr/fb:.2f}+-{(fr/fb)*np.hypot(sr/fr,sb/fb):.2f}')
        res[f'{mol}|{wn}']=dict(ridge=fr,ridge_e=sr,lobe=fb,lobe_e=sb,mirror1p85=frs,blue_at_N_mirror=fbN)
    # ---- wing/core at the two positions (core |v|<150), r=1"
    Mc=Wm(-150,150); core_r=apf(Mc,RED); core_l=apf(Mc,LOBE); Mr=Wm(180,320); Mb=Wm(-320,-180)
    P(f'{mol} core(-150..150) ridge position {core_r:.2f}, lobe position {core_l:.3f} Jy km/s; wing/core ridge {apf(Mr,RED)/core_r*100:.1f}%, lobe {apf(Mb,LOBE)/core_l*100:.0f}% (wing |v|=180..320)')
    res[f'{mol}|core']=dict(ridge=core_r,lobe=core_l)
    # ---- spectra at peak positions, native-resolution not used here (1.05" beam): print spectra at r=1" as Jy/beam-mean
    for nm,pos in (('ridge',RED),('lobe',LOBE),('N',NPOS),('S',SPOS)):
        m=apmask(n,ps,R['pxy'](*pos),1.0); sp=c[:,m].mean(1); e=float(np.std(sp[lf]))
        DAT[(mol,nm)]=(sp,e)
        # flux-weighted mean velocity and sigma in the wing 150..400 with MC error (white noise, channel-independent; real noise is correlated by Hanning smoothing -> errors x1.5 fudge)
        side=1 if nm in ('ridge','N') else -1
        for sd,tag in ((1,'red'),(-1,'blue')):
            k=(v*sd>=150)&(v*sd<=400); w=sp[k]; vv=v[k]; F=w.sum()*dv
            mv=(w*vv).sum()/w.sum() if w.sum()>0 else np.nan
            mc=[]
            for _ in range(2000):
                w2=w+rng.normal(0,e,k.sum()); 
                if w2.sum()>0: mc.append((w2*vv).sum()/w2.sum())
            mc=np.array(mc); sig=np.sqrt(max((w*(vv-mv)**2).sum()/w.sum(),0)) if w.sum()>0 else np.nan
            P(f'  {mol} {nm:5s} r=1" {tag} wing 150..400: F={F:.3f} Jy/b km/s (mean over aperture, 1.05" beam), <v>={mv:+.0f} (MC sd {mc.std()*1.5:.0f}, frac MC>0 sums {len(mc)/2000:.2f}), rms-width about <v> (truncated at 150/400) {sig:.0f} km/s, noise {e*1e3:.2f} mJy/b/chan')
            res[f'{mol}|{nm}|{tag}']=dict(F=float(F),mv=float(mv),mv_err=float(mc.std()*1.5),sig=float(sig))
    # ---- 2D centroid fit of the red ridge (two beam-shaped components: free ridge + N fixed) and of the blue lobe (free) in the mirrored 180..320 windows
    def fit2d(M,pos0,fixed=None,rad=4.5,intr=True):
        xs,ys=R['pxy'](*SPOS); Y,X=np.mgrid[:n,:n]; sx=sy=0.425*1.05/ps  # beam sigma in px
        sel=np.hypot(X-xs,Y-ys)*ps<7.5
        # N component position and S component position (fixed at nuclei)
        x0,y0=R['pxy'](*pos0); nxy=R['pxy'](*NPOS)
        def model(p):
            A1,dx,dy,si,A2,A3=p
            out=A1*np.exp(-0.5*((X-(x0+dx))**2+(Y-(y0+dy))**2)/(sx**2+si**2))+A2*np.exp(-0.5*((X-nxy[0])**2+(Y-nxy[1])**2)/sx**2)+A3*np.exp(-0.5*((X-xs)**2+(Y-ys)**2)/sx**2)
            return out
        f=lambda p:((model(p)-M)[sel])/sig_
        sig_=sigfield(M)/ (np.pi*1.0**2/R['ba'])   # per-pixel-ish scale: aperture sigma / number of beams? only used as relative weight
        p0=[max(M[int(y0),int(x0)],1e-3),0,0,0.5/ps,max(M[int(nxy[1]),int(nxy[0])]*0.5,1e-3),1e-3]
        lo=[0,-3/ps,-3/ps,0,0,0]; hi=[10,3/ps,3/ps,2.5/ps,10,10]
        r=least_squares(f,p0,bounds=(lo,hi)); J=r.jac; chi2=np.sum(r.fun**2); dof=sel.sum()-6
        cov=np.linalg.pinv(J.T@J)*max(chi2/dof,1.0); e=np.sqrt(np.diag(cov))
        return r.x,e,xs,ys,nxy,ps
    for nm,M,pos0 in (('ridge180-320',Wm(180,320),RED),('lobe-320..-180',Wm(-320,-180),LOBE)):
        p,e,xs,ys,nxy,ps_=fit2d(M,pos0)
        x0,y0=R['pxy'](*pos0); cx=x0+p[1]; cy=y0+p[2]
        dE=-(cx-xs)*ps; dN=(cy-ys)*ps; sE=e[1]*ps; sN=e[2]*ps
        P(f'{mol} 2D fit {nm}: centroid dE={dE:+.2f}+-{sE:.2f}" dN={dN:+.2f}+-{sN:.2f}" from S (r={np.hypot(dE,dN):.2f}", PA={pa_deg(dE,dN):.0f} deg); intrinsic sigma {p[3]*ps:.2f}" (FWHM {2.355*p[3]*ps:.2f}" +-{2.355*e[3]*ps:.2f}); amps ridge {p[0]*1e3:.1f} N {p[4]*1e3:.1f} S {p[5]*1e3:.1f} mJy/b km/s')
        dEn=dE-rel(NPOS)[0]; dNn=dN-rel(NPOS)[1]; P(f'     relative to N: dE={dEn:+.2f}" dN={dNn:+.2f}" r={np.hypot(dEn,dNn):.2f}" PA from N {pa_deg(dEn,dNn):.0f} deg; N position rel S: dE={rel(NPOS)[0]:+.2f} dN={rel(NPOS)[1]:+.2f}')
        res[f'{mol}|fit2d|{nm}']=dict(dE=dE,dN=dN,sE=sE,sN=sN,r=float(np.hypot(dE,dN)),PA=pa_deg(dE,dN),fwhm_int=float(2.355*p[3]*ps),fwhm_err=float(2.355*e[3]*ps))
    # ---- narrow-peak velocity field near the ridge (|v|<150, Xb4d 1.05") : range of the line-peak velocity of pixels with >5 sigma core
    cm=c[np.abs(v)<150]; vm=v[np.abs(v)<150]; pk=vm[np.argmax(cm,axis=0)]; amp=cm.max(0); e_=np.std(c[lf],axis=0).mean()
    xs,ys=R['pxy'](*RED); Y,X=np.mgrid[:n,:n]; sel=(np.hypot(X-xs,Y-ys)*ps<1.5)&(amp>5*e_)
    P(f'{mol} core line-peak velocity at the ridge (r<1.5", >5 sigma, {sel.sum()} px): median {np.median(pk[sel]):+.0f}, range {pk[sel].min():+.0f}..{pk[sel].max():+.0f} km/s (cf. ridge wing <v> ~+250; s77g rotation null predicts <=0.04 vs observed 0.2 Jy km/s)')
    res[f'{mol}|core_vfield_ridge']=dict(median=float(np.median(pk[sel])),lo=float(pk[sel].min()),hi=float(pk[sel].max()))
json.dump(res,open(OUT+'/s104b.json','w'),indent=1,default=float); open(OUT+'/s104b.txt','w').write('\n'.join(txt))
# save spectra for plotting
np.savez(OUT+'/s104b_spectra.npz',v_HCN=VM['HCN'],**{'v_HCO+':VM['HCO+']},**{f'{m}_{n_}':DAT[(m,n_)][0] for (m,n_) in DAT},**{f'{m}_{n_}_e':DAT[(m,n_)][1] for (m,n_) in DAT})
