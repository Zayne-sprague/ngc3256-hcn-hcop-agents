# Q20 step d: independent lines available at the 1" (Xb4d, 1.05" circular) resolution of the ridge analysis: SO 6_5-5_4 (251.826 GHz, spw31) and HC3N 28-27 (254.699 GHz, spw29).
# Same r=1" aperture convention as s106a (mean over aperture, Jy km/s, baseline = median(-900..-450) [+450..600 if present]); noise = sd of 1" grid apertures within 10" of the S-N midpoint (excluding r<1.5" around N,S,ridge,lobe),
# plus HCN 3-2 reference at the same apertures.  3-sigma UL = max(obs,0)+3 sd.   usage: python3 -I s110d_q20_so_hc3n_1arcsec.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.signal import fftconvolve
OUT='/workspace/out/q20'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
RIDGE=offset_pos(S0,0.15,3.3); LOBE_=offset_pos(S0,0.15,-1.85); MID=offset_pos(S0,0.1,2.55)
def win(v,dv,a,b): return np.clip(np.minimum(v+dv/2,b)-np.maximum(v-dv/2,a),0,dv)/dv
for nm,ous,spw,f0 in (('HCN32','Xb4d','spw27',265.88618e9),('SO65','Xb4d','spw31',251.82577e9),('HC3N28-27','Xb4d','spw29',254.69945e9)):
    c,v,ps,wc,x0,y0,b0=load_cube(ous,spw,f0,-1500,1500,half=24.0); dv=abs(np.median(np.diff(v))); n=c.shape[1]
    c=conv(c,ps,b0,(1.05,1.05,0.0)); ba=1.1331*1.05**2/ps**2
    m=((v>-900)&(v<-450))|((v>450)&(v<600)); c=c-np.median(c[m],axis=0)
    h=(np.hypot(*np.mgrid[-int(1.2/ps):int(1.2/ps)+1,-int(1.2/ps):int(1.2/ps)+1])*ps<=1.0).astype(float)/ba
    pxy=pxy_fn(wc,x0,y0); Y,X=np.mgrid[:n,:n]; dist=lambda p:(lambda xy:np.hypot(X-xy[0],Y-xy[1])*ps)(pxy(*p))
    dR,dL,dN,dS,dM=dist(RIDGE),dist(LOBE_),dist(N0),dist(S0),dist(MID)
    step=int(round(1.0/ps)); grid=np.zeros((n,n),bool); grid[::step,::step]=True
    ok=grid&(dM<10)&(dR>1.5)&(dL>1.5)&(dN>1.5)&(dS>1.5)
    ap=lambda a,b:fftconvolve(np.tensordot(win(v,dv,a,b),c,axes=(0,0))*dv,h,mode='same')
    core=ap(-150,150); red=ap(180,320); blue=ap(-320,-180); redf=ap(400,540); bluef=ap(-540,-400)
    ir=tuple(int(round(t)) for t in pxy(*RIDGE)[::-1]); il=tuple(int(round(t)) for t in pxy(*LOBE_)[::-1]); iN=tuple(int(round(t)) for t in pxy(*N0)[::-1]); iS=tuple(int(round(t)) for t in pxy(*S0)[::-1])
    P(f'== {nm} ({ous} {spw}), 1.05" beam, r=1" apertures; channel dv {dv:.2f} km/s; control apertures N={ok.sum()}')
    P(f'   core (-150..150) mean-aperture Jy km/s: S {core[iS]:.3f} N {core[iN]:.3f} ridge {core[ir]:.3f} lobe {core[il]:.3f}')
    for wn,W,pos,ix in (('red 180-320',red,'ridge',ir),('blue -320..-180',blue,'lobe',il),('blue -320..-180',blue,'ridge (mirror)',ir),('red 180-320',red,'lobe (mirror)',il),('red 180-320',red,'N',iN),('red 180-320',red,'S',iS)):
        sd=np.std(W[ok]); mn=np.mean(W[ok]); val=W[ix]
        P(f'   {wn} @ {pos}: {val:+.4f} Jy km/s  (controls mean {mn:+.4f}, sd {sd:.4f}) -> {(val-mn)/sd:+.1f} sigma; 3-sigma UL {max(val-mn,0)+3*sd:.4f}')
        res[f'{nm}|{wn}|{pos}']=dict(val=float(val),ctrl_mean=float(mn),ctrl_sd=float(sd))
    res[f'{nm}|core']=dict(S=float(core[iS]),N=float(core[iN]),ridge=float(core[ir]),lobe=float(core[il]))
    res[f'{nm}|far']=dict(redf_ridge=float(redf[ir]),bluef_lobe=float(bluef[il]))
json.dump(res,open(OUT+'/s110d.json','w'),indent=1); open(OUT+'/s110d.txt','w').write('\n'.join(txt))
