# Q20 step c: forward-model fit of the 4-3 / 7-6 / SO / H13CN wing maps at the common 4.5"x3.0" beam, using the 3-2 ridge/lobe wing EXCESS (observed at 1.05", core-regression residual
# as in s106a) as spatial template.  4-3 window map = a*core + b + s_R*T_R + s_L*T_L, T = 3-2 excess tapered (Gaussian sigma 1.5") around ridge / lobe, convolved to the common beam.
# s = (4-3 flux)/(3-2 flux) of the excess at matched beam (Jy ratio; thermal Tb ratio r = s/1.777 for 4-3/3-2, s/(342.9/265.9)^2... computed in code).
# Null: same fit with templates displaced to other positions (ring 4.5-8" from the S-N midpoint, matched shifts) -> empirical sd of s.  usage: python3 -I s110c_q20_templatefit.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.ndimage import map_coordinates, shift as ndshift
OUT='/workspace/out/q20'; CB=(4.5,3.0,87.0); txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
d=np.load(OUT+'/s110b_maps.npz'); ps7=float(d['ps']); wra=d['wra']; wdec=d['wdec']; n7=wra.shape[0]
cosd=np.cos(np.deg2rad(S0[1]))
def win(v,dv,a,b): return np.clip(np.minimum(v+dv/2,b)-np.maximum(v-dv/2,a),0,dv)/dv
# --- 3-2 templates at 1.05" beam
MID=offset_pos(S0,0.1,2.55); TPL={}
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,1.05,24.0); v=R['v']; dv=R['dv']; ps=R['ps']; n=R['n']
    c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0)
    core=np.tensordot(win(v,dv,-150,150),c,axes=(0,0))*dv; red=np.tensordot(win(v,dv,180,320),c,axes=(0,0))*dv; blue=np.tensordot(win(v,dv,-320,-180),c,axes=(0,0))*dv
    Y,X=np.mgrid[:n,:n]; dist=lambda p:(lambda xy:np.hypot(X-xy[0],Y-xy[1])*ps)(R['pxy'](*p))
    dR,dL,dM=dist(RED),dist(LOBE),dist(MID)
    ok=(dM<10)&(dR>1.5)&(dL>1.5)&(core>0.1)
    out={}
    for wn,W in (('red',red),('blue',blue)):
        x=core[ok]; y=W[ok]; m=np.ones(len(x),bool)
        for _ in range(3):
            A=np.vstack([x[m],np.ones(m.sum())]).T; co=np.linalg.lstsq(A,y[m],rcond=None)[0]; r_=y-(co[0]*x+co[1]); s=np.std(r_[m]); m=np.abs(r_)<3*s
        resid=W-(co[0]*core+co[1]); out[wn]=resid
        P(f'{mol} 3-2 @1.05": {wn} = {co[0]:.4f}*core{co[1]:+.3f}, resid sd {s:.3f}')
    for sig in (1.5,2.5):
        TR=out['red']*np.exp(-dR**2/2/sig**2); TL=out['blue']*np.exp(-dL**2/2/sig**2)
        for nm,T in (('TR',TR),('TL',TL)):
            Tc=conv(T[None],ps,R['beam'],CB)[0]; px=R['pxy'](wra,wdec); TPL[(mol,sig,nm)]=map_coordinates(Tc,[px[1],px[0]],order=1,cval=0.0)
    ir=tuple(int(round(t)) for t in R['pxy'](*RED)[::-1]); il=tuple(int(round(t)) for t in R['pxy'](*LOBE)[::-1])
    P(f'  template excess at 1.05" at ridge pixel (red) {out["red"][ir]:.3f}, lobe pixel (blue) {out["blue"][il]:.3f} Jy/beam km/s')
# --- B7 grid geometry
c0,v0,ps7_,wc7,xx,yy,b07=load_cube('Xb47','spw16',354.50547e9,-100,100,half=24.0); px7=pxy_fn(wc7,xx,yy)
Y,X=np.mgrid[:n7,:n7]; mx,my=px7(*MID); rM=np.hypot(X-mx,Y-my)*ps7
FREQ={'HCN32':265.886,'HCO+32':267.558,'HCN43':354.505,'HCO+43':356.734,'CS76':342.883,'SO65':251.826,'H13CN43':345.340}
def fit(y,core,TR,TL,mask,extra=None):
    cols=[core[mask],np.ones(mask.sum()),TR[mask],TL[mask]]
    A=np.vstack(cols).T; yy_=y[mask]; msk=np.ones(len(yy_),bool)
    for _ in range(3):
        co=np.linalg.lstsq(A[msk],yy_[msk],rcond=None)[0]; r_=yy_-A@co; s=np.std(r_[msk]); msk=np.abs(r_)<3.5*s
    neff=max(6,mask.sum()*ps7**2/(1.1331*CB[0]*CB[1]))
    cov=np.linalg.inv(A[msk].T@A[msk])*s**2*(msk.sum()/neff)
    return co,np.sqrt(np.diag(cov)),s
mask=(rM<9.0)&(np.isfinite(d['HCN43|core']))
def shifted(T,dx,dy): return ndshift(T,(dy/ps7,dx/ps7),order=1,cval=0.0)
SH=[(dx,dy) for dx in np.arange(-9,9.1,1.5) for dy in np.arange(-9,9.1,1.5) if 4.5<=np.hypot(dx,dy)<=10]
sig=float(sys.argv[1]) if len(sys.argv)>1 else 1.5
def run(y,core,T1,T2,msk,which):
    """fit y=a*core+b+s1*T1+s2*T2; returns s_which, ols err, null sd (T_which displaced; shifts that keep >=70% of its power in the mask only), coefficients"""
    co,se,s=fit(y,core,T1,T2,msk); sn=[]
    Tw=T1 if which==0 else T2; p0=np.sum(Tw[msk]**2)
    for dx,dy in SH:
        Ts=shifted(Tw,dx,dy)
        if np.sum(Ts[msk]**2)<0.7*p0: continue
        try: c2,_,_=fit(y,core,Ts if which==0 else T1,T2 if which==0 else Ts,msk); sn.append(c2[2+which])
        except Exception: pass
    sn=np.array(sn); return co[2+which],se[2+which],(sn.mean(),sn.std(),len(sn)),co,s
def report(tag,val,se,nl,pk,ratio_to_T=None):
    nsd=max(se,nl[1]); return f'{tag}: s = {val:+.2f} +- {se:.2f}(OLS) / null sd {nl[1]:.2f} (mean {nl[0]:+.2f}, N={nl[2]}) -> s = {val:+.2f} +- {nsd:.2f}, flux at template peak {val*pk:+.3f} +- {nsd*pk:.3f} Jy/beam km/s'
LINES=('HCN32','HCO+32','HCN43','HCO+43','CS76','SO65','H13CN43')
for line in LINES:
    mol='HCN' if line.startswith(('HCN','CS','SO','H13')) else 'HCO+'
    core=d[f'{line}|core']; TRm,TLm=TPL[(mol,sig,'TR')],TPL[(mol,sig,'TL')]; pkR,pkL=TRm.max(),TLm.max()
    # E1: separate windows, free core coefficient
    for wn,which,T1,T2 in (('red',0,TRm,TLm),('blue',1,TRm,TLm)):
        if line=='H13CN43' and wn=='blue': continue
        y=d[f'{line}|{wn}']; msk=mask&np.isfinite(y)
        v,se,nl,co,rs=run(y,core,T1,T2,msk,which); pk=pkR if which==0 else pkL
        P(f'E1 {line:8s} {wn:4s} {"ridge" if which==0 else "lobe "}: '+report('',v,se,nl,pk)+f' | core coeff a={co[0]:+.4f}, resid sd {rs:.3f}')
        res[f'E1|{line}|{wn}']=dict(s=float(v),se=float(max(se,nl[1])),pk=float(pk),a=float(co[0]),rms=float(rs),null_mean=float(nl[0]))
    # E2: asymmetry map red-blue, templates TR and -TL
    if line!='H13CN43':
        y=d[f'{line}|red']-d[f'{line}|blue']; msk=mask&np.isfinite(y)
        for which in (0,1):
            v,se,nl,co,rs=run(y,core,TRm,-TLm,msk,which); pk=pkR if which==0 else pkL
            P(f'E2 {line:8s} red-blue {"ridge" if which==0 else "lobe "}: '+report('',v,se,nl,pk)+f' | a={co[0]:+.4f}, resid sd {rs:.3f}')
            res[f'E2|{line}|{"ridge" if which==0 else "lobe"}']=dict(s=float(v),se=float(max(se,nl[1])),pk=float(pk),a=float(co[0]),rms=float(rs))
# profile shape: per 40 km/s bin (red bins at ridge, blue bins at lobe), amplitude = coefficient x template peak, Jy/beam km/s at 4.5x3.0
for mol,l32,l43 in (('HCN','HCN32','HCN43'),('HCO+','HCO+32','HCO+43')):
    for side,which,pk in (('r',0,TPL[(mol,sig,'TR')].max()),('b',1,TPL[(mol,sig,'TL')].max())):
        rowtxt=[]
        for a in range(120,400,40):
            vals={}
            for l in (l32,l43):
                y=d[f'{l}|{side}{a}']; msk=mask&np.isfinite(y); v,se,nl,co,rs=run(y,d[f'{l}|core'],TPL[(mol,sig,'TR')],TPL[(mol,sig,'TL')],msk,which); vals[l]=(v*pk,max(se,nl[1])*pk)
            r=vals[l43][0]/vals[l32][0] if abs(vals[l32][0])>1e-3 else np.nan
            rowtxt.append(f'{a}-{a+40}: 3-2 {vals[l32][0]:+.3f}+-{vals[l32][1]:.3f}, 4-3 {vals[l43][0]:+.3f}+-{vals[l43][1]:.3f}')
            res[f'bin|{mol}|{"ridge" if which==0 else "lobe"}|{a}']=dict(f32=vals[l32],f43=vals[l43])
        P(f'bins {mol} {"ridge (red)" if which==0 else "lobe (blue)"} Jy/beam km/s @4.5x3.0: '+' | '.join(rowtxt))
json.dump(res,open(OUT+('/s110c.json' if sig==1.5 else f'/s110c_sig{sig}.json'),'w'),indent=1); open(OUT+('/s110c.txt' if sig==1.5 else f'/s110c_sig{sig}.txt'),'w').write('\n'.join(txt))
