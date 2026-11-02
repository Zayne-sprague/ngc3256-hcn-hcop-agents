# Q17 round 2, step a: core-brightness-matched null for the red ridge wing (reviewer fix 1).
# Maps (Xb4d pbcor, circular 1.05" beam, flux rescaled): core (-150..150), red wing (+180..320), mirrored blue wing (-320..-180), far-control wings (+/-400..540) and sub-windows,
# each as the mean over an r=1" aperture (Jy km/s, same convention as s104b) evaluated on a 1" grid (~independent apertures) within 10" of the S-N midpoint.
# Model: wing = a*core + b (fit on grid points excluding r<1.5" around the ridge [and lobe]); excess = observed - model; error = scatter of residuals of grid points with matched core (0.5x..2x ridge core), not field sigma.
# Baseline variants: A median(-900..-450) [as before], B linear fit to -900..-450 & +450..+650, C median of (-900..-450 & +450..+650).
# usage: python3 -I s106a_q17_core_null.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.signal import fftconvolve
OUT='/workspace/out/q17r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
cosd=np.cos(np.deg2rad(S0[1]))
def xy_rel(pos): return ((pos[0]-S0[0])*3600*cosd,(pos[1]-S0[1])*3600)
RIDGE=offset_pos(S0,0.15,3.3); LOBE_=offset_pos(S0,0.15,-1.85)
MID=offset_pos(S0,0.1,2.55)
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,1.05,24.0); v=R['v']; dv=R['dv']; ps=R['ps']; n=R['n']; raw=R['c']
    lo_=(v>-900)&(v<-450); hi_=(v>450)&(v<650)
    def base(kind):
        if kind=='A': return raw-np.median(raw[lo_],axis=0)
        if kind=='C': return raw-np.median(raw[lo_|hi_],axis=0)
        if kind=='B':
            m=lo_|hi_; X=np.vstack([np.ones(m.sum()),v[m]]).T; co=np.linalg.lstsq(X,raw[m].reshape(m.sum(),-1),rcond=None)[0]
            return raw-(np.outer(np.ones(len(v)),co[0])+np.outer(v,co[1])).reshape(raw.shape)
    h=(np.hypot(*np.mgrid[-int(1.2/ps):int(1.2/ps)+1,-int(1.2/ps):int(1.2/ps)+1])*ps<=1.0).astype(float)/R['ba']
    fm=field_pixels(R,rmin=3.5,rmax=22.0)
    Y,X=np.mgrid[:n,:n]
    mx,my=R['pxy'](*MID); sel=np.hypot(X-mx,Y-my)*ps<10.0
    # 1" grid
    step=max(1,int(round(1.0/ps))); grid=np.zeros((n,n),bool); grid[::step,::step]=True
    def dist(pos): x,y=R['pxy'](*pos); return np.hypot(X-x,Y-y)*ps
    dR,dL,dN,dS=dist(RIDGE),dist(LOBE_),dist(N0),dist(S0)
    for kind in ('A','B','C'):
        c=base(kind)
        def Wm(a,b):
            k=np.clip(np.minimum(v+dv/2,b)-np.maximum(v-dv/2,a),0,dv)/dv; return fftconvolve(np.tensordot(k,c,axes=(0,0))*dv,h,mode='same')
        core=Wm(-150,150); red=Wm(180,320); blue=Wm(-320,-180); redf=Wm(400,540); bluef=Wm(-540,-400)
        ir=tuple(int(round(t)) for t in R['pxy'](*RIDGE)[::-1]); il=tuple(int(round(t)) for t in R['pxy'](*LOBE_)[::-1])
        obsR=dict(core=core[ir],red=red[ir],blue=blue[ir]); obsL=dict(core=core[il],red=red[il],blue=blue[il])
        P(f'== {mol} baseline {kind}: ridge ap core {obsR["core"]:.2f} red180-320 {obsR["red"]:.3f} blue(mirror pos) {obsR["blue"]:.3f}; lobe ap core {obsL["core"]:.2f} blue {obsL["blue"]:.3f} red {obsL["red"]:.3f}')
        res[f'{mol}|{kind}|obs']=dict(ridge=float(red[ir]),lobe=float(blue[il]),core_r=float(core[ir]),core_l=float(core[il]))
        for excl_nm,excl in (('excl ridge+lobe (r<1.5")',(dR<1.5)|(dL<1.5)),('excl ridge+lobe+N+S cores (r<1.5")',(dR<1.5)|(dL<1.5)|(dN<1.5)|(dS<1.5))):
            ok=sel&grid&~excl&(core>0.1)
            for wn,W,ref,lab in (('red',red,obsR['red'],'ridge'),('blue',blue,obsL['blue'],'lobe'),('red@lobe-mirror',red,None,''),('redfar',redf,None,''),('bluefar',bluef,None,'')):
                if ref is None: continue
                cr=obsR['core'] if lab=='ridge' else obsL['core']
                # linear fit weighted equally; robust (Huber-like) via 2 iterations clipping 3 sigma
                x=core[ok]; y=W[ok]; msk=np.ones(len(x),bool)
                for _ in range(3):
                    A=np.vstack([x[msk],np.ones(msk.sum())]).T; co=np.linalg.lstsq(A,y[msk],rcond=None)[0]; r_=y-(co[0]*x+co[1]); s=np.std(r_[msk]); msk=np.abs(r_)<3*s
                pred=co[0]*cr+co[1]
                # matched-core scatter
                mc=(x>0.5*cr)&(x<2*cr); sm=np.std(r_[mc]) if mc.sum()>=5 else np.nan
                se=np.sqrt(np.mean(r_[mc]**2)) if mc.sum()>=5 else np.nan
                # also direct control apertures (no model): mean and sd of wing at matched core
                mean_mc=np.mean(y[mc]) if mc.sum()>=5 else np.nan; sd_mc=np.std(y[mc]) if mc.sum()>=5 else np.nan
                # leave-one-out style error on the model prediction: scatter, ignoring prediction error (fit has 2 params over many points)
                exc=ref-pred; sg=sm if np.isfinite(sm) else s
                P(f'  [{excl_nm}] {wn} wing at {lab}: obs {ref:.3f}; model wing={co[0]:.4f}*core{co[1]:+.3f} (N={ok.sum()} grid pts, resid sd {s:.3f}) -> pred {pred:.3f}; excess {exc:+.3f} +- {sg:.3f} (matched-core 0.5-2x: N={mc.sum()}, resid sd {sm:.3f}, mean wing {mean_mc:.3f} sd {sd_mc:.3f}) = {exc/sg:.1f} sigma_scatter ({ref/np.hypot(1,0)/0.025:.1f} sigma field)')
                res[f'{mol}|{kind}|{excl_nm}|{wn}']=dict(obs=float(ref),slope=float(co[0]),icpt=float(co[1]),pred=float(pred),excess=float(exc),err=float(sg),npts=int(ok.sum()),nmatch=int(mc.sum()),mean_matched=float(mean_mc),sd_matched=float(sd_mc))
        if kind=='A':
            # core-bin table, all pixels in sel (grid), excluding N/S/ridge/lobe discs
            ok=sel&grid&~((dR<1.5)|(dL<1.5)|(dN<1.5)|(dS<1.5))
            for lo,hi in ((0.3,1),(1,2),(2,4),(4,8)):
                m=ok&(core>lo)&(core<=hi)
                if m.sum()>=3: P(f'  core bin {lo}-{hi}: N={m.sum()} red {red[m].mean():+.3f}+-{red[m].std():.3f}  blue {blue[m].mean():+.3f}+-{blue[m].std():.3f}  redfar {redf[m].mean():+.3f} bluefar {bluef[m].mean():+.3f}')
            # sub-window wing fluxes at ridge / lobe, mirrored, with matched-core scatter of controls (grid pts, core 0.5-2x)
            for a_,b_ in ((150,200),(180,230),(230,280),(280,330),(330,400)):
                Wr=Wm(a_,b_); Wb=Wm(-b_,-a_); okc=sel&grid&~((dR<1.5)|(dL<1.5)|(dN<1.5)|(dS<1.5))
                mcr=okc&(core>0.5*obsR['core'])&(core<2*obsR['core']); mcl=okc&(core>0.5*obsL['core'])&(core<2*obsL['core'])
                P(f'  sub-window |v|={a_}-{b_}: ridge {Wr[ir]:.3f} (matched-core controls: mean {Wr[mcr].mean():.3f}, sd {Wr[mcr].std():.3f}, N={mcr.sum()}) | lobe {Wb[il]:.3f} (controls blue: mean {Wb[mcl].mean():.3f}, sd {Wb[mcl].std():.3f}, N={mcl.sum()}) | red@lobe-core controls: {Wr[mcl].mean():.3f}+-{Wr[mcl].std():.3f}')
                res[f'{mol}|sub|{a_}-{b_}']=dict(ridge=float(Wr[ir]),ridge_ctrl_mean=float(Wr[mcr].mean()),ridge_ctrl_sd=float(Wr[mcr].std()),lobe=float(Wb[il]),lobe_ctrl_mean=float(Wb[mcl].mean()),lobe_ctrl_sd=float(Wb[mcl].std()))
            np.savez(f'{OUT}/s106a_maps_{mol}.npz',core=core,red=red,blue=blue,redfar=redf,bluefar=bluef,ps=ps,sel=sel,grid=grid,ir=ir,il=il,dR=dR,dL=dL,dN=dN,dS=dS)
json.dump(res,open(OUT+'/s106a.json','w'),indent=1); open(OUT+'/s106a.txt','w').write('\n'.join(txt))
