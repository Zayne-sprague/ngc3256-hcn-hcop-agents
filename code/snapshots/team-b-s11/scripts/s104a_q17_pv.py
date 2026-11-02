# Q17 step a: kinematics of the red ridge vs the blue lobe along the S-N axis (PA 185 through S), HCN/HCO+ 3-2.
# Xb4d @1.05" circular (flux-rescaled), pbcor, per-pixel baseline median(-900..-450). Apertures r=0.75" (mean Jy/beam over the aperture) every 0.5" along the axis.
# The wing is fit with ONE Gaussian over |v|=120..450 km/s only (red: +120..+450; blue: -450..-120) with zero baseline: the narrow Bridge/disk core (sigma<=40 km/s) is
# <1% of its peak beyond 120 km/s, except at the nuclei (two-peaked cores, +-60 km/s, plus unmodelled N/S wings) where the numbers are flagged.
# Noise per channel from line-free channels (-900..-450, +450..+650) of the same aperture spectrum.  Errors from the covariance (absolute sigma).
# Also model-independent: wing flux in sub-windows 120-200,200-280,280-360 (red) and mirrored (blue), and flux-weighted <v> in 150..400.
# usage: python3 -I s104a_q17_pv.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.optimize import curve_fit
OUT='/workspace/out/q17'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
PA=185.0; th=np.deg2rad(PA); ex=np.array([np.sin(th),np.cos(th)])
def g(v,a,c,s): return a*np.exp(-0.5*((v-c)/s)**2)
def pos_at(o): return offset_pos(S0,-o*ex[0],-o*ex[1])        # o>0 north of S along the axis
def wfl(v,dv,sp,lo,hi):
    ov=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return float((sp*ov).sum()*dv)
def wingfit(v,sp,e,side):
    k=(v>=120)&(v<=450) if side>0 else (v<=-120)&(v>=-450)
    best=None
    for c0 in (side*180,side*250,side*320):
        try:
            p,cv=curve_fit(g,v[k],sp[k],p0=[max(sp[k].max(),3*e),c0,70],sigma=np.full(k.sum(),e),absolute_sigma=True,bounds=([0,side*600,25],[1,side*600+(1200 if side>0 else -1200) if False else (600 if side>0 else -120),150]) if False else ([0,(100 if side>0 else -600),25],[1,(600 if side>0 else -100),150]),maxfev=20000)
            chi=float(np.sum(((sp[k]-g(v[k],*p))/e)**2))
            if best is None or chi<best[0]: best=(chi,p,np.sqrt(np.diag(cv)))
        except Exception: pass
    chi0=float(np.sum((sp[k]/e)**2)); 
    if best is None: return None
    chi,p,pe=best; return dict(A=p[0],c=p[1],fwhm=2.355*p[2],Ae=pe[0],ce=pe[1],fwhme=2.355*pe[2],flux=p[0]*p[2]*np.sqrt(2*np.pi),dchi2_vs0=chi0-chi,chi2=chi,n=int(k.sum()))
res={}; rng_=np.random.default_rng(11)
offs=np.arange(-4.0,7.01,0.5)
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,1.05,24.0); v=R['v']; dv=R['dv']; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0); n=R['n']
    lf=((v>-900)&(v<-450))|((v>450)&(v<650)); rows=[]; pv=[]
    for o in offs:
        m=apmask(n,R['ps'],R['pxy'](*pos_at(o)),0.75); sp=c[:,m].mean(1); e=float(np.std(sp[lf])); pv.append(sp)
        row=dict(off=float(o),rms=e)
        for side,nm in ((1,'red'),(-1,'blue')):
            f=wingfit(v,sp,e,side); row[nm]=f
            sub=[wfl(v,dv,sp,side*a,side*b) if side>0 else wfl(v,dv,sp,side*b,side*a) for a,b in ((120,200),(200,280),(280,360))]
            kk=(v*side>=150)&(v*side<=400); w=sp[kk]; row[nm+'_sub']=sub; row[nm+'_mean_v']=float((w*v[kk]).sum()/w.sum()) if w.sum()>3*e*np.sqrt(kk.sum()) else None
            row[nm+'_flux150_400']=float(w.sum()*dv)
            if row[nm+'_mean_v'] is not None:
                mc=[(lambda w2:(w2*v[kk]).sum()/w2.sum())(w+np.random.default_rng(5).normal(0,e,kk.sum()) if False else w+rng_.normal(0,e,kk.sum())) for _ in range(500)]
                row[nm+'_mean_v_err']=float(1.5*np.nanstd(mc))      # x1.5: channel noise is correlated (Hanning), white-noise MC underestimates
        rows.append(row)
        def s(f): return 'nofit' if f is None else f'A={f["A"]*1e3:.2f}+-{f["Ae"]*1e3:.2f} mJy/b c={f["c"]:+.0f}+-{f["ce"]:.0f} FWHM={f["fwhm"]:.0f}+-{f["fwhme"]:.0f} F={f["flux"]:.3f} dchi2={f["dchi2_vs0"]:.0f}'
        P(f'{mol} off {o:+.1f}" rms {e*1e3:.2f} | RED {s(row["red"])} | red sub(120-200,200-280,280-360) {np.round(row["red_sub"],3).tolist()} <v>150-400={row["red_mean_v"]} | BLUE {s(row["blue"])} blue sub {np.round(row["blue_sub"],3).tolist()} <v>={row["blue_mean_v"]}')
    res[mol]=rows; np.savez(f'{OUT}/s104a_pv_{mol}.npz',v=v,offs=offs,pv=np.array(pv))
def slope(rows,nm,lo,hi):
    r=[x for x in rows if lo<=x['off']<=hi and x[nm+'_mean_v'] is not None]; o=np.array([x['off'] for x in r]); y=np.array([x[nm+'_mean_v'] for x in r]); e=np.array([x[nm+'_mean_v_err'] for x in r])
    A=np.vstack([o-o.mean(),np.ones_like(o)]).T; W=1/e**2; C=np.linalg.inv(A.T@(A*W[:,None])); p=C@(A.T@(W*y)); chi=float(np.sum(W*(y-A@p)**2))
    return dict(slope=float(p[0]),slope_err=float(np.sqrt(C[0,0])),mean=float(p[1]),mean_err=float(np.sqrt(C[1,1])),chi2=chi,n=len(o))
for mol in list(res):
    rr=slope(res[mol],'red',1.5,4.0); bb=slope(res[mol],'blue',-3.0,-1.0); res[mol+'_slopes']=dict(red_1p5_4p0=rr,blue_m3_m1=bb)
    P(f'{mol}: <v>(150..400) vs offset: RED +1.5..+4.0" mean {rr["mean"]:+.0f}+-{rr["mean_err"]:.0f} km/s, slope {rr["slope"]:+.1f}+-{rr["slope_err"]:.1f} km/s/arcsec (chi2 {rr["chi2"]:.1f}/{rr["n"]-2}); BLUE -3.0..-1.0" mean {bb["mean"]:+.0f}+-{bb["mean_err"]:.0f}, slope {bb["slope"]:+.1f}+-{bb["slope_err"]:.1f} (chi2 {bb["chi2"]:.1f}/{bb["n"]-2}); (slope sign: per arcsec toward N)')
json.dump(res,open(OUT+'/s104a.json','w'),indent=1,default=float); open(OUT+'/s104a.txt','w').write('\n'.join(txt))
