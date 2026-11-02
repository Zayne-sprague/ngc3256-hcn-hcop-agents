# Q4 round 2, fix 3: injection-recovery for the H13CN 4-3 limit (Xb47 spw20, N nucleus beam aperture) with the SAME pipeline as s41 (template = HCN 4-3 profile, CO 3-2 wing nuisance family, profile likelihood, 3-sigma one-sided limit, noise x sqrt(chi2_red>=1)).
# (A) coverage fractions: fraction of the H13CN template flux inside each fit window and inside the band (edge 342.396 GHz); Fisher information (blue vs red of line centre); blue-half-only fit.
# (B) injection into the REAL spectrum: alpha_inj = 0.03,0.1,0.2,0.3; recovery bias = alpha_rec - (alpha_real + alpha_inj).
# (C) injection on SIMULATED spectra: truth = best-fit wing family (alpha=0 fit of the real spectrum, per family) + template*alpha_inj + noise segments drawn from the line-free 340.7-341.6 GHz spectrum (random start, random sign, reversed at random; mean-subtracted per segment);
#     recovered with all 5 wing families; report bias, scatter, 3-sigma-limit coverage (fraction of alpha=0 sims with UL<0.03... i.e. false-limit rate: UL < alpha_inj) and detection rate (best-fit > 3 sigma).
import sys,warnings,json; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.optimize import least_squares
from s41_h13cn import template,z20,nu20,F13,FCO,C,VSYS
from scipy.interpolate import interp1d
FM={ # family -> (function of (p,u,x,f), start list)
 'exp-quad':(lambda p,u,x,f:np.exp(p[0]+p[1]*u+p[2]*u**2),[[2,.5,0],[0,1,.05],[3,.5,.02]]),
 'exp-linear':(lambda p,u,x,f:np.exp(p[0]+p[1]*u),[[2,.5],[0,1],[3,.8]]),
 'power-law':(lambda p,u,x,f:p[0]*np.abs(x)**(-p[1]),[[200,2],[100,3],[500,1.5]]),
 'exp-cubic':(lambda p,u,x,f:np.exp(p[0]+p[1]*u+p[2]*u**2+p[3]*u**3),[[2,.5,0,0],[0,1,.05,0],[3,.5,.02,.01]])}
AL=np.arange(0,0.5001,0.02)
def setup(reg,fmin,fmax):
    Ti,I=template(reg,'ap'); m=(nu20>fmin)&(nu20<fmax); f=nu20[m]; vh=(1-f/F13)*C-VSYS
    vv=np.linspace(-600,600,2401); area=np.trapezoid(Ti(vv),vv); t=Ti(vh)*(I/area)
    x=(f-FCO*(1-VSYS/C))/0.1; u=(f-342.0)/0.1; return m,f,t,x,u,Ti,I,area
def prof_fit(y,t,x,u,f,sg,fam,al_grid=AL):
    fn,starts=FM[fam]; best_prev=None; chis=[]
    for al in al_grid:       # warm-started sweep in alpha
        cand=[]
        for p0 in ([best_prev] if best_prev is not None else starts):
            try: r=least_squares(lambda p:(y-al*t-fn(p,u,x,f))/sg,p0,max_nfev=300); cand.append((r.cost,r.x))
            except Exception: pass
        if best_prev is None or len(cand)==0 or (al==al_grid[0]):
            for p0 in starts[1:]:
                try: r=least_squares(lambda p:(y-al*t-fn(p,u,x,f))/sg,p0,max_nfev=300); cand.append((r.cost,r.x))
                except Exception: pass
        if not cand: chis.append(np.inf); continue
        c,p=min(cand,key=lambda z:z[0]); best_prev=p; chis.append(2*c)
    ch=np.array(chis); i=ch.argmin(); dof=len(y)-len(starts[0])-1; s2=max(ch[i]/dof,1.0); d=(ch-ch[i])/s2
    ok=np.where((d>=9)&(AL>=AL[i]))[0]; ul=AL[ok[0]] if len(ok) else np.nan
    sig=np.sqrt(max(d[0],0)) if i>0 else 0.0
    # 1-sigma interval and significance of alpha_best>0
    sel=AL[d<1]; return dict(best=AL[i],ul=ul,chi2red=ch[i]/dof,sig0=np.sqrt(max((ch[0]-ch[i])/s2,0)),lo=sel.min(),hi=sel.max())
def job(a):
    reg,win,fam_truth,fam_rec,alpha_inj,seed,mode=a
    rng=np.random.default_rng(seed); m,f,t,x,u,Ti,I,area=setup(reg,*win)
    y_all=z20[f'{reg}_ap']*1e3; lf=(nu20>340.7)&(nu20<341.6); sg0=np.std(y_all[lf]); y_real=y_all[m]
    if mode=='real': y=y_real+alpha_inj*t
    else:
        # truth wing from the alpha=0 fit of the real spectrum with family fam_truth
        fn,starts=FM[fam_truth]; best=min([least_squares(lambda p:(y_real-fn(p,u,x,f))/sg0,p0,max_nfev=500) for p0 in starts],key=lambda r:r.cost)
        wing=fn(best.x,u,x,f); seg=y_all[lf]; n=len(y_real); st=rng.integers(0,len(seg)-n); e=seg[st:st+n]-seg[st:st+n].mean()
        e=e*rng.choice([-1,1]); e=e[::-1] if rng.random()<.5 else e
        y=wing+alpha_inj*t+e
    r=prof_fit(y,t,x,u,f,sg0,fam_rec); r.update(reg=reg,win=f'{win[0]:.2f}-{win[1]:.2f}',fam_truth=fam_truth,fam_rec=fam_rec,alpha_inj=alpha_inj,mode=mode,seed=seed); return r
if __name__=='__main__':
    # ---------- (A) coverage
    print('(A) coverage of the H13CN template (N, beam aperture)')
    reg='N'; Ti,I=template(reg,'ap'); vv=np.linspace(-700,700,5601); T=Ti(vv); T=np.clip(T,0,None); cum=np.cumsum(T); cum/=cum[-1]
    vedge=(1-342.396/F13)*C-VSYS; print(f'  H13CN line centre (v=0) sky {F13*(1-VSYS/C):.3f} GHz; band edge 342.396 GHz = {vedge:+.0f} km/s from v=0')
    fr_band=np.interp(vedge,vv,cum); print(f'  fraction of the (HCN-4-3-profile) H13CN flux BEYOND the band edge (blue/high-frequency side, v<{vedge:+.0f} km/s): {fr_band:.3f}; inside the band {1-fr_band:.3f}')
    cov={'band':float(fr_band)}
    for win in ((341.80,342.20),(341.85,342.12),(341.80,342.28)):
        lo=(1-win[1]/F13)*C-VSYS; hi=(1-win[0]/F13)*C-VSYS; fw=np.interp(lo,vv,cum)  # nu high -> v low
        fw=np.interp(hi,vv,cum)-np.interp(lo,vv,cum); cov[f'{win[0]}-{win[1]}']=float(fw)
        m,f,t,x,u,_,_,_=setup(reg,*win); v=(1-f/F13)*C-VSYS; sg=np.std(z20[f'{reg}_ap'][(nu20>340.7)&(nu20<341.6)]*1e3)
        w=t**2/sg**2; blue=w[v<0].sum()/w.sum()   # v<0 = higher frequency (blue side, toward the band edge and the CO core)
        print(f'  window {win[0]}-{win[1]} GHz (v {lo:+.0f}..{hi:+.0f}): template flux fraction inside window {fw:.2f}; Fisher info sum(t^2/sigma^2): v<0 side {blue:.2f}, v>0 (red, lower-frequency) side {1-blue:.2f}')
        cov[f'fisher_blue_{win[0]}-{win[1]}']=float(blue)
    # v>0 (red-half) window 341.80-342.14: UL (the fit windows exclude v<-54 or -124 km/s because of the CO 3-2 core)
    rows=[]
    for fam in FM:
        m,f,t,x,u,_,_,_=setup('N',341.80,342.14); y=z20['N_ap'][m]*1e3; sg=np.std(z20['N_ap'][(nu20>340.7)&(nu20<341.6)]*1e3)
        r=prof_fit(y,t,x,u,f,sg,fam); r['fam']=fam; rows.append(r)
    print('  red-half-only window 341.80-342.14 GHz (v>0) fits, UL:'); print(pd.DataFrame(rows).round(3).to_string())
    cov['blue_only']=pd.DataFrame(rows).round(3).to_dict('records')
    json.dump(cov,open('/workspace/out/q4r2/s49_coverage.json','w'),indent=1,default=float)
    # ---------- (B)+(C)
    jobs=[]
    wins=[(341.80,342.20),(341.85,342.12)]
    for win in wins:
        for fam in ('exp-quad','exp-linear','power-law','exp-cubic'):
            for al in (0.03,0.1,0.2,0.3): jobs.append(('N',win,fam,fam,al,0,'real'))
            # recovered with the other family also
        for ft in ('exp-quad','exp-linear','power-law'):
            for al in (0.0,0.03,0.1,0.2,0.3):
                for fr in ('exp-quad','exp-linear','power-law','exp-cubic'):
                    for k in range(12): jobs.append(('N',win,ft,fr,al,1000*k+int(al*100)+7,'sim'))
    print(len(jobs),'jobs'); 
    with Pool(4) as p: out=p.map(job,jobs,chunksize=4)
    df=pd.DataFrame(out); df.to_csv('/workspace/out/q4r2/s49_injection.csv',index=False)
    # baseline: real spectrum alpha_best per (win,fam) with alpha_inj=0
    base=[]
    for win in wins:
        for fam in FM:
            m,f,t,x,u,_,_,_=setup('N',*win); y=z20['N_ap'][m]*1e3; sg=np.std(z20['N_ap'][(nu20>340.7)&(nu20<341.6)]*1e3); r=prof_fit(y,t,x,u,f,sg,fam); r.update(win=f'{win[0]:.2f}-{win[1]:.2f}',fam=fam); base.append(r)
    b=pd.DataFrame(base); print('\nreal spectrum (alpha_inj=0):'); print(b.round(3).to_string())
    bb={(r.win,r.fam):r.best for r in b.itertuples()}
    R=df[df['mode']=='real'].copy(); R['alpha_real']=[bb[(w,f)] for w,f in zip(R.win,R.fam_rec)]; R['bias']=R.best-(R.alpha_real+R.alpha_inj)
    print('\n(B) injection into the REAL spectrum: bias = alpha_rec - (alpha_real + alpha_inj)'); print(R.groupby(['win','alpha_inj']).bias.agg(['mean','min','max']).round(3).to_string())
    S=df[df['mode']=='sim'].copy(); S['bias']=S.best-S.alpha_inj; S['miss']=S.ul<S.alpha_inj; S['det3']=S.sig0>=3
    S['match']=S.fam_truth==S.fam_rec
    print('\n(C) simulated spectra: by window, injected alpha:  bias mean+-std (matched wing family | all recovery families), 3-sigma-UL false-limit rate (UL<alpha_inj), UL median, rate of >=3sigma "detection", chi2red median')
    g=S.groupby(['win','alpha_inj']); T=pd.DataFrame({'bias_all':g.bias.mean(),'std_all':g.bias.std(),'bias_match':S[S.match].groupby(['win','alpha_inj']).bias.mean(),'false_limit':g.miss.mean(),'UL_med':g.ul.median(),'det3':g.det3.mean(),'chi2red':g.chi2red.median(),'n':g.size()})
    print(T.round(3).to_string())
    T0=S[S.alpha_inj==0].groupby('win').apply(lambda d:pd.Series({'UL_median':d.ul.median(),'UL_p84':d.ul.quantile(.84),'UL_max':d.ul.max(),'best_median':d.best.median(),'n':len(d)})); print('\nalpha_inj=0 (pure wing+noise): UL distribution'); print(T0.round(3).to_string())
    # sensitivity: injected alpha at which 3-sigma UL > ... ; minimum tau detectable: ratio=(1-e^-tau)/(1-e^-tau/X)=1/alpha  with alpha=3*sigma_alpha
    from scipy.optimize import brentq
    print('\nSensitivity: alpha=I(H13CN)/I(HCN)=(1-e^-tau/X)/(1-e^-tau). Optically thin alpha=1/X (0.033-0.01).')
    for X in (30,50,100):
        for a in (0.075,0.105,0.12):
            fnc=lambda t:(1-np.exp(-t/X))/(1-np.exp(-t))-a
            print(f'  X={X}: alpha={a}: tau(HCN)>= {brentq(fnc,1e-3,5e3) if fnc(1e-3)<0 else float("nan"):.1f} would be needed for H13CN at the limit' )
