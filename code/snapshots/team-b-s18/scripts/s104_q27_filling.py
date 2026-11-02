# Q27: propagate p0026 beam filling (T_b 49 -> 250 pc falls 3-8x) into the 10" RADEX ladder fits (s62_lib.GEngine, same data/calibration/priors as Q3r4 s64).
# The Q3 fits already marginalise a shared ln(eta) (flat, eta<=1, no lower bound). Here eta is varied as a PRIOR / constraint, one position at a time (N, S, Bridge):
#   fixed1   eta=1 (no dilution)            flat0   flat ln eta in [ln 1e-3, 0]  (~ s64 baseline)       f05   flat ln eta in [0.05,1] (agenda prior)
#   f01      flat ln eta in [0.01,1]         
#   floor49  eta free (<=1) AND model T_R(3-2, HCN and HCO+) >= 0.85 x peak T_b at 49 pc (p0026 s160.json, native 0.24"): all of the 10" gas as bright as the cores
#   floor250 same with the T_b at 250 pc (1.29", Xb4d-resolution-limited smooth gas): the envelope scenario
# T_R(model) = exp(M)/(1.0645*150): RADEX LVG dv=150, comparison of integrated I with peak T_b is approximate (line FWHM of the 0.24" peaks ~ 150 km/s).
# Reports Tk, n, N, eta posteriors, P(Tk>50), corr(ln eta, log n), corr(ln eta, log Tk), min chi2.  Output out/q27/s104.txt/json
import sys,json,os,numpy as np
sys.path.insert(0,'/workspace/scripts')
from scipy.stats import norm
from scipy.special import log_ndtr
import s53_ladder_fit as L, s62_lib as S
os.makedirs('/workspace/out/q27',exist_ok=True); log=[]
def P(s=''): print(s,flush=True); log.append(s)
F=S.Grid(); T=10**F.lT
TR=np.exp(F.M)/(1.0645*150)       # (models,6): HCN2,3,4,HCO2,3,4
s160=json.load(open('/workspace/out/q6r2/s160.json'))
def floors(pos,key):
    return s160['HCN3-2'][pos][key]['Tb'],s160['HCO+3-2'][pos][key]['Tb']
def lnl(E,mode,lo=None,mask=None):
    l,eh,sd,chi=E.lnlike(); eh=eh[:,0]; sd=sd[0]; chi=chi
    if mode=='fixed1': chi=chi+(eh/sd)**2; ll=-0.5*chi      # chi2 at eta=1 exactly (profile min otherwise)
    else:
        hi=0.0; ll=-0.5*chi+np.log(np.maximum(norm.cdf((hi-eh)/sd)-norm.cdf((lo-eh)/sd),1e-300))-np.log(hi-lo)   # flat-in-ln-eta prior normalised -> lnZ comparable (Occam factor)
    if mask is not None: ll=np.where(mask,ll,-np.inf)
    return ll,eh,sd,chi
def q(x,p,qs=(16,50,84)):
    o=np.argsort(x); c=np.cumsum(p[o]); return [float(x[o][min(np.searchsorted(c,a/100),len(x)-1)]) for a in qs]
def run(pos):
    I,e=L.DATA[pos]; E=S.GEngine(F,[dict(y=np.log(I),rel=e,cls='centre')],sig_b=(.1,.1,.1))
    fH,fO=floors({'N':'N','S':'S','Bridge':'Bridge'}[pos],'native 0.24"'); eH,eO=floors(pos,'250pc(1.29")')
    m49=(TR[:,1]>=0.85*fH)&(TR[:,4]>=0.85*fO); m250=(TR[:,1]>=0.85*eH)&(TR[:,4]>=0.85*eO)
    tc=F.tauH[:,2]<np.log10(6.4)
    V=[('fixed1 (eta=1)','fixed1',None,None),('flat ln eta [1e-3,1] (s64-like)','flat',np.log(1e-3),None),('flat ln eta [0.01,1]','flat',np.log(.01),None),('flat ln eta [0.05,1] (agenda prior)','flat',np.log(.05),None),
       (f'eta<=1 + T_R>=0.85 x T_b(49 pc) [{fH:.2f},{fO:.2f} K]','flat',np.log(1e-3),m49),('flat ln eta [1e-3,1] + H13CN limit tau(HCN 4-3)<6.4 (p0007, X=50)','flat',np.log(1e-3),tc),('eta<=1 + T_R>=0.85 x T_b(250 pc) + tau(HCN4-3)<6.4','flat',np.log(1e-3),m250&tc),
       (f'eta<=1 + T_R>=0.85 x T_b(250 pc) [{eH:.2f},{eO:.2f} K]','flat',np.log(1e-3),m250)]
    res={}; P(f'\n########## {pos}: 10" ladder HCN/HCO+ 2-1,3-2,4-3 (K km/s) = {np.round(I,2).tolist()}')
    for name,mode,lo,mk in V:
        ll,eh,sd,chi=lnl(E,mode,lo,mk)
        for pr in ('logT300','linn7'):
            w=S.prior_w(F,pr); lp=ll+np.log(np.where(w>0,w,1e-300)); lp=np.where(w>0,lp,-np.inf); p=np.exp(lp-lp.max()); lnZ=lp.max()+np.log(p.sum()); p/=p.sum()
            # eta posterior given model: truncated normal, sample
            rng=np.random.default_rng(1); idx=rng.choice(len(p),40000,p=p); lo_=lo if lo is not None else 0
            if mode=='fixed1': le=np.zeros(len(idx))
            else:
                a=norm.cdf((lo-eh[idx])/sd); b=norm.cdf((0-eh[idx])/sd); u=rng.random(len(idx)); le=eh[idx]+sd*norm.ppf(np.clip(a+u*(b-a),1e-300,1-1e-16))
            ln_=F.ln[idx]; lt=F.lT[idx]; ok=np.isfinite(le)
            def c(a,b): 
                return float(np.corrcoef(a[ok],b[ok])[0,1]) if ok.sum()>10 and np.std(a[ok])>0 and np.std(b[ok])>0 else float('nan')
            r=dict(lnZ=float(lnZ),chi2min=float(chi[np.isfinite(ll)].min()) if np.isfinite(ll).any() else None,
                 logn=q(F.ln,p),Tk=q(T,p),Tk5=q(T,p,(5,))[0],logN=q(F.lNO-np.log10(150),p),logA=q(F.lA,p),tauH=[10**x for x in q(F.tauH[:,1],p)],
                 PTk50=float(p[T>50].sum()),PTk100=float(p[T>100].sum()),Pn5=float(p[F.ln>5].sum()),Pn4=float(p[F.ln<4].sum()),
                 eta=[10**(x/np.log(10)) for x in np.percentile(le[ok],[16,50,84])] if ok.sum() else None,
                 corr_eta_logn=c(le,ln_),corr_eta_logTk=c(le,lt),corr_logn_logTk=float(np.corrcoef(F.ln[idx],lt)[0,1]))
            res[f'{name}|{pr}']=r
            P(f' [{pr}] {name}\n     log n {r["logn"][1]:.2f} [{r["logn"][0]:.2f},{r["logn"][2]:.2f}]  Tk {r["Tk"][1]:.0f} [{r["Tk"][0]:.0f},{r["Tk"][2]:.0f}] (5%: {r["Tk5"]:.0f})  P(Tk>50)={r["PTk50"]:.2f} P(Tk>100)={r["PTk100"]:.2f} P(n>1e5)={r["Pn5"]:.2f} P(n<1e4)={r["Pn4"]:.2f}  log N(HCO+)/dv {r["logN"][1]:.2f}  tau(HCN3-2) {r["tauH"][1]:.1f} [{r["tauH"][0]:.1f},{r["tauH"][2]:.1f}]  eta {("%.3f [%.3f,%.3f]"%(r["eta"][1],r["eta"][0],r["eta"][2])) if r["eta"] else "-"}  min chi2 {r["chi2min"]}  lnZ {r["lnZ"]:.1f}\n     corr(ln eta,log n)={r["corr_eta_logn"]:.2f} corr(ln eta,log Tk)={r["corr_eta_logTk"]:.2f} corr(log n,log Tk)={r["corr_logn_logTk"]:.2f}')
    return res
out={pos:run(pos) for pos in ('N','S','Bridge')}
json.dump(out,open('/workspace/out/q27/s104.json','w'),indent=1,default=float); open('/workspace/out/q27/s104.txt','w').write('\n'.join(log))
