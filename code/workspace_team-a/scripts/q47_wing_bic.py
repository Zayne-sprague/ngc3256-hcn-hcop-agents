# Q47 fix 2: H13CN 4-3 alpha limit from chi2-ACCEPTABLE wing models only, alpha>=0 enforced (grid 0..0.6), BIC model selection, BIC-weighted model-averaged limit.
# Uses model definitions of h13cn.py (same data, same free-baseline wings). Windows 341.80-342.28 and -342.20 GHz. Output out/q47/wing_bic.csv, wing_bic_summary.csv
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np,pandas as pd
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
from scipy.stats import chi2 as X2
C=299792.458; F13=345.339760; FCO=345.795990; FHCN=354.505478
d=np.load('/workspace/out/region_spectra.npz'); fr=d['Xb47__spw20__freq']; f16=d['Xb47__spw16__freq']; VLOC={'N':2778.,'S':2764.,'Bridge':2777.}
AL=np.linspace(0,0.6,121)
def run(reg,fmax,fmin=341.80):
    s=d[f'Xb47__spw20__{reg}__pk']*1e3; v0=VLOC[reg]; sg=np.std(s[(fr>340.7)&(fr<341.7)])
    T=d[f'Xb47__spw16__{reg}__pk']*1e3; vT=(1-f16/FHCN)*C; Ti=interp1d(vT-v0,T,bounds_error=False,fill_value=0.)
    m=(fr>fmin)&(fr<fmax); f=fr[m]; y=s[m]; t=Ti((1-f/F13)*C-v0); fco=FCO*(1-v0/C); x=(f-fco)/0.1; u=(f-342.0)/0.1
    M={'zero CO':(lambda p:0*f,0,[[]]),'exp':(lambda p:np.exp(p[0]+p[1]*u),2,[[2.,.5],[0.,1.],[3.,.8]]),
       'log-quad':(lambda p:np.exp(p[0]+p[1]*u+p[2]*u**2),3,[[2.,.5,0.],[0.,1.,.05],[3,.5,.02]]),
       'log-cubic':(lambda p:np.exp(p[0]+p[1]*u+p[2]*u**2+p[3]*u**3),4,[[2.,.5,0.,0.],[0.,1.,.05,0.],[3,.5,.02,.01]]),
       'Lorentzian@CO':(lambda p:p[0]/(1+(x/p[1])**2),2,[[150.,3.],[300.,2.],[50.,5.]]),
       'power-law@CO':(lambda p:p[0]*np.abs(x)**(-p[1]),2,[[200.,2.],[100.,3.],[500.,1.5]])}
    out={}; n=len(f)
    for nm,(fn,npar,st) in M.items():
        ch=[]
        for al in AL:
            b=np.inf
            for p0 in st:
                if npar==0: b=np.sum(((y-al*t)/sg)**2); break
                try: r=least_squares(lambda p:(y-al*t-fn(p))/sg,p0,max_nfev=300); b=min(b,2*r.cost)
                except Exception: pass
            ch.append(b)
        out[nm]=(np.array(ch),npar,n)
    return out
rows=[];summ=[]
for reg in ('N','S','Bridge'):
  for fmax in (342.28,342.20):
    o=run(reg,fmax); n=list(o.values())[0][2]; wts={};post={}
    for nm,(ch,npar,n) in o.items():
        k=npar+1; cm=ch.min(); dof=n-k; p=X2.sf(cm,dof); bic=cm+k*np.log(n); aic=cm+2*k
        L=np.exp(-0.5*(ch-cm)); Lp=L/np.trapezoid(L,AL); cdf=np.cumsum(Lp)*(AL[1]-AL[0])
        ul=lambda q: np.interp(q,cdf/cdf[-1],AL)
        ok=np.where((ch-cm>=9)&(AL>=AL[ch.argmin()]))[0]; ul3=AL[ok[0]] if len(ok) else np.nan   # profile dchi2=9 (alpha>=0)
        rows.append(dict(region=reg,fmax=fmax,model=nm,k=k,chi2=cm,dof=dof,p_chi2=p,acc=p>0.01,AIC=aic,BIC=bic,alpha_best=AL[ch.argmin()],ul_prof_3sig=ul3,ul_post_95=ul(.95),ul_post_9987=ul(.9987)))
        post[nm]=Lp; wts[nm]=(bic,p>0.01)
    acc=[m for m in wts if wts[m][1]]; b0=min(wts[m][0] for m in acc); w=np.array([np.exp(-0.5*(wts[m][0]-b0)) for m in acc]); w/=w.sum()
    # BMA: weight by exp(-BIC/2) (Laplace evidence, flat alpha prior ~ equal), mixture of normalised alpha posteriors
    mix=sum(wi*post[m] for wi,m in zip(w,acc)); cdf=np.cumsum(mix)*(AL[1]-AL[0]); cdf/=cdf[-1]
    summ.append(dict(region=reg,fmax=fmax,accepted='+'.join(acc),weights='/'.join(f'{x:.2f}' for x in w),best_BIC=acc[int(np.argmax(w))],alpha95=np.interp(.95,cdf,AL),alpha9987=np.interp(.9987,cdf,AL),
                     HCN_over_H13CN_ge_95=1/np.interp(.95,cdf,AL),HCN_over_H13CN_ge_3sig=1/np.interp(.9987,cdf,AL)))
R=pd.DataFrame(rows); S=pd.DataFrame(summ); R.round(4).to_csv('/workspace/out/q47/wing_bic.csv',index=False); S.round(4).to_csv('/workspace/out/q47/wing_bic_summary.csv',index=False)
pd.set_option('display.width',250); print(R.round(3).to_string()); print(S.round(3).to_string())
