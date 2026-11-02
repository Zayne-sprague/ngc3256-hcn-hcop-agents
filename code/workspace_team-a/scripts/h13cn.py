# H13CN 4-3 on the CO 3-2 wing (Xb47 spw20 upper edge). Template = HCN 4-3 profile of the same region (shape fixed, scale alpha = H13CN/HCN integrated ratio), plus an explicit CO-wing model with free nuisance
# parameters (exponential; log-quadratic; Lorentzian centred on CO 3-2 at the local velocity; two-Gaussian), fitted over 341.80-342.28 GHz simultaneously. Profile-likelihood (delta chi2=9) one-sided 3-sigma upper limit on alpha for every wing model;
# zero-CO case for comparison.  Output out/h13cn_limits.csv, out/h13cn_fits.png
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
C=299792.458; F13=345.339760; FCO=345.795990; FHCN=354.505478
d=np.load('/workspace/out/region_spectra.npz'); fr=d['Xb47__spw20__freq']; f16=d['Xb47__spw16__freq']
VLOC={'N':2778.,'S':2764.,'Bridge':2777.}
def run(reg,fmax=342.28,fmin=341.80):
    s=d[f'Xb47__spw20__{reg}__pk']*1e3; v0=VLOC[reg]; sg=np.std(s[(fr>340.7)&(fr<341.7)])      # mJy/beam, line-free part of the spw
    T=d[f'Xb47__spw16__{reg}__pk']*1e3; vT=(1-f16/FHCN)*C; Ti=interp1d(vT-v0,T,bounds_error=False,fill_value=0.)
    IT=np.sum(Ti(np.linspace(-500,500,2001)))*0.5   # mJy km/s
    m=(fr>fmin)&(fr<fmax); f=fr[m]; y=s[m]; vh=(1-f/F13)*C-v0; t=Ti(vh)      # template (mJy per unit alpha)
    fco=FCO*(1-v0/C); x=(f-fco)/0.1       # GHz offset from CO centre in units of 0.1 GHz (neg: higher v)
    models={'zero CO (previous)':(lambda p: 0*f,0,[]),
            'exp wing':(lambda p: np.exp(p[0]+p[1]*(f-342.0)/0.1),2,[(0,1.),(.3,1.)]),
            'log-quadratic wing':(lambda p: np.exp(p[0]+p[1]*(f-342.0)/0.1+p[2]*((f-342.0)/0.1)**2),3,[]),
            'log-cubic wing':(lambda p: np.exp(p[0]+p[1]*(f-342.0)/0.1+p[2]*((f-342.0)/0.1)**2+p[3]*((f-342.0)/0.1)**3),4,[]),
            'Lorentzian wing @CO':(lambda p: p[0]/(1+(x/ p[1])**2),2,[]),
            'power-law wing @CO':(lambda p: p[0]*np.abs(x)**(-p[1]),2,[])}
    rows=[]; best={}
    frange=f'{fmin:.2f}-{fmax:.2f}'
    for nm,(fn,npar,_) in models.items():
        def res(p,al):
            return (y-al*t-fn(p))/sg
        starts={'zero CO (previous)':[[]],'exp wing':[[2.,0.5],[0.,1.],[3.,.8]],'log-quadratic wing':[[2.,.5,0.],[0.,1.,0.05],[3,.5,.02]],'log-cubic wing':[[2.,.5,0.,0.],[0.,1.,0.05,0.],[3,.5,.02,0.01]],'Lorentzian wing @CO':[[150.,3.],[300.,2.],[50.,5.]],'power-law wing @CO':[[200.,2.],[100.,3.],[500.,1.5]]}[nm]
        def prof(al):
            b=None
            for p0 in starts:
                if npar==0: r=np.sum(res([],al)**2); return r,[]
                try: rr=least_squares(res,p0,args=(al,),max_nfev=400)
                except Exception: continue
                if b is None or 2*rr.cost<b[0]: b=(2*rr.cost,rr.x)
            return b if b else (np.inf,None)
        als=np.linspace(-0.05,0.6,131); chi=np.array([prof(a)[0] for a in als]); i=chi.argmin(); a_best=als[i]
        # upper limit: smallest alpha>=a_best with chi2 = min+9 (profile), also report best-fit and its 1-sigma
        ok=np.where((chi-chi[i]>=9)&(als>=a_best))[0]; ul=als[ok[0]] if len(ok) else np.nan
        lo1=als[np.where((chi-chi[i]<1))[0]]; 
        rows.append(dict(region=reg,fit_range_GHz=frange,model=nm,npar=npar,alpha_best=a_best,alpha_1sig_lo=lo1.min(),alpha_1sig_hi=lo1.max(),alpha_3sig_upper=ul,chi2_min=chi[i],ndata=int(m.sum()),
                         HCN_over_H13CN_lower=1/ul if ul>0 else np.nan,I_template_mJykms=IT))
        best[nm]=(a_best,prof(a_best)[1],fn)
    return pd.DataFrame(rows),(f,y,t,best,sg)
if __name__=='__main__':
    fig,ax=plt.subplots(1,3,figsize=(18,5)); allr=[]
    for a,reg in zip(ax,['N','S','Bridge']):
        df,(f,y,t,best,sg)=run(reg); allr.append(df); allr.append(run(reg,fmax=342.20)[0])
        a.errorbar(f,y,sg,fmt='k.'); 
        for nm,(al,p,fn) in best.items(): a.plot(f,al*t+(fn(p) if len(p) else 0*f),label=f'{nm} (a={al:.2f})')
        a.set_title(f'{reg}: H13CN 4-3 on CO 3-2 wing'); a.set_xlabel('sky GHz'); a.set_ylabel('mJy/beam'); a.legend(fontsize=7)
    plt.tight_layout(); plt.savefig('/workspace/out/h13cn_fits.png',dpi=60)
    R=pd.concat(allr); R.round(4).to_csv('/workspace/out/h13cn_limits.csv',index=False); pd.set_option('display.width',250); print(R.round(3).to_string())
