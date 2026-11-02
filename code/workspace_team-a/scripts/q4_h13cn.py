# Q4a: H13CN 4-3 (345.340 GHz) at a COMMON Band-7 beam, N / S / Bridge (+ r<6" sum), with five CO(3-2)-wing models plus a new kinematic-template wing (HCO+ 4-3 profile broadened by a free Gaussian kernel, placed at CO 3-2).
# H13CN sits 400 km/s (red side) of the CO 3-2 line (345.796 GHz), whose core is out of band (spw20 ends at 342.40 GHz sky): the limit is wing-model dependent, so the full range is reported.
# HCN 4-3 (spw16) and HCO+ 4-3 (spw18) are convolved with the same kernel as H13CN (spw20). Output: out/q4/h13cn_common_beam.csv, h13cn_tau.csv, h13cn_common_beam.png
import sys; sys.path.insert(0,'/workspace/scripts')
from q4_lib import *
from scipy.optimize import least_squares, brentq
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
F13,FCO,FHCN,FHCO=345.339760,345.795990,354.505478,356.734223
TB=(4.6,3.1,87.)           # common target beam (native: 4.14x2.70 / 4.26x2.77 / 4.13x2.70 / 4.31x2.79, PA 86-88)
VLOC={'N':2778.,'S':2764.,'Bridge':2777.}
def getspec(reg):
    f20,s20,_,_=spectrum('Xb47','spw20',reg,TB); f16,s16,_,_=spectrum('Xb47','spw16',reg,TB); f18,s18,_,_=spectrum('Xb47','spw18',reg,TB)
    return f20,s20*1e3,f16,s16*1e3,f18,s18*1e3
def prep(reg,fmin=341.80,fmax=342.28):
    f20,s20,f16,s16,f18,s18=getspec(reg); v0=VLOC[reg]; sg=robust_sigma(s20[(f20>340.7)&(f20<341.7)])
    Ti=interp1d(vel(f16,FHCN)-v0,s16,bounds_error=False,fill_value=0.); Ho=interp1d(vel(f18,FHCO)-v0,s18,bounds_error=False,fill_value=0.)
    m=(f20>fmin)&(f20<fmax); return dict(f=f20[m],y=s20[m],sg=sg,v0=v0,T=Ti,H=Ho,IT=np.sum(Ti(np.linspace(-500,500,2001)))*0.5,f20=f20,s20=s20)
def models(P):
    f=P['f']; x=(f-FCO*(1-P['v0']/C))/0.1; u=(f-342.0)/0.1
    vco=(1-f/FCO)*C-P['v0']        # velocity relative to the CO 3-2 line centre at the local velocity (H13CN is at about +400 km/s)
    vv=np.linspace(-1500,1500,3001)
    def tmpl(sk):   # HCO+ 4-3 profile, same shape broadened by sk (km/s)
        h=P['H'](vv); h=gaussian_filter1d(h,sk/1.0) if sk>1 else h; return np.interp(vco,vv,h)
    return {'zero CO':(lambda p: 0*f,0,[[]]),
     'exp wing':(lambda p: np.exp(p[0]+p[1]*u),2,[[2.,.5],[0.,1.],[3.,.8]]),
     'log-quadratic wing':(lambda p: np.exp(p[0]+p[1]*u+p[2]*u**2),3,[[2.,.5,0.],[0.,1.,.05],[3,.5,.02]]),
     'log-cubic wing':(lambda p: np.exp(p[0]+p[1]*u+p[2]*u**2+p[3]*u**3),4,[[2.,.5,0.,0.],[0.,1.,.05,0.],[3,.5,.02,.01]]),
     'Lorentzian @CO':(lambda p: p[0]/(1+(x/p[1])**2),2,[[150.,3.],[300.,2.],[50.,5.]]),
     'power-law @CO':(lambda p: p[0]*np.abs(x)**(-p[1]),2,[[200.,2.],[100.,3.],[500.,1.5]]),
     'CO=broadened HCO+43 template':(lambda p: p[0]*tmpl(p[1]),2,[[2.,40.],[2.,100.],[5.,200.],[10.,60.]])}
def fit_reg(reg,fmin=341.80,fmax=342.28):
    P=prep(reg,fmin,fmax); f=P['f']; y=P['y']; sg=P['sg']; vh=(1-f/F13)*C-P['v0']; t=P['T'](vh); rows=[]; best={}
    for nm,(fn,npar,starts) in models(P).items():
        def prof(al):
            if npar==0: return np.sum(((y-al*t)/sg)**2),[]
            b=None
            for p0 in starts:
                try: rr=least_squares(lambda p:(y-al*t-fn(p))/sg,p0,max_nfev=400,bounds=((-np.inf,np.inf) if nm!='CO=broadened HCO+43 template' else ([0,5],[1e3,300])))
                except Exception: continue
                if b is None or 2*rr.cost<b[0]: b=(2*rr.cost,rr.x)
            return b if b else (np.inf,None)
        als=np.linspace(-0.05,0.6,131); chi=np.array([prof(a)[0] for a in als]); i=chi.argmin(); a0=als[i]
        ok=np.where((chi-chi[i]>=9)&(als>=a0))[0]; ul=als[ok[0]] if len(ok) else np.nan; lo1=als[chi-chi[i]<1]
        rows.append(dict(region=reg,fit_GHz=f'{fmin:.2f}-{fmax:.2f}',model=nm,npar=npar,chi2=chi[i],ndata=int(len(f)),alpha_best=a0,alpha_1s_lo=lo1.min(),alpha_1s_hi=lo1.max(),alpha_3s_upper=ul,
                         HCN_over_H13CN_flux_lower=1/ul if ul>0 else np.nan))
        best[nm]=(a0,prof(a0)[1],fn)
    return pd.DataFrame(rows),(P,best,t)
if __name__=='__main__':
    allr=[]; fig,ax=plt.subplots(1,3,figsize=(18,5))
    for a,reg in zip(ax,['N','S','Bridge']):
        for fx in (342.28,342.20):
            df,(P,best,t)=fit_reg(reg,341.80,fx); allr.append(df)
        df,(P,best,t)=fit_reg(reg)
        a.errorbar(P['f'],P['y'],P['sg'],fmt='k.')
        for nm,(al,p,fn) in best.items(): a.plot(P['f'],al*t+(fn(p) if len(p) else 0*P['f']),label=f'{nm} (a={al:.2f})')
        a.set_title(f'{reg}: H13CN 4-3 at {TB[0]}x{TB[1]}"'); a.set_xlabel('sky GHz'); a.set_ylabel('mJy/beam'); a.legend(fontsize=7)
    plt.tight_layout(); plt.savefig(OUT+'h13cn_common_beam.png',dpi=60)
    R=pd.concat(allr); R=R[~R.duplicated(['region','fit_GHz','model'])]; R.round(4).to_csv(OUT+'h13cn_common_beam.csv',index=False)
    pd.set_option('display.width',250); print(R.round(3).to_string())
