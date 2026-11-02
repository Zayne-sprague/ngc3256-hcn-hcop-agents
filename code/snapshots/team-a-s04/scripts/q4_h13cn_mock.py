# Q4b-5: cross-model mock test of the H13CN fit. Truth = best-fit CO-wing of model M_true fitted to the N spectrum (with the H13CN template removed) + injected H13CN (alpha_inj of the N HCN template,
# nominal shape or width x1.25 / +30 km/s offset = shape mismatch) + Gaussian noise (rms of the N spectrum); fitted with every wing model M_fit (alpha free). 20 noise realisations, windows to 342.28 and 342.20.
# Output out/q4b/h13cn_mock.csv (bias = mean(alpha_rec - alpha_inj) and scatter per (M_true, M_fit, case)); q4_h13cn_mock_summary below.
import sys; sys.path.insert(0,'/workspace/scripts')
from q4_h13cn_null import *
from concurrent.futures import ProcessPoolExecutor
def fit_full(P,y,t,nm):
    fn,npar,starts=models(P)[nm]; sg=P['sg']; best=None
    for p0 in starts:
        try: rr=least_squares(lambda p:(y-p[0]*t-fn(p[1:]))/sg,[0.02]+list(p0),max_nfev=300,bounds=(([-1]+[-np.inf]*npar,[1]+[np.inf]*npar) if nm!='CO=broadened HCO+43 template' else ([-1,0,5],[1,1e3,300])))
        except Exception: continue
        if best is None or rr.cost<best.cost: best=rr
    return best
def job(args):
    mt,fx=args; rng=np.random.default_rng(abs(hash((mt,fx)))%2**31); P=prep('N',341.80,fx); vh=(1-P['f']/F13)*C-P['v0']; t=P['T'](vh)
    b=fit_full(P,P['y'],t,mt); fn=models(P)[mt][0]; wing=fn(b.x[1:]); sg=P['sg']
    t_wide=P['T'](vh/1.25); t_shift=P['T'](vh-30.)
    cases=[('nominal',0.,t),('nominal',0.03,t),('nominal',0.075,t),('width x1.25',0.075,t_wide),('shift +30 km/s',0.075,t_shift)]
    out=[]
    for cname,ai,tt in cases:
        recs={m:[] for m in MODELS}
        for k in range(20):
            y=wing+ai*tt+rng.normal(0,sg,len(wing))
            for m in MODELS:
                r=fit_full(P,y,t,m); recs[m].append(r.x[0] if r is not None else np.nan)
        for m in MODELS: out.append(dict(M_true=mt,fx=fx,case=cname,alpha_inj=ai,M_fit=m,bias=np.nanmean(recs[m])-ai,scatter=np.nanstd(recs[m]),rms_noise=sg))
    return out
if __name__=='__main__':
    with ProcessPoolExecutor(8) as ex: res=list(ex.map(job,[(m,f) for m in MODELS for f in (342.28,342.20)]))
    M=pd.DataFrame([r for o in res for r in o]); M.round(4).to_csv(OUTB+'h13cn_mock.csv',index=False)
    pd.set_option('display.width',250)
    print(M[(M.case=='nominal')&(M.alpha_inj==0.075)].pivot_table(index=['fx','M_true'],columns='M_fit',values='bias').round(3).to_string())
    print(M[(M.case=='nominal')&(M.alpha_inj==0.)].pivot_table(index=['fx','M_true'],columns='M_fit',values='bias').round(3).to_string())
    print(M.groupby(['case','alpha_inj']).agg(bias_mean=('bias','mean'),bias_min=('bias','min'),bias_max=('bias','max'),scatter=('scatter','median')).round(3).to_string())
