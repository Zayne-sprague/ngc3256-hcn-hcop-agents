# Q12c summary: (1) fit statistics with 7 independent data (chi2_prof, dof, p, AIC/BIC, lnZ, dlnZ vs one-zone); (2) H13CN survival fraction = Z(with H13CN)/Z(without) per model and its effect on Tc, n_c, tau32;
# (3) N-vs-S Delta T (independent posteriors, Gaussian-copula rho=-0.5/+0.5, and joint draw-by-draw over the 16 deblend template variants); outputs in out/q12c/.
import pandas as pd, numpy as np, glob, os
from scipy.stats import norm
pd.set_option('display.width',250,'display.max_columns',60,'display.max_rows',500)
D='/workspace/out/q12c/'
fs=[f for f in glob.glob(D+'runs/*.csv') if '_sim' not in f and os.path.basename(f).split('_')[0][:3]!='var' and not os.path.basename(f).startswith('var')]
R=pd.concat([pd.read_csv(f) for f in fs],ignore_index=True)
R['h13on']=R.h13>0
g=R.groupby(['tag','region','grid','data','model'])
S=g.agg(chi2_prof=('chi2_prof','min'),dof=('dof','first'),ndata=('ndata','first'),neff=('neff','first'),p_prof=('p_prof','mean'),logZ=('logZ','mean'),logZ_sd=('logZ','std'),AIC=('AIC','min'),BIC=('BIC','min'),n=('seed','count'),
   Tc_16=('Tc_16','mean'),Tc_50=('Tc_50','mean'),Tc_84=('Tc_84','mean'),nc_50=('nc_50','mean'),tau32_50=('tau32H_c_50','mean'),tau32_16=('tau32H_c_16','mean'),tau32_84=('tau32H_c_84','mean'),P80=('P(Tc>=80)','mean')).reset_index()
S.round(3).to_csv(D+'q12c_summary.csv',index=False)
T=S[(S.data=='r7')&(S.grid=='_tbg2.73')]
piv=T.pivot_table(index=['region','tag'],columns='model',values='logZ')
for m in ('2cT','2cI'): piv['dlnZ_'+m]=piv[m]-piv['1z']
piv.round(2).to_csv(D+'lnZ_table.csv'); print('--- lnZ (r7 = 7 independent data)'); print(piv.round(2).to_string())
print('--- one-zone fit statistics (profile chi2 incl. g,h nuisance penalties)')
print(T[T.model=='1z'][['tag','region','chi2_prof','dof','p_prof','AIC','BIC','logZ','Tc_16','Tc_50','Tc_84','nc_50','P80']].round(3).to_string())
# --- H13CN survival
rows=[]
for (reg,tag,m),a in S[(S.data=='r7')&(S.grid=='_tbg2.73')].groupby(['region','tag','model']):
    if tag.startswith('h_'): continue
    b=S[(S.region==reg)&(S.tag=='h_'+tag)&(S.model==m)&(S.data=='r7')&(S.grid=='_tbg2.73')]
    if len(b)==0: continue
    a=a.iloc[0]; b=b.iloc[0]
    rows.append(dict(region=reg,tag=tag,model=m,survival=float(np.exp(b.logZ-a.logZ)),dlnZ=b.logZ-a.logZ,Tc_noH=f'{a.Tc_50:.0f} [{a.Tc_16:.0f},{a.Tc_84:.0f}]',Tc_H=f'{b.Tc_50:.0f} [{b.Tc_16:.0f},{b.Tc_84:.0f}]',nc_noH=a.nc_50,nc_H=b.nc_50,tau32_noH=f'{a.tau32_50:.1f}',tau32_H=f'{b.tau32_50:.1f}'))
H=pd.DataFrame(rows); H.round(3).to_csv(D+'h13_survival.csv',index=False); print('--- H13CN survival (Z_with/Z_without)'); print(H.round(3).to_string())
# --- posterior samples
def samp(tag,model,reg,ds='r7',suf='tbg2.73',col='lTc'):
    X=[]
    for s in (1,2,3):
        f=D+f'runs/{tag}_{model}_{reg}_{ds}_{suf}_s{s}.npz'
        if os.path.exists(f): z=np.load(f); nm=list(z['names']); X.append(z['X'][:,nm.index(col)] if col in nm else None)
    X=[x for x in X if x is not None]; return np.concatenate(X) if X else None
def tau_stats(tag,model,reg):
    out={}
    for s in (1,2,3):
        f=D+f'runs/{tag}_{model}_{reg}_r7_tbg2.73_s{s}.npz'
        if os.path.exists(f): z=np.load(f); nm=list(z['names']); X=z['X']; out.setdefault('lT',[]).append(10**X[:,nm.index('lTc')]); out.setdefault('nc',[]).append(X[:,nm.index('nc')]); out.setdefault('NHc',[]).append(X[:,nm.index('NHc')])
    return {k:np.concatenate(v) for k,v in out.items()}
rng=np.random.default_rng(0); rows=[]
def dT_draws(N,Sx,rho,n=200000):
    z=rng.multivariate_normal([0,0],[[1,rho],[rho,1]],n); u=norm.cdf(z)
    return np.quantile(N,u[:,0])-np.quantile(Sx,u[:,1])
for tag in ('base','us','us30','geo0','geo','geo30','geo_floor25','floor25','h_base','h_geo','h_geo30','h_geo_floor25'):
    for model in ('1z','2cT'):
        N=samp(tag,model,'N'); Sx=samp(tag,model,'S')
        if N is None or Sx is None: continue
        r=dict(tag=tag,model=model,N16=np.percentile(N,16),N50=np.median(N),N84=np.percentile(N,84),S16=np.percentile(Sx,16),S50=np.median(Sx),S84=np.percentile(Sx,84))
        for rho in (-0.5,0,0.5):
            d=dT_draws(10**N,10**Sx,rho); r[f'P_NgtS_rho{rho:+.1f}']=(d>0).mean()
            if rho==0: r.update(dT16=np.percentile(d,16),dT50=np.median(d),dT84=np.percentile(d,84),P_gt10=(d>10).mean())
        rows.append(r)
DT=pd.DataFrame(rows); DT.round(3).to_csv(D+'dT_table.csv',index=False); print('--- Delta T (independent / copula)'); print(DT.round(2).to_string())
# --- joint over variants
rows=[]
for st in ('base_','geo','h_base_','h_geo'):
  for model in ('1z','2cT'):
    allN=[];allS=[];perv=[]
    for k in range(16):
        tg=f'var{k}_{st}'; N=[];Sx=[]
        for reg,L in (('N',N),('S',Sx)):
            f=D+f'runs/{tg}_{model}_{reg}_r7_tbg2.73_s1.npz'
            if os.path.exists(f): z=np.load(f); L.append(10**z['X'][:,list(z['names']).index('lTc')])
        if N and Sx:
            n=rng.choice(N[0],20000); s=rng.choice(Sx[0],20000); allN.append(n); allS.append(s); perv.append((np.median(n),np.median(s),(n>s).mean()))
    if not allN: continue
    n=np.concatenate(allN); s=np.concatenate(allS); d=n-s; pv=np.array(perv)
    rows.append(dict(setting=st,model=model,nvar=len(allN),N16=np.percentile(n,16),N50=np.median(n),N84=np.percentile(n,84),S16=np.percentile(s,16),S50=np.median(s),S84=np.percentile(s,84),dT16=np.percentile(d,16),dT50=np.median(d),dT84=np.percentile(d,84),P_NgtS=(d>0).mean(),P_NgtS_minvar=pv[:,2].min(),P_NgtS_maxvar=pv[:,2].max()))
JV=pd.DataFrame(rows); JV.round(3).to_csv(D+'dT_joint_variants.csv',index=False); print('--- joint N-S over 16 deblend variants'); print(JV.round(2).to_string())
# --- posterior shape with/without H13CN: P(tau32_c>10), logX
rows=[]
for reg in 'NS':
  for tag in ('base','geo'):
    for pre in ('','h_'):
      for model in ('1z','2cT'):
        X=[]
        for s in (1,2,3):
            f=D+f'runs/{pre}{tag}_{model}_{reg}_r7_tbg2.73_s{s}.csv'
        c=R[(R.region==reg)&(R.tag==pre+tag)&(R.model==model)&(R.data=='r7')&(R.grid=='_tbg2.73')]
        if len(c)==0: continue
        rows.append(dict(region=reg,tag=pre+tag,model=model,Tc=f"{c.Tc_50.mean():.0f} [{c.Tc_16.mean():.0f},{c.Tc_84.mean():.0f}]",nc=f"{c.nc_50.mean():.2f}",tau32=f"{c.tau32H_c_50.mean():.1f} [{c.tau32H_c_16.mean():.1f},{c.tau32H_c_84.mean():.1f}]",logXc=f"{c.logXc_50.mean():.2f}",lfc=f"{c.lfc_50.mean():.2f}"))
print('--- posterior summaries'); print(pd.DataFrame(rows).to_string())
pd.DataFrame(rows).to_csv(D+'posterior_summaries.csv',index=False)
