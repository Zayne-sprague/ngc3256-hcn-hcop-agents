# Q12b summary: aggregates out/q12b/runs (non-sim) -> out/q12b/q12b_summary.csv, lnZ table (prior/error sensitivity), N-vs-S Delta T table (out/q12b/dT_table.csv), Tc-vs-nc conditional medians.
import pandas as pd, numpy as np, glob, os
pd.set_option('display.width',250,'display.max_columns',60,'display.max_rows',400)
fs=[f for f in glob.glob('/workspace/out/q12b/runs/*.csv') if '_sim' not in f]
R=pd.concat([pd.read_csv(f) for f in fs],ignore_index=True)
g=R.groupby(['tag','region','grid','data','model'])
S=g.agg(chi2_prof=('chi2_prof','min'),chi2_prof_max=('chi2_prof','max'),dof=('dof','first'),p_prof=('p_prof','mean'),logZ=('logZ','mean'),logZ_sd=('logZ','std'),AIC=('AIC','min'),BIC=('BIC','min'),n=('seed','count'),
        Tc_16=('Tc_16','mean'),Tc_50=('Tc_50','mean'),Tc_84=('Tc_84','mean'),nc_50=('nc_50','mean'),**{'P(Tc>=80)':('P(Tc>=80)','mean')})
S.round(3).to_csv('/workspace/out/q12b/q12b_summary.csv')
print(S.xs(('_tbg2.73','r8'),level=('grid','data')).round(2).to_string())
# lnZ table
T=S.reset_index(); T=T[(T.data=='r8')&(T.grid=='_tbg2.73')]
piv=T.pivot_table(index=['region','tag'],columns='model',values='logZ'); 
for m in ('2cT','2cI'): piv['dlnZ_'+m]=piv[m]-piv['1z']
print(piv.round(2).to_string()); piv.round(3).to_csv('/workspace/out/q12b/lnZ_table.csv')
# posterior samples
def samp(tag,model,reg,ds,suf='tbg2.73'):
    X=[];
    for s in (1,2):
        f=f'/workspace/out/q12b/runs/{tag}_{model}_{reg}_{ds}_{suf}_s{s}.npz'
        if os.path.exists(f): z=np.load(f); X.append(10**z['X'][:,list(z['names']).index('lTc')])
    return np.concatenate(X) if X else None
rng=np.random.default_rng(0); rows=[]
for tag in ('base','us','us30','us_nc7','us_nc6','nc7','nc6','floor25','hsig50','flo2'):
  for ds in ('r8','r8a','full'):
    for model in ('1z','2cT','2cI'):
        for suf in ('tbg2.73','_tbg5'):
          if suf=='_tbg5' and tag not in('base','us'): continue
          if ds!='r8' and tag not in ('base','us'): continue
          N=samp(tag,model,'N',ds,suf); Sx=samp(tag,model,'S',ds,suf)
          if N is None or Sx is None: continue
          a=rng.choice(N,200000); b=rng.choice(Sx,200000); dT=a-b
          rows.append(dict(tag=tag,data=ds,model=model,tbg=suf,N16=np.percentile(N,16),N50=np.median(N),N84=np.percentile(N,84),S16=np.percentile(Sx,16),S50=np.median(Sx),S84=np.percentile(Sx,84),
             dT16=np.percentile(dT,16),dT50=np.median(dT),dT84=np.percentile(dT,84),P_NgtS=(dT>0).mean(),P_gt10=(dT>10).mean(),P_gt20=(dT>20).mean(),PN80=(N>=80).mean(),PS80=(Sx>=80).mean()))
D=pd.DataFrame(rows); D.round(3).to_csv('/workspace/out/q12b/dT_table.csv',index=False); print(D.round(2).to_string())
# conditional medians Tc vs nc, vs NHc (base, r8)
out=[]
for reg in 'NS':
  for tag in ('base','us'):
    for model in ('1z','2cT'):
        X=[];
        for s in (1,2):
            f=f'/workspace/out/q12b/runs/{tag}_{model}_{reg}_r8_tbg2.73_s{s}.npz'
            if os.path.exists(f): z=np.load(f); X.append(z['X']); nm=list(z['names'])
        if not X: continue
        X=np.concatenate(X); T=10**X[:,nm.index('lTc')]
        for var in ('nc','NHc'):
            v=X[:,nm.index(var)]
            edges=np.percentile(v,[0,20,40,60,80,100])
            for i in range(5):
                m=(v>=edges[i])&(v<=edges[i+1])
                out.append(dict(region=reg,tag=tag,model=model,var=var,lo=edges[i],hi=edges[i+1],T16=np.percentile(T[m],16),T50=np.median(T[m]),T84=np.percentile(T[m],84)))
        out.append(dict(region=reg,tag=tag,model=model,var='corr(lT,nc)',lo=np.corrcoef(X[:,0],X[:,nm.index('nc')])[0,1]))
        out.append(dict(region=reg,tag=tag,model=model,var='P(nc>edge-0.5)',lo=(X[:,nm.index('nc')]>X[:,nm.index('nc')].max()-0.5).mean()))
C=pd.DataFrame(out); C.round(3).to_csv('/workspace/out/q12b/Tc_conditional.csv',index=False); print(C.round(2).to_string())
