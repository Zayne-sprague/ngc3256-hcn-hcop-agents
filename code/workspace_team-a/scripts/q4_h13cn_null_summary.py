# Q4b-4b: summary of the H13CN control / injection test (reads out/q4b/h13cn_null_*.csv): alpha at N vs control distribution, recovery bias per model, empirical significance and lower limit on HCN/H13CN.
import numpy as np, pandas as pd
O='/workspace/out/q4b/'; C=pd.read_csv(O+'h13cn_null_controls.csv'); I=pd.read_csv(O+'h13cn_null_injection.csv')
rows=[]
for (fx,m),g in C.groupby(['fx','model']):
    a=g.alpha.values; mu,sd=a.mean(),a.std(ddof=1)
    n=I[(I.region=='N')&(I.fx==fx)&(I.model==m)&(I.alpha_inj==0)].iloc[0]
    # recovery: slope and offset of alpha_rec vs alpha_inj at N, S, Bridge and controls (difference to the alpha_inj=0 fit at the same position)
    rec={}
    for ai in (0.03,0.075):
        d=[]
        for (r,p),gg in I[(I.fx==fx)&(I.model==m)].groupby(['region','pos']):
            b=gg[gg.alpha_inj==0].alpha_rec.iloc[0]; x=gg[gg.alpha_inj==ai].alpha_rec.iloc[0]; d.append(x-b)
        rec[ai]=(np.mean(d),np.std(d))
    Iflux=1.0
    rows.append(dict(fit_upper_GHz=fx,model=m,alpha_N=n.alpha_rec,err_N_stat=n.err,chi2red_N=n.chi2/(n.ndata-1-(2 if m!='log-quadratic wing' else 3)),ctrl_mean=mu,ctrl_std=sd,ctrl_min=a.min(),ctrl_max=a.max(),n_ctrl=len(a),
        z_emp=(n.alpha_rec-mu)/sd,z_emp_stat_only=n.alpha_rec/n.err if n.err>0 else np.nan,frac_ctrl_above_N=np.mean(a>=n.alpha_rec),
        recovered_for_0p03=rec[0.03][0],recovered_for_0p03_scatter=rec[0.03][1],recovered_for_0p075=rec[0.075][0],recovered_for_0p075_scatter=rec[0.075][1]))
S=pd.DataFrame(rows); 
# alpha corrected for control offset: alpha_N - ctrl_mean; 3-sigma empirical upper limit alpha_N - ctrl_mean + 3 ctrl_std (>=0)
S['alpha_N_minus_ctrl']=S.alpha_N-S.ctrl_mean; S['alpha_3sig_upper_emp']=np.clip(S.alpha_N-S.ctrl_mean,0,None)+3*S.ctrl_std
S['HCN_over_H13CN_flux_lower_emp']=1/S.alpha_3sig_upper_emp
S.round(4).to_csv(O+'h13cn_null_summary.csv',index=False)
pd.set_option('display.width',300,'display.max_columns',40); print(S.round(3).to_string())
print('\nS/Bridge alpha at injection 0:'); print(I[(I.alpha_inj==0)&I.region.isin(['S','Bridge'])].pivot_table(index=['region','fx'],columns='model',values='alpha_rec').round(3).to_string())
