# Q12 summary: aggregates out/q12/runs/*.csv -> out/q12/q12_summary.csv (+ printed tables). Seeds are averaged for chi2_min (min) and logZ (mean, with seed spread).
import pandas as pd, numpy as np, glob
R=pd.concat([pd.read_csv(f) for f in glob.glob('/workspace/out/q12/runs/*.csv')],ignore_index=True)
def neff(r):
    if r.model=='1z': return 4 if r.data=='r8' else 5
    base=r.npar-(1 if r.data=='r8' else 0)   # r8: fc,we enter only through their ratio
    return base
R['npar_eff']=R.apply(neff,axis=1); R['BIC_eff']=R.chi2_min+R.npar_eff*np.log(R.ndata)
g=R.groupby(['region','grid','data','model'])
keys=['Tc_16','Tc_50','Tc_84','Te_16','Te_50','Te_84','nc_50','ne_50','logXc_50','logXe_50','P(Tc>Te)','P(Tc<=40)','P(Tc>=80)','compact_frac32_HCN_50','tau32H_c_50','tau32O_c_50','lfc_50','lwe_50']
S=g.agg(chi2_min=('chi2_min','min'),chi2_seedmax=('chi2_min','max'),logZ=('logZ','mean'),logZ_sd=('logZ','std'),ndata=('ndata','first'),npar_eff=('npar_eff','first'),BIC=('BIC_eff','min'),nseeds=('seed','count')).join(g[[k for k in keys if k in R.columns]].mean())
S['dof']=S.ndata-S.npar_eff
S.round(3).to_csv('/workspace/out/q12/q12_summary.csv')
pd.set_option('display.width',250,'display.max_columns',40,'display.max_rows',200)
print(S[['chi2_min','chi2_seedmax','dof','logZ','logZ_sd','BIC','Tc_16','Tc_50','Tc_84','Te_50','nc_50','ne_50','compact_frac32_HCN_50','P(Tc>=80)']].round(2).to_string())
