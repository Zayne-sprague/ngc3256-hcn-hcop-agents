"""Q55h: tables from out/q55/joint_summary.json: lnBF per hypothesis and data set, posterior probability with flat prior over the 8 non-null hypotheses.
Output out/q55/hyp_table.csv, hyp_table.md"""
import json,numpy as np,pandas as pd
A=json.load(open('/workspace/out/q55/joint_summary.json'))
rows=[]
for k,v in A.items():
    if 'free' in k: continue
    nm=k.split('|')[0]; hy=[d for d in v if d['hyp']!='all shared']; z=np.array([d['lnZ'] for d in hy]); p=np.exp(z-z.max()); p/=p.sum()
    for d,pp in zip(hy,p): rows.append(dict(dataset=nm,hyp=d['hyp'],lnZ=round(d['lnZ'],2),lnBF_vs_best=round(d['lnBF_vs_best'],2),P_flat8=round(pp,3)))
    sh=[d for d in v if d['hyp']=='all shared'][0]; rows.append(dict(dataset=nm,hyp='all shared',lnZ=round(sh['lnZ'],2),lnBF_vs_best=round(sh['lnBF_vs_best'],2),P_flat8=0.0))
T=pd.DataFrame(rows); T.to_csv('/workspace/out/q55/hyp_table.csv',index=False)
P=T.pivot(index='hyp',columns='dataset',values='lnBF_vs_best'); order=['all shared','T differs','n differs','T,n differ (excitation)','HCN/HCO+ column ratio X differs','N(HCN) differs','N(HCN),X differ (column+abundance)','T,n,N(HCN) differ','all four differ']
P=P.loc[order]; open('/workspace/out/q55/hyp_table.md','w').write(P.to_string()); pd.set_option('display.width',250); print(P.to_string())
Q=T.pivot(index='hyp',columns='dataset',values='P_flat8').loc[order]; print(Q.to_string())
