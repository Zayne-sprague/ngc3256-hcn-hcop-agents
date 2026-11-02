# Q4c: sensitivity of the class chi2 to data subset and prior (out/q4c/joint_sensitivity.csv)
import sys; sys.path.insert(0,'/workspace/scripts')
import q4c_joint as J, pandas as pd
rows=[]
for c in ('S','N'):
    for mode in ('all_no21','all','no_ladder_21_and_10'):
        for prior in ('grid','physT30'):
            r,d,u=J.run(c,mode,prior)
            for x in r: rows.append((c,mode,prior,x['cls'][:2],round(x['chi2_min'],1),x['ndata'],round(x['BIC'],1)))
df=pd.DataFrame(rows,columns=['reg','mode','prior','cls','chi2','nd','BIC']); df.to_csv('/workspace/out/q4c/joint_sensitivity.csv',index=False)
print(df.pivot_table(index=['reg','mode','prior'],columns='cls',values='chi2').to_string())
