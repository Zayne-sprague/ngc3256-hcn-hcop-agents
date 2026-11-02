# Q12c-4: per-template-variant deblended ladders (N and S from the SAME deblend fit = draw-by-draw), from out/q1/q1_deblend_amplitudes.csv (16 'nomiss' variants).
# Ladders as in q1_deblend.py: 4-3/3-2[Xb4f] = a(4-3)/a(3-2@Xb4f); 3-2[Xb4f]/2-1 = a(3-2@Xb4f)/a(2-1). Output out/q12c/variant_ladders.csv (row = variant; cols N_/S_ x 4 data order of q12 d[0:4]) and N-S correlation.
import pandas as pd, numpy as np
A=pd.read_csv('/workspace/out/q1/q1_deblend_amplitudes.csv'); A=A[~A.variant.str.contains(r'\+miss')]
def a(line,var,c):
    r=A[(A.line==line)&(A.variant==var)&(A.comp==c)]; return r.a.iloc[0] if len(r) else np.nan
rows=[]
for v in sorted(A.variant.unique()):
    row=dict(variant=v)
    for c in 'NS':
        row[f'{c}_0']=a('HCN 3-2@Xb4f',v,c)/a('HCN 2-1',v,c)      # HCN 3-2/2-1
        row[f'{c}_1']=a('HCN 4-3',v,c)/a('HCN 3-2@Xb4f',v,c)
        row[f'{c}_2']=a('HCO+ 3-2@Xb4f',v,c)/a('HCO+ 2-1',v,c)
        row[f'{c}_3']=a('HCO+ 4-3',v,c)/a('HCO+ 3-2@Xb4f',v,c)
    rows.append(row)
R=pd.DataFrame(rows).dropna(); R.to_csv('/workspace/out/q12c/variant_ladders.csv',index=False)
print(len(R),'variants'); 
for k in range(4):
    x,y=R[f'N_{k}'],R[f'S_{k}']; print(k,'N med %.3f sd %.3f | S med %.3f sd %.3f | corr(N,S)=%.2f'%(x.median(),x.std(),y.median(),y.std(),np.corrcoef(np.log(x),np.log(y))[0,1]))
