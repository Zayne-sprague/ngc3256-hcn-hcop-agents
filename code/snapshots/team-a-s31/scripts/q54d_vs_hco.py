"""Q54: which species differs between N and S? X/HCO+ and X/HCN for X=HCN,HNC,CN at 1-0 (2.5") and 3-2 (1.4, 2.5"). Run: python3 -I q54d_vs_hco.py (imports only functions)."""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q54_chem as c
src=open('/workspace/scripts/q54b_cn_ratios.py').read(); ns={'__name__':'lib'}
exec(src[:src.index('rows=[]')],ns); cn=ns['cn']; W=ns['W']
def get(k,bm,reg):
    if k in W: return cn(k,bm,reg)
    I,e,*_=c.meas(k,bm,reg); return I,e
rows=[]
for J,hco,sp,beams in [('1-0','HCO+1-0_00993',['HCN1-0_00993','HNC1-0_00993','CN1-0_X4d','CN1-0_main'],[2.5]),('3-2','HCO+3-2_Xb4d',['HCN3-2_Xb4d','HNC3-2_00412','CN3-2_X1df1'],[1.4,2.5])]:
  for bm in beams:
    H={r:get(hco,bm,r) for r in 'NS'}
    for k in sp:
      R={}
      for r in 'NS':
        A=get(k,bm,r); v=A[0]/H[r][0]; R[r]=(v,v*np.hypot(A[1]/A[0],H[r][1]/H[r][0]))
      D=R['N'][0]/R['S'][0]; eD=D*np.hypot(R['N'][1]/R['N'][0],R['S'][1]/R['S'][0])
      rows.append(dict(J=J,beam=bm,num=k,den=hco,RN=R['N'][0],eN=R['N'][1],RS=R['S'][0],eS=R['S'][1],NoverS=D,eD=eD)); print(rows[-1],flush=True)
pd.DataFrame(rows).round(3).to_csv('/workspace/out/q54/vs_hco.csv',index=False)
