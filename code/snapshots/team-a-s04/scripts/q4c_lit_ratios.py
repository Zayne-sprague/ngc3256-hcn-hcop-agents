# Q4c-9: further Harada+2018 (arXiv:1801.05941) Table 13 ratios at N and S (1.7" beam, K km/s) parsed from the arXiv HTML, for comparison with our SO / CH3OH / HOC+ / CS / SiO results. Output out/q4c/harada_ratios.csv
import re, pandas as pd, numpy as np
tx=open('/workspace/lit/harada/harada.txt').read(); t=tx[tx.index('Table 13: Observed peak'):][:60000].replace('−','-'); s=' '.join(t.split())
def get(pos,sp,rf):
    m=[x for x in re.finditer(r'(?:^| )(N|S|C|TNE|TSE|TSW|OS|SW) (%s) .{0,60}? %s ([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+)'%(sp,re.escape(rf)),s) if x.group(1)==pos]
    assert len(m)==1,(pos,sp,rf,len(m)); g=m[0].groups(); return float(g[4]),float(g[5])
L={'HCO+ 1-0':('HCO [+]','89.189'),'HCN 1-0':('HCN','88.634'),'CS 2-1':('CS','97.981'),'SO 2_1-3_2 (99.3)':('SO','99.3'),'SO 6_5-5_4 (219.949)':('SO','219.949'),'CH3OH 2k-1k (96.74)':('CH 3 OH','96.741'),
   'HOC+ 1-0':('HOC [+]','89.487'),'SiO 2-1':('SiO','86.847'),'HC3N 10-9':('HC 3 N','90.979'),'HNC 1-0':('HNC','90.664'),'HCO+ 3-2':('HCO [+]','267.558'),'HCN 3-2':('HCN','265.886'),'CH3OH 5k-4k (241.79)':('CH 3 OH','241.791')}
rows=[]
for pos in ('N','S'):
    I={}
    for k,(sp,rf) in L.items():
        try: I[k]=get(pos,sp,rf)
        except AssertionError: pass
    for k in I:
        if k in ('HCO+ 1-0','HCN 1-0','HCO+ 3-2','HCN 3-2'): continue
        for den in ('HCO+ 1-0','HCO+ 3-2'):
            a,ea=I[k]; b,eb=I[den]; r=a/b; rows.append(dict(region=pos,numerator=k,denominator=den,ratio=r,err=r*np.hypot(ea/a,eb/b),I_num=a,I_den=b))
    rows.append(dict(region=pos,numerator='HCN 1-0',denominator='CS 2-1',ratio=I['HCN 1-0'][0]/I['CS 2-1'][0],err=np.nan,I_num=I['HCN 1-0'][0],I_den=I['CS 2-1'][0]))
R=pd.DataFrame(rows); R.round(4).to_csv('/workspace/out/q4c/harada_ratios.csv',index=False); pd.set_option('display.width',200); print(R[R.denominator=='HCO+ 1-0'].round(3).to_string())
