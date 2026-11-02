# Q4c-7: SO scope check. (1) beam matching: q4_so.py convolves SO 4_5-3_4 (Xb55, native 9.4x6.3"), SO 6_5-5_4 (Xb4f, 5.7x3.5") and SO 8_8-7_7 (Xb47, 4.1x2.7") channel-wise to the single beam 9.4x6.3" PA-88 (TB in the script)
# and takes all at the field-centre position: the measurement is a central-beam (N+S+Bridge) AVERAGE, not nuclear. (2) source-size assumption: if the 6_5-5_4 (Eu 51 K) emission fills only a fraction r of the area
# over which 4_5-3_4 (Eu 24 K) emits (compact warm component), the intrinsic Tb ratio is the observed one divided by r. Trot limit (thin LTE, 3 sigma, 14% calibration) vs (f, r). Output out/q4c/so_scope.csv
import numpy as np, pandas as pd, re
src=open('/workspace/scripts/q4_so.py').read(); print('beam definition in q4_so.py:',re.search(r'^TB=.*$',src,re.M).group(0))
S=pd.read_csv('/workspace/out/q4/so_matched.csv').set_index('line'); print(S[['ous','beam','I_Kkms','eI_Kkms']])
LINES={'SO 4_5-3_4':dict(f=178.605403,Eu=24.42876,logA=-4.15357,J=5),'SO 6_5-5_4':dict(f=251.825770,Eu=50.66106,logA=-3.71555,J=5),'SO 8_8-7_7':dict(f=344.310612,Eu=87.48155,logA=-3.28519,J=8)}
def Wp(k,T): L=LINES[k]; return (2*L['J']+1)*np.exp(-L['Eu']/T)*10**L['logA']/L['f']**2
Tg=np.logspace(np.log10(3),np.log10(300),400)
def Tmax(ratio,key): R=Wp(key,Tg)/Wp('SO 4_5-3_4',Tg); return np.where(ratio<=R[0],Tg[0],np.where(ratio>=R[-1],Tg[-1],np.interp(ratio,R,Tg)))
rng=np.random.default_rng(7); N=200000; rows=[]
W45,e45=S.loc['SO 4_5-3_4','I_Kkms'],S.loc['SO 4_5-3_4','eI_Kkms']
for key in ('SO 6_5-5_4','SO 8_8-7_7'):
    W,e=S.loc[key,'I_Kkms'],S.loc[key,'eI_Kkms']
    for f in (1.0,0.7,0.5):
        for r in (1.0,0.5,0.25,0.1):
            cal=np.exp(rng.normal(0,0.14,N)); w45=np.clip(rng.normal(W45,e45,N),0.05,None); w=rng.normal(W,e,N)/(f*r)/cal
            T=Tmax(w/w45,key); rows.append(dict(line=key,f=f,area_ratio_r=r,Trot_max_95=np.percentile(T,95),Trot_max_997=np.percentile(T,99.7)))
R=pd.DataFrame(rows); R.round(1).to_csv('/workspace/out/q4c/so_scope.csv',index=False); print(R.round(1).to_string())
