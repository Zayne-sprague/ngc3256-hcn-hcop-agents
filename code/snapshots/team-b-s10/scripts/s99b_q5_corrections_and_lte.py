# Q5 round 2: (1) baseline-offset (off-source median of the matched amplitude) correction of the CS 7-6 and SO 6_5-5_4 significances; (2) N/S CS ratios from the common-beam decomposition (s98) with errors;
# (3) LTE extrapolation of Harada's SO 219.949 (J_N 6_5-5_4, Eu 35 K; the line the proposal calls SO 6_5-5_4) at N to our 251.826 (J_N 5_6-4_5, Eu 50.7 K, Band 6) and 344.311 (8_8-7_7, Eu 87.5 K, Band 7) lines (A, g_u, E_u from LAMDA so@lique.dat), optically thin, same beam filling.
# Output out/q5/q5r2_summary.json    usage: python3 -I s99b_q5_corrections_and_lte.py
import numpy as np, pandas as pd, json
F=json.load(open('/workspace/out/q5/field_noise.json')); B=pd.read_csv('/workspace/out/q5/b7_fits.csv'); out={}
cs={}
for reg in ('N','S','Bridge','Sblue','E1','E2','W'):
    r=B[(B.region==reg)&(B.model=='CS1')&(B['shape']=='profile')].iloc[0]; f=F[f'Xb47|CS 7-6|{reg}']
    cs[reg]=dict(W=r.W,sig_formal=r.sig,field_mad=f['mad'],offset=f['median'],snr_raw=r.W/max(r.sig,f['mad']),W_corr=r.W-f['median'],snr_corr=(r.W-f['median'])/max(r.sig,f['mad']))
out['CS_offset_correction']=cs
so={}
for key,lab in (('Xb4f|SO 6_5-5_4|S','Xb4f S SO-only'),('Xb4f|SO 6_5-5_4|Bridge','Xb4f Bridge SO-only'),('Xb4f|SO 6_5-5_4|N','Xb4f N')):
    f=F[key]; so[key]=dict(offset=f['median'],mad=f['mad'])
rt=pd.read_csv('/workspace/out/q5/ratios.csv')
for ous,reg,note in (('Xb4f','S','SO only (CH3OH ladder ignored: upper bound) shape=profile'),('Xb4f','S','SO + CH3OH A-ladder (LTE 10 K shape, free amplitude) shape=profile'),('Xb4f','Bridge','SO only (CH3OH ladder ignored: upper bound) shape=profile'),('Xb4f','N','SO only (CH3OH ladder ignored: upper bound) shape=profile')):
    r=rt[(rt.ratio=='SO 6_5-5_4 / HCO+ 3-2')&(rt.ous==ous)&(rt.region==reg)&(rt.note==note)].iloc[0]; off=F[f'{ous}|SO 6_5-5_4|{reg}']['median']
    so[f'{ous} {reg} {note[:22]}']=dict(W=r.W_num,sig=r.sig_used,snr_raw=r.W_num/r.sig_used,W_corr=r.W_num-off,snr_corr=(r.W_num-off)/r.sig_used,R_K_corr=(r.W_num-off)/r.W_den*(267.55763/251.82577)**2)
out['SO_offset_correction']=so
# decomposition ratios
J=json.load(open('/workspace/out/q5/cs_matched_beam.json'))['common_CSbeam']['decomp']
fk=(353.4/339.77)**2   # nu_obs(HCO+4-3)/nu_obs(CS)  -> K ratio factor for Jy ratios;  HCN: (354.5*0.99077.. ) use obs freqs
nuobs={'CS 7-6':339.77,'HCN 4-3':351.2,'HCO+ 4-3':353.4}
dr={}
for model in ('N+S','N+S+Bridge'):
    for den in ('HCN 4-3','HCO+ 4-3'):
        a=J[f'CS 7-6|{model}']; b=J[f'{den}|{model}']; f=(nuobs[den]/nuobs['CS 7-6'])**2
        R={c:(a['amp'][c]/b['amp'][c]*f,a['amp'][c]/b['amp'][c]*f*np.hypot(a['err'][c]/a['amp'][c],b['err'][c]/b['amp'][c])) for c in a['amp']}
        ns=R['N'][0]/R['S'][0]; ens=ns*np.hypot(R['N'][1]/R['N'][0],R['S'][1]/R['S'][0])
        dr[f'{model} CS/{den}']=dict(R_K={c:[round(v[0],3),round(v[1],3)] for c,v in R.items()},N_over_S=[round(ns,2),round(ens,2)])
out['decomposition_ratios']=dr
# LTE
lev={'219.949':(1.335e-4,13,35.0,219.949),'251.826':(1.925e-4,11,50.7,251.826),'344.311':(5.186e-4,17,87.5,344.311)}
def w(l,T): A,g,E,nu=lev[l]; return g*A*np.exp(-E/T)/nu**2
lte={}
for T in (10,15,20,30,40,60,100):
    lte[T]={t:w(t,T)/w('219.949',T) for t in ('251.826','344.311')}
out['LTE_ratio_to_219.949']={str(k):{a:round(b,3) for a,b in v.items()} for k,v in lte.items()}
H_N=(2.42,0.74); out['Harada_N_SO219.949_K']=H_N
out['pred_SO251.826_N_Kkms']={str(T):[round(H_N[0]*lte[T]['251.826'],2),round(H_N[1]*lte[T]['251.826'],2)] for T in lte}
out['pred_SO344.311_N_Kkms']={str(T):[round(H_N[0]*lte[T]['344.311'],2),round(H_N[1]*lte[T]['344.311'],2)] for T in lte}
json.dump(out,open('/workspace/out/q5/q5r2_summary.json','w'),indent=1,default=float)
for k,v in out.items(): print(k); print(json.dumps(v,default=float)[:1800] if not isinstance(v,dict) else json.dumps(v,indent=0,default=lambda x:round(float(x),3))[:2600])
