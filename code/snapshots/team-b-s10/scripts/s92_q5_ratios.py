# Q5: ratio tables (CS 7-6/HCN 4-3, SO/HCO+, SO/CS, ...) in Jy/beam km/s and in brightness-temperature (K km/s) units, with sigma = max(formal, off-source field MAD);
# detections >=5 sigma, tentative 3-5 sigma, otherwise 3 and 5 sigma limits (k*sigma_emp, measured from zero). Uses out/q5/{b7,b6}_fits.csv, field_noise.json.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, json
b7=pd.read_csv('/workspace/out/q5/b7_fits.csv'); b6=pd.read_csv('/workspace/out/q5/b6_fits.csv'); fn=json.load(open('/workspace/out/q5/field_noise.json'))
CAL=0.05   # relative flux-scale uncertainty between spws of the same OUS (added to the ratio)
REG=['N','S','Bridge','E1','E2','Sblue','W']
def get(df,reg,model,line,shape='profile',ous=None):
    q=df[(df.region==reg)&(df.model==model)&(df.line==line)&(df['shape']==shape)]
    if ous is not None: q=q[q.ous==ous]
    return None if q.empty else q.iloc[0]
def status(W,s):
    z=W/s; return 'detected' if z>=5 else ('tentative' if z>=3 else 'not detected')
rows=[]
def entry(tag,ous,reg,num_line,den_line,Wn,sn_formal,fkey,Wd,sd,nu_n,nu_d,extra=''):
    sf=fn.get(fkey,{}).get('mad',np.nan); s=max(sn_formal,sf) if np.isfinite(sf) else sn_formal
    z=Wn/s; R=Wn/Wd; Rk=R*(nu_d/nu_n)**2
    dR=R*np.sqrt((s/Wd)**2+(abs(Wn)*sd/Wd**2/R if Wn!=0 else 0)**2+CAL**2) if Wn!=0 else s/Wd
    rows.append(dict(ratio=tag,ous=ous,region=reg,W_num=Wn,sig_formal=sn_formal,sig_field=sf,sig_used=s,snr=z,status=status(Wn,s),W_den=Wd,
        R_Jy=R,R_Jy_err=(dR if Wn!=0 else np.nan),R_K=Rk,R_K_err=dR*(nu_d/nu_n)**2,UL3_Jy=3*s,UL5_Jy=5*s,R_K_UL3=3*s/Wd*(nu_d/nu_n)**2,R_K_UL5=5*s/Wd*(nu_d/nu_n)**2,note=extra))
for sh in ('profile','gauss'):
  for reg in REG:
    hcn=get(b7,reg,'ref','HCN 4-3',sh); hco=get(b7,reg,'ref','HCO+ 4-3',sh); cs1=get(b7,reg,'CS1','CS 7-6',sh); cs3=get(b7,reg,'CS3','CS 7-6',sh); so=get(b7,reg,'SO1','SO 8_8-7_7',sh)
    for nm,c in (('CS1 (all blend as CS)',cs1),):
        entry('CS 7-6 / HCN 4-3','Xb47',reg,'CS 7-6','HCN 4-3',c.W,c.sig,f'Xb47|CS 7-6|{reg}',hcn.W,hcn.sig,342.88285,354.50547,nm+f' shape={sh}')
        entry('CS 7-6 / HCO+ 4-3','Xb47',reg,'CS 7-6','HCO+ 4-3',c.W,c.sig,f'Xb47|CS 7-6|{reg}',hco.W,hco.sig,342.88285,356.73424,nm+f' shape={sh}')
    entry('SO 8_8-7_7 / HCO+ 4-3','Xb47',reg,'SO','HCO+ 4-3',so.W,so.sig,f'Xb47|SO 8_8-7_7|{reg}',hco.W,hco.sig,344.31061,356.73424,f'shape={sh}')
    entry('SO 8_8-7_7 / CS 7-6','Xb47',reg,'SO','CS 7-6',so.W,so.sig,f'Xb47|SO 8_8-7_7|{reg}',cs1.W,cs1.sig,344.31061,342.88285,f'shape={sh}')
    for ous in ('Xb4f','Xb4d','Xb4b'):
        hco3=get(b6,reg,'ref','HCO+ 3-2',sh,ous); hcn3=get(b6,reg,'ref','HCN 3-2',sh,ous)
        if hco3 is None: continue
        for model,lab in (('SO_only','SO only (CH3OH ladder ignored: upper bound)'),('SO+ladder10K','SO + CH3OH A-ladder (LTE 10 K shape, free amplitude)')):
            s=get(b6,reg,model,'SO 6_5-5_4',sh,ous)
            if s is None: continue
            entry('SO 6_5-5_4 / HCO+ 3-2',ous,reg,'SO','HCO+ 3-2',s.W,s.sig,f'{ous}|SO 6_5-5_4|{reg}',hco3.W,hco3.sig,251.82577,267.55763,lab+f' shape={sh}')
        s=get(b6,reg,'SO_only','SO 6_5-5_4',sh,ous)
        entry('SO 6_5-5_4 / HCN 3-2',ous,reg,'SO','HCN 3-2',s.W,s.sig,f'{ous}|SO 6_5-5_4|{reg}',hcn3.W,hcn3.sig,251.82577,265.88618,f'SO only shape={sh}')
        for ln,f0 in (('SO2 6(3,3)-6(2,4)',254.2805),('CH3OH 2_0-1_-1 E',254.01538),('SiS 14-13',254.1031)):
            s=get(b6,reg,'spw20joint',ln,sh,ous)
            entry(f'{ln} / HCO+ 3-2',ous,reg,ln,'HCO+ 3-2',s.W,s.sig,f'{ous}|{ln}|{reg}',hco3.W,hco3.sig,f0,267.55763,f'joint with neighbours shape={sh}')
R=pd.DataFrame(rows); R.to_csv('/workspace/out/q5/ratios.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',500)
P=R[R.note.str.contains('shape=profile')]
print(P[P.ratio.isin(['CS 7-6 / HCN 4-3','CS 7-6 / HCO+ 4-3','SO 8_8-7_7 / HCO+ 4-3','SO 8_8-7_7 / CS 7-6'])][['ratio','region','W_num','sig_formal','sig_field','snr','status','W_den','R_K','R_K_err','R_K_UL3','R_K_UL5']].round(4).to_string())
print(P[P.ratio.str.startswith('SO 6_5')&P.note.str.startswith('SO only')][['ratio','ous','region','W_num','sig_formal','sig_field','snr','status','R_K','R_K_UL3','R_K_UL5']].round(4).to_string())
print(P[P.ratio.str.contains('SO2|SiS|CH3OH')][['ratio','ous','region','W_num','sig_used','snr','status','R_K','R_K_UL3']].round(4).to_string())
