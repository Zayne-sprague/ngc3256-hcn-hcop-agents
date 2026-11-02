# Q5 round 2: Harada+2018 (ApJ 855, 49; arXiv:1801.05941) Table 13 rows (parsed from the arXiv PDF text: lit/harada2018_1801.05941.txt -> lit/harada_table13_parsed.csv,
# cross-checked against team-a's /snapshots/team-a-s04/out/q4c/harada_table13_rows.csv: identical N,S rows) vs this work at the same 1.7" resolution (out/q5/harada_compare_1p7.csv).
# Ratios are K km/s ratios; ours noise = max(formal, empirical sliding-template MAD) ; 10% flux scale (each of Harada, ALMA) added in quadrature to line-ratio enhancement factors.
import pandas as pd, numpy as np, json
H=pd.read_csv('/workspace/lit/harada_table13_parsed.csv'); O=pd.read_csv('/workspace/out/q5/harada_compare_1p7.csv').set_index('pos')
def h(pos,sp,f): 
    r=H[(H.pos==pos)&(H.species_transition.str.replace(' ','').str.startswith(sp))&(abs(H.restfreq-f)<0.01)]
    return (float(r.I_Kkms.iloc[0]),float(r.eI_Kkms.iloc[0])) if len(r) else (np.nan,np.nan)
def ratio(a,b): 
    if np.isnan(a[0]) or np.isnan(b[0]): return (np.nan,np.nan)
    return (a[0]/b[0],abs(a[0]/b[0])*np.hypot(a[1]/a[0] if a[0]>0 else np.nan,b[1]/b[0]))
rows=[]
for pos in ['N','S','C','TNE','TSE','TSW','OS','SW']:
    hc=h(pos,'HCO+',267.558); hn=h(pos,'HCN',265.886); hc1=h(pos,'HCO+',89.189)
    for nm,sp,f in (('CH3OH 241.79 5k-4k','CH3OH',241.791),('CH3OH 254.015 2_0-1_-1 E','CH3OH',254.015),('CH3OH 261.806','CH3OH',261.806),('CS 5-4','CS',244.936),('SO 219.949','SO',219.949)):
        v=h(pos,sp,f); rr=ratio(v,hc); rows.append(dict(src='Harada+2018 (1.7")',pos=pos,numerator=nm,I_num=v[0],eI_num=v[1],den='HCO+ 3-2',I_den=hc[0],R_K=rr[0],eR=rr[1],note='not in Table 13 (<2sigma peak)' if np.isnan(v[0]) else ''))
    for nm,sp,f in (('CH3OH 96.74 2k-1k','CH3OH',96.741),('SO 99.3 2_1-3_2','SO',99.3),('SiO 2-1','SiO',86.847),('CS 2-1','CS',97.981)):
        v=h(pos,sp,f); rr=ratio(v,hc1); rows.append(dict(src='Harada+2018 (1.7")',pos=pos,numerator=nm,I_num=v[0],eI_num=v[1],den='HCO+ 1-0',I_den=hc1[0],R_K=rr[0],eR=rr[1],note='not in Table 13 (<2sigma peak)' if np.isnan(v[0]) else ''))
conv=1.222e6/(1.7**2)
for pos in ['N','S','Bridge','Sblue','W']:
    o=O.loc[pos]; hco=o['HCO+ 3-2 K']; nu_h=265.1; 
    for nm,W,sf,mad,med,nu in (('CH3OH 254.015 2_0-1_-1 E',o['CH3OH254 W'],o['CH3OH254 sf'],o['CH3OH254 mad'],o['CH3OH254 median'],251.7),('SO 251.826 (SO only)',o['SO W'],o['SO sf'],o['SO mad'],o['SO median'],249.5),('SO 251.826 (+CH3OH A ladder)',o['SO+lad W'],o['SO+lad sf'],o['SO mad'],o['SO median'],249.5)):
        s=max(sf,mad); Ik=conv*W/nu**2; eIk=conv*s/nu**2; rows.append(dict(src='this work (Xb4d smoothed to 1.7")',pos=pos,numerator=nm,I_num=Ik,eI_num=eIk,den='HCO+ 3-2',I_den=hco,R_K=Ik/hco,eR=eIk/hco,note=f'snr={W/s:.1f}; 3sigma UL R_K<{(max(W,0)+3*s)*conv/nu**2/hco:.4f}; empirical off-line median {med:.3f} Jy/beam km/s'))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q5/harada_vs_ours.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
print(df[df.den=='HCO+ 3-2'].round(4).to_string())
# ours vs Harada HCN/HCO+ 3-2 and HCO+ / HCN K km/s
for pos in ['N','S']:
    o=O.loc[pos]; print(pos,'HCO+3-2 K',round(o['HCO+ 3-2 K'],1),'vs',h(pos,'HCO+',267.558),'; HCN3-2',round(o['HCN 3-2 K'],1),'vs',h(pos,'HCN',265.886),'; HCN/HCO+ ours',round(o['HCN 3-2 K']/o['HCO+ 3-2 K'],3),'H',round(h(pos,'HCN',265.886)[0]/h(pos,'HCO+',267.558)[0],3))
