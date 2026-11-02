# Q1: matched-beam line fluxes and ratios at the p0001 regions.
#  sets: 1.05" circular (Xb4d, Xb4b), 6.5" (Xb4f, Xb4d, Xb4b, Xb47), 10" (Xb55, Xb4f/d/b, Xb47). Moment-0 over VSYS+-350 km/s, pbcor, convolved with the
#  image-plane kernel (matched.kernel), value at the region position (= beam-averaged brightness of the common beam), noise from empirical null windows (q1_lib).
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import pandas as pd
reg=pd.read_csv('/workspace/out/regions.csv'); reg=reg[reg.region!='Disk']
POS={r.region:(r.dRA_arcsec,r.dDec_arcsec) for r in reg.itertuples()}
SETS={1.05:['HCN 3-2@Xb4d','HCO+ 3-2@Xb4d','SO 6_5-5_4@Xb4d','HCN 3-2@Xb4b','HCO+ 3-2@Xb4b','SO 6_5-5_4@Xb4b'],
      6.5:['HCN 3-2@Xb4f','HCO+ 3-2@Xb4f','SO 6_5-5_4@Xb4f','HCN 3-2@Xb4d','HCO+ 3-2@Xb4d','SO 6_5-5_4@Xb4d','HCN 3-2@Xb4b','HCO+ 3-2@Xb4b','HCN 4-3','HCO+ 4-3','CS 7-6'],
      10.0:['HCN 2-1','HCO+ 2-1','HCN 3-2@Xb4f','HCO+ 3-2@Xb4f','HCN 3-2@Xb4d','HCO+ 3-2@Xb4d','HCN 4-3','HCO+ 4-3','CS 7-6']}
CROP=11.0    # arcsec half-size kept from the 1.3 GB Xb4b cubes (>= all regions + smoothing margin)
def nu_sky(key): return LINES[key][2]/(1+VSYS/C)   # approx sky frequency for the Tb conversion
rows=[]; cache={}
def get(key):
    if key not in cache: cache[key]=load(key,crop=CROP if 'Xb4b' in key else None)
    return cache[key]
for fw,keys in SETS.items():
    for key in keys:
        d=get(key); sm=smooth(d,fw)
        if sm is None: print('skip',key,fw); continue
        prof=noise_profile(d,fw)
        for rn,(x,y) in POS.items():
            r=np.hypot(x,y); I=sample(sm,d,x,y)
            if r>np.hypot(d['dx'].max(),d['dy'].max())-1 or abs(x)>np.abs(d['dx']).max() or abs(y)>np.abs(d['dy']).max(): continue
            e=sigma_at(prof,r); kT=Tb_factor(nu_sky(key),fw)
            disc=np.hypot(d['dx']-x,d['dy']-y)<fw/2; Id=float(np.nanmean(sm[disc])) if disc.sum()>0 else np.nan
            for ap,Iv in (('point',I),('disc',Id)):
                rows.append(dict(beam=fw,key=key,ous=LINES[key][0],line=key.split('@')[0],region=rn,ap=ap,I=Iv,eI=e,snr=Iv/e,I_K=Iv*kT,eI_K=e*kT,lim3=3*e,lim3_K=3*e*kT,n_null=len(d['nulls'])))
T=pd.DataFrame(rows); T.to_csv('/workspace/out/q1/q1_matched_fluxes.csv',index=False); print(T[T.beam==6.5].round(3).head(30).to_string())

# ---------- ratio table ----------
SYS_SAME=0.05; SYS_X=0.14   # window/baseline+calibration systematic (same OUS, HCN/HCO+ in the same execution) ; cross-band (2 x 10% flux scale)
def val(beam,key,region,ap):
    q=T[(T.beam==beam)&(T.key==key)&(T.region==region)&(T.ap==ap)]
    return None if q.empty else q.iloc[0]
AP='point'
def ratio(beam,nk,dk,region,tb=True,name=''):
    ap=AP
    a,b=val(beam,nk,region,ap),val(beam,dk,region,ap)
    if a is None or b is None: return None
    ia,ea,ib,eb=(a.I_K,a.eI_K,b.I_K,b.eI_K) if tb else (a.I,a.eI,b.I,b.eI)
    sa,sb=ia/ea,ib/eb; sys_=SYS_SAME if LINES[nk][0]==LINES[dk][0] else SYS_X
    row=dict(beam=beam,ap=ap,ratio=name,ous_num=LINES[nk][0],ous_den=LINES[dk][0],num=nk,den=dk,region=region,num_snr=sa,den_snr=sb)
    if sa>=3 and sb>=3:
        R=ia/ib; es=abs(R)*np.hypot(ea/ia,eb/ib); row.update(R=R,eR_stat=es,eR_tot=np.hypot(es,R*sys_),kind='value')
    elif sb>=3: row.update(R=3*ea/ib,eR_stat=np.nan,eR_tot=np.nan,kind='upper limit (3sig)')
    elif sa>=3: row.update(R=ia/(3*eb),eR_stat=np.nan,eR_tot=np.nan,kind='lower limit (3sig)')
    else: row.update(R=np.nan,kind='undetermined')
    return row
RR=[]
for b,ap_ in [(b,a) for b in (1.05,6.5,10.0) for a in ('point','disc')]:
    AP=ap_
    for rn in POS:
        for ous in ('Xb4f','Xb4d','Xb4b'):
            RR.append(ratio(b,f'HCN 3-2@{ous}',f'HCO+ 3-2@{ous}',rn,name='HCN/HCO+ 3-2 (Tb)'))
            RR.append(ratio(b,f'SO 6_5-5_4@{ous}',f'HCO+ 3-2@{ous}',rn,name='SO 6_5-5_4/HCO+ 3-2 (Tb)'))
        RR.append(ratio(b,'HCN 2-1','HCO+ 2-1',rn,name='HCN/HCO+ 2-1 (Tb)'))
        RR.append(ratio(b,'HCN 4-3','HCO+ 4-3',rn,name='HCN/HCO+ 4-3 (Tb)'))
        RR.append(ratio(b,'CS 7-6','HCN 4-3',rn,name='CS 7-6/HCN 4-3 (Tb)'))
        RR.append(ratio(b,'CS 7-6','HCO+ 4-3',rn,name='CS 7-6/HCO+ 4-3 (Tb)'))
        for sp in ('HCN','HCO+'):
            RR.append(ratio(b,f'{sp} 3-2@Xb4f',f'{sp} 2-1',rn,name=f'{sp} 3-2/2-1 (Tb)'))
            RR.append(ratio(b,f'{sp} 4-3',f'{sp} 3-2@Xb4f',rn,name=f'{sp} 4-3/3-2 (Tb)'))
            RR.append(ratio(b,f'{sp} 4-3',f'{sp} 2-1',rn,name=f'{sp} 4-3/2-1 (Tb)'))
R=pd.DataFrame([r for r in RR if r is not None]); R=R[R.kind!='undetermined'] if 'kind' in R else R
R.to_csv('/workspace/out/q1/q1_ratios.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',400)
print(R[R.ratio.str.contains('HCN/HCO')&(R.ap=='disc')].round(3)[['beam','ratio','region','R','eR_stat','eR_tot','kind','num_snr','den_snr']].to_string())
