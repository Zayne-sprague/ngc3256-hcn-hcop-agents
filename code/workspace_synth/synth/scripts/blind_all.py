# Blind line search on every region spectrum (peak pixel, native resolution), boxcar matched filter at 3 widths, threshold 5 sigma (robust MAD noise of filtered spectrum)
import sys,warnings,json; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from scipy.ndimage import uniform_filter1d
from curated import curated
from catalog import load
C=299792.458
d=np.load('/workspace/out/region_spectra.npz')
keys=sorted({k.rsplit('__',1)[0] for k in d.files if k.endswith('__pk')})
cat=load(300); cat=cat[~((cat.Eu==0)&(cat.logA==0))&(cat.logA>-4.3)]
cur=curated(300)
KNOWN={'HCN':[265.886434,354.505478,177.261111],'HCO+':[267.557626,356.734223,178.375056]}
rows=[]
for k in keys:
    ous,spw,reg=k.split('__'); s=d[k+'__pk']; fr=d[f'{ous}__{spw}__freq'] if f'{ous}__{spw}__freq' in d.files else None
    if fr is None: continue
    ok=np.isfinite(s); s=np.where(ok,s,0.)
    edge=(fr<fr.min()+0.03)|(fr>fr.max()-0.03)   # exclude 30 MHz at spw edges
    dv=abs(np.median(np.diff(fr)))/np.mean(fr)*C
    for wk in (100,200,300):
        nb=max(1,int(round(wk/dv))); sm=uniform_filter1d(s,nb,mode='nearest')
        m=np.ones_like(sm,bool)
        for _ in range(10):
            sg=1.4826*np.median(np.abs(sm[m]-np.median(sm[m]))); m=np.abs(sm-np.median(sm[m]))<2.5*sg
        base=np.median(sm[m]); sg=np.std(sm[m]); snr=(sm-base)/sg
        pk=[i for i in range(2,len(snr)-2) if (not edge[i]) and snr[i]>5 and snr[i]>=snr[i-1] and snr[i]>=snr[i+1]]
        for i in pk:
            fsky=fr[i]; rows.append(dict(ous=ous,spw=spw,region=reg,width=wk,fsky=fsky,peak_mJy=(sm[i]-base)*1e3,snr=snr[i],sig_mJy=sg*1e3))
df=pd.DataFrame(rows)
# merge features: same ous/region within 0.3*... take strongest per 0.08*fsky/c*300 km/s cluster
out=[]
for (ous,spw,reg),g in df.groupby(['ous','spw','region']):
    g=g.sort_values('snr',ascending=False); used=[]
    for _,r in g.iterrows():
        if any(abs(r.fsky-u)<400/C*r.fsky for u in used): continue
        used.append(r.fsky); out.append(r)
o=pd.DataFrame(out)
def ident(fsky,v0=2780.,dvtol=120.):
    f0lo=fsky/(1-(v0-dvtol)/C); f0hi=fsky/(1-(v0+dvtol)/C)
    s=cur[(cur.f0>=min(f0lo,f0hi))&(cur.f0<=max(f0lo,f0hi))&(cur.logA>-3.9)&(cur.Eu<200)].copy()
    s['v']=(1-fsky/s.f0)*C
    s=s.sort_values('logA',ascending=False).drop_duplicates('sp')
    return '; '.join(f'{r.sp} {str(r.qn).strip()[:22]} ({r.f0:.4f}, Eu={r.Eu:.0f}K, v={r.v:.0f})' for r in s.head(4).itertuples())
o['candidates']=[ident(f) for f in o.fsky]
o=o.sort_values(['ous','fsky','region'])
o.to_csv('/workspace/out/blind_detections_all.csv',index=False)
pd.set_option('display.width',300,'display.max_rows',500,'display.max_colwidth',140)
print(o[['ous','spw','region','fsky','width','peak_mJy','snr','candidates']].round(3).to_string())
