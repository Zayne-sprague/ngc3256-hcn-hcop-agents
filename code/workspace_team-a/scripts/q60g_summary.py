"""Q60 step 5: per-project summary of uv-plane N/S fits (out/q60/uvfit.csv): HCN N,S; HCO+ N,S; N/S; R43=HCN/HCO+ Tb ratio (flux ratio*(nu_HCO+/nu_HCN)^2);
contrast R43_N/R43_S. Xb47 EB X8300 is excluded (HCN N/S inverted, N flux 2-5x low: failed fit) but reported. Run: python3 -I q60g_summary.py"""
import pandas as pd, numpy as np, itertools, json
d=pd.read_csv('/workspace/out/q60/uvfit.csv',keep_default_na=False,na_values=['nan'])
K=(356.734223/354.505477)**2; BAD=['X8300']; out={}
def get(proj,line,us,mod,reg,hw=350):
    s=d[(d.proj==proj)&(d.line==line)&(d.uvset==us)&(d['model']==mod)&(d.comp==reg)&(d.hw==hw)&(~d.eb.isin(BAD))]
    return s.set_index('eb').I
lines=[]
for us,mod in [('all',2),('all',3),('B',2),('C',2),('D',2)]:
  for hw in (250,350,500):
    res={}
    for tag,(ph,po) in {'Xb47':('Xb47','Xb47'),'00493+00223':('00493','00223')}.items():
        R={}
        for reg in 'NS':
            h=get(ph,'HCN',us,mod,reg,hw); o=get(po,'HCO+',us,mod,reg,hw)
            if ph==po: r=(h/o).dropna()*K
            else: r=pd.Series([a/b*K for a,b in itertools.product(h.values,o.values)])
            R[reg]=r
        ns=get(ph,'HCN',us,mod,'N',hw)/get(ph,'HCN',us,mod,'S',hw)
        res[tag]=dict(HCN_N=get(ph,'HCN',us,mod,'N',hw).median(),HCN_S=get(ph,'HCN',us,mod,'S',hw).median(),
            HCNNS=ns.median(),HCNNS_sd=ns.std(),
            HCOp_N=get(po,'HCO+',us,mod,'N',hw).median(),HCOp_S=get(po,'HCO+',us,mod,'S',hw).median(),
            R43_N=R['N'].median(),R43_S=R['S'].median(),R43_sdN=R['N'].std(),R43_sdS=R['S'].std(),
            contrast=(R['N']/R['S'].median()).median() if ph==po else R['N'].median()/R['S'].median(),
            contrast_sd=(R['N']/R['S']).std() if ph==po else np.nan,nHCN=len(h),nHCO=len(o))
    for t,v in res.items(): lines.append(dict(uvset=us,model=mod,hw=hw,proj=t,**v))
T=pd.DataFrame(lines); T.to_csv('/workspace/out/q60/summary.csv',index=False)
pd.set_option('display.width',250)
print(T[T.hw==350].round(2).to_string())
print(T[(T.uvset=='all')&(T['model']==2)].round(2)[['hw','proj','contrast','R43_N','R43_S']].to_string())
# X8300 itself and 12m HCN (3.9-1.0)
print(d[(d.eb.isin(BAD))&(d.hw==350)&(d.uvset=='all')&(d['model']==2)][['eb','line','comp','I','err']])
