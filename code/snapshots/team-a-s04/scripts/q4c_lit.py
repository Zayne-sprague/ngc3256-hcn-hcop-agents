# Q4c-1: published 1-0 data (Harada+2018, ApJ 855, 49, arXiv:1801.05941, Table 13: 1.7" common beam, 12m only, pbcor) parsed from the arXiv HTML
# (/workspace/lit/harada/arxiv.html, downloaded with urllib) and compared with OUR 3-2 at the same 1.7" circular beam (Xb4d, convolved channel-wise, at the N and S positions).
# Outputs out/q4c/harada_table13_rows.csv, out/q4c/lit_vs_ours.csv
import sys,re; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from q4_lib import spectrum,vel,robust_sigma
tx=open('/workspace/lit/harada/harada.txt').read(); i0=tx.index('Table 13: Observed peak')
t=tx[i0:i0+60000]
# rows: position species transition restfreq Ipeak eIpeak I eI ; species/transition tokens appear as "HCN 1 − 0 1-0 88.634 1111.22 10.64 197.23 0.96"
rows=[]
for pos in ('N','S'):
    for sp,tr,rf in (('HCN','3-2','265.886'),('HCO + ','3-2','267.558'),('HCN','1-0','88.634'),('HCO + ','1-0','89.189')):
        m=re.search(r'(?<![A-Za-z0-9])%s %s [^ ]* ?%s? ?[^ ]* ?(?:%s[^0-9])?'%(pos,re.escape(sp),'',''),t)  # placeholder, replaced below
rows=[]
cur=None
toks=t.replace('−','-').split()
# find each "<pos> <species...> <J> <J> <restfreq> <peak> <epeak> <I> <eI>" by restfreq anchor
for pos,spn,rf in (('N','HCN 3-2','265.886'),('N','HCO+ 3-2','267.558'),('S','HCN 3-2','265.886'),('S','HCO+ 3-2','267.558'),
                   ('N','HCN 1-0','88.634'),('N','HCO+ 1-0','89.189'),('S','HCN 1-0','88.634'),('S','HCO+ 1-0','89.189')):
    sp=spn.split()[0].replace('+',' +'); 
    pat=r'\b%s (?:HCN|HCO \+) (?:\d) - (?:\d) \d-\d %s ([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+)'%(pos,re.escape(rf))
    s=' '.join(toks)
    cand=[m for m in re.finditer(r'(?:^| )(N|S|C|TNE|TSE|TSW|OS|SW) (HCN|HCO \+) (\d) - (\d) \d-\d %s ([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+)'%re.escape(rf),s) if m.group(1)==pos and m.group(2).replace(' ','')==spn.split()[0]]
    assert len(cand)==1,(pos,spn,len(cand))
    g=cand[0].groups(); rows.append(dict(pos=pos,line=spn,restfreq=rf,Ipeak_mK=float(g[4]),eIpeak_mK=float(g[5]),I_Kkms=float(g[6]),eI_Kkms=float(g[7])))
H=pd.DataFrame(rows); H.to_csv('/workspace/out/q4c/harada_table13_rows.csv',index=False); print(H.to_string())
# ours at 1.7" circular: integrated over +-165 km/s around v_sys=2780 (Harada convention: fixed window, no blend correction), baseline from 300-700 km/s off-line channels
TB=(1.7,1.7,0.); NU={'HCN 3-2':265.886434,'HCO+ 3-2':267.557626}; SPW={'HCN 3-2':'spw27','HCO+ 3-2':'spw25'}
out=[]
for reg in ('N','S'):
    for ln in NU:
        fr,s,tb,nb=spectrum('Xb4d',SPW[ln],reg,tb=TB); v=vel(fr,NU[ln]); dv=abs(np.median(np.diff(v)))
        win=np.abs(v-2780)<165; off=(np.abs(v-2780)>300)&(np.abs(v-2780)<700)
        c=np.polyfit(v[off],s[off],1); sd=robust_sigma(s[off]-np.polyval(c,v[off]))
        k=1.222e6/(NU[ln]**2*TB[0]*TB[1])   # K per Jy/beam at 1.7" (RJ)
        I=np.sum(s[win]-np.polyval(c,v[win]))*dv*k; eI=sd*np.sqrt(win.sum())*dv*k
        pk=(s-np.polyval(c,v))[win].max()*k
        out.append(dict(region=reg,line=ln,I_Kkms=I,eI=eI,Tpeak_K=pk,rms_K=sd*k,dv=dv))
O=pd.DataFrame(out); print(O.round(3).to_string())
def R(df,reg,a,b,col): 
    x=df[(df[df.columns[0]]==reg)]; return None
res=[]
for reg in ('N','S'):
    h=lambda ln:H[(H.pos==reg)&(H.line==ln)].iloc[0]; o=lambda ln:O[(O.region==reg)&(O.line==ln)].iloc[0]
    rr={}
    rr['H R10 (K km/s)']=(h('HCN 1-0').I_Kkms/h('HCO+ 1-0').I_Kkms, None)
    rr['H R32 (K km/s)']=(h('HCN 3-2').I_Kkms/h('HCO+ 3-2').I_Kkms,None)
    rr['ours R32 (K km/s, 1.7")']=(o('HCN 3-2').I_Kkms/o('HCO+ 3-2').I_Kkms,None)
    rr['H R10 (peak)']=(h('HCN 1-0').Ipeak_mK/h('HCO+ 1-0').Ipeak_mK,None)
    rr['H R32 (peak)']=(h('HCN 3-2').Ipeak_mK/h('HCO+ 3-2').Ipeak_mK,None)
    rr['ours R32 (peak, 1.7")']=(o('HCN 3-2').Tpeak_K/o('HCO+ 3-2').Tpeak_K,None)
    for sp in ('HCN','HCO+'):
        rr[f'H {sp} 3-2/1-0 (K km/s)']=(h(f'{sp} 3-2').I_Kkms/h(f'{sp} 1-0').I_Kkms,None)
        rr[f'ours3-2/H1-0 {sp} (K km/s)']=(o(f'{sp} 3-2').I_Kkms/h(f'{sp} 1-0').I_Kkms,None)
    for k,(v,_) in rr.items(): res.append(dict(region=reg,quantity=k,value=v))
    # errors: propagate stat errors of integrated intensities
    def e(r1,r2,e1,e2): return r1/r2*np.hypot(e1/r1,e2/r2)
    errs={'H R10 (K km/s)':e(h('HCN 1-0').I_Kkms,h('HCO+ 1-0').I_Kkms,h('HCN 1-0').eI_Kkms,h('HCO+ 1-0').eI_Kkms),
          'H R32 (K km/s)':e(h('HCN 3-2').I_Kkms,h('HCO+ 3-2').I_Kkms,h('HCN 3-2').eI_Kkms,h('HCO+ 3-2').eI_Kkms),
          'ours R32 (K km/s, 1.7")':e(o('HCN 3-2').I_Kkms,o('HCO+ 3-2').I_Kkms,o('HCN 3-2').eI,o('HCO+ 3-2').eI)}
    for sp in ('HCN','HCO+'):
        errs[f'H {sp} 3-2/1-0 (K km/s)']=e(h(f'{sp} 3-2').I_Kkms,h(f'{sp} 1-0').I_Kkms,h(f'{sp} 3-2').eI_Kkms,h(f'{sp} 1-0').eI_Kkms)
        errs[f'ours3-2/H1-0 {sp} (K km/s)']=e(o(f'{sp} 3-2').I_Kkms,h(f'{sp} 1-0').I_Kkms,o(f'{sp} 3-2').eI,h(f'{sp} 1-0').eI_Kkms)
    for r in res:
        if r['region']==reg: r['err_stat']=errs.get(r['quantity'],np.nan)
R=pd.DataFrame(res); R.round(4).to_csv('/workspace/out/q4c/lit_vs_ours.csv',index=False); O.round(4).to_csv('/workspace/out/q4c/ours_1p7arcsec.csv',index=False)
pd.set_option('display.width',200); print(R.round(3).to_string())
