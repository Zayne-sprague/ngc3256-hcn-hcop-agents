"""Q56 r2: summarise s344.json -> error budget for R32/10, S/N, contrast. usage: python3 -I s345_q56_summary.py -> out/q56/s345.txt/json"""
import json, numpy as np
R=json.load(open('/workspace/out/q56/s344.json')); out=[]; sm={}
def P(t): print(t); out.append(t)
def r(th,w,b,uv,rn,lab,ln): a=R[f'{th}|{w}|{b}|{uv}|{rn}|{lab}|{ln}']; return a[1]/a[0]
def stat(vals): v=np.array(vals); return v.mean(),v.std(ddof=1) if len(v)>1 else 0,v.min(),v.max()
P('== A. reproduction + estimator (uv unfiltered, lin baseline, w=100): R32/10 [N, S, S/N]')
for th in (2.5,3.5):
  for lab in ('nom','peak','ap0.5','ap1'):
    row=f'{th}" {lab:6s}'
    for ln in ('HCN','HCO+'):
        n,s=r(th,100,'lin',0,'N',lab,ln),r(th,100,'lin',0,'S',lab,ln); row+=f' | {ln}: {n:.3f} {s:.3f} S/N {s/n:.3f}'
    P(row)
P('\n== B. pixel-position jackknife (3x3 grid, +-0.25 theta), uv unfiltered, w=100 lin: S/N and contrast')
for th in (2.5,3.5):
  for ln in ('HCN','HCO+'):
    sn=[];rn=[];rs=[]
    for ix in (-1,0,1):
      for iy in (-1,0,1):
        n=r(th,100,'lin',0,'N',f'jk{ix}{iy}',ln); s=r(th,100,'lin',0,'S',f'jk{ix}{iy}',ln); rn.append(n);rs.append(s);sn.append(s/n)
    P(f'{th}" {ln}: R32/10 N {np.mean(rn):.3f}+-{np.std(rn,ddof=1):.3f} [{min(rn):.3f},{max(rn):.3f}]  S {np.mean(rs):.3f}+-{np.std(rs,ddof=1):.3f} [{min(rs):.3f},{max(rs):.3f}]  S/N {np.mean(sn):.3f}+-{np.std(sn,ddof=1):.3f} [{min(sn):.3f},{max(sn):.3f}]')
P('\n== C. full variant table: S/N (R32/10) and D=(S/N)_HCO+/(S/N)_HCN, by uv filter')
for uv in (0,7.7,13.0):
  for th in (2.5,3.5):
    for est in ('nom','peak','ap0.5','ap1'):
      for w in (100,150,350):
        for b in ('lin','const'):
          snH=r(th,w,b,uv,'S',est,'HCN')/r(th,w,b,uv,'N',est,'HCN'); snO=r(th,w,b,uv,'S',est,'HCO+')/r(th,w,b,uv,'N',est,'HCO+')
          sm[f'{uv}|{th}|{est}|{w}|{b}']=(snH,snO,snO/snH)
P('uv  th   est  : HCN S/N mean[min,max] over (w,baseline) | HCO+ S/N | D')
for uv in (0,7.7,13.0):
  for th in (2.5,3.5):
    for est in ('nom','peak','ap0.5','ap1'):
      k=[v for kk,v in sm.items() if kk.startswith(f'{uv}|{th}|{est}|')]
      a=np.array(k); P(f'{uv:5.1f} {th} {est:6s}: HCN {a[:,0].mean():.3f} [{a[:,0].min():.3f},{a[:,0].max():.3f}] | HCO+ {a[:,1].mean():.3f} [{a[:,1].min():.3f},{a[:,1].max():.3f}] | D {a[:,2].mean():.3f} [{a[:,2].min():.3f},{a[:,2].max():.3f}]')
P('\n== D. absolute R32/10 N,S under uv filters (nom, w100 lin)')
for uv in (0,7.7,13.0):
  for th in (2.5,3.5):
    P(f'uv>{uv} {th}": '+' | '.join(f'{ln} N {r(th,100,"lin",uv,"N","nom",ln):.3f} S {r(th,100,"lin",uv,"S","nom",ln):.3f}' for ln in ('HCN','HCO+')))
open('/workspace/out/q56/s345.txt','w').write('\n'.join(out)); json.dump(sm,open('/workspace/out/q56/s345.json','w'))
# ---- E. pooled error budget and contrast table
P('\n== E. pooled systematic budget (all variants: theta, estimator [nom,peak,ap0.5,ap1], window, baseline; uv filter separately)')
def contrast(th,w,b,uv,est,J):
    k=0 if J=='10' else 1
    g=lambda ln,rn: R[f'{th}|{w}|{b}|{uv}|{rn}|{est}|{ln}'][k]
    return (g('HCN','N')/g('HCO+','N'))/(g('HCN','S')/g('HCO+','S'))
tab={}
for uvs,lab in (((0,),'full uv'),((7.7,13.0),'uv-filtered 7.7/13 kl'),((0,7.7,13.0),'all')):
    a=np.array([v for kk,v in sm.items() if float(kk.split('|')[0]) in uvs])
    P(f'{lab}: HCN S/N {a[:,0].mean():.3f}+-{a[:,0].std(ddof=1):.3f} (range {a[:,0].min():.2f}-{a[:,0].max():.2f}); HCO+ S/N {a[:,1].mean():.3f}+-{a[:,1].std(ddof=1):.3f}; D {a[:,2].mean():.3f}+-{a[:,2].std(ddof=1):.3f} (range {a[:,2].min():.2f}-{a[:,2].max():.2f}); n={len(a)}')
    tab[lab]=[a[:,0].mean(),a[:,0].std(ddof=1),a[:,1].mean(),a[:,1].std(ddof=1),a[:,2].mean(),a[:,2].std(ddof=1)]
P('N/S contrast of HCN/HCO+ (mean over window/baseline; spread in brackets), uv unfiltered:')
for th in (2.5,3.5):
  for est in ('nom','peak','ap0.5','ap1'):
    c10=[contrast(th,w,b,0,est,'10') for w in (100,150,350) for b in ('lin','const')]; c32=[contrast(th,w,b,0,est,'32') for w in (100,150,350) for b in ('lin','const')]
    P(f'{th}" {est:6s}: 1-0 {np.mean(c10):.2f} [{min(c10):.2f},{max(c10):.2f}]  3-2 {np.mean(c32):.2f} [{min(c32):.2f},{max(c32):.2f}]  3-2/1-0 {np.mean(c32)/np.mean(c10):.3f}')
# jackknife-informed S/N error for 2.5 and 3.5 (nominal, w100)
json.dump(tab,open('/workspace/out/q56/s345_budget.json','w'))
open('/workspace/out/q56/s345.txt','w').write('\n'.join(out))
