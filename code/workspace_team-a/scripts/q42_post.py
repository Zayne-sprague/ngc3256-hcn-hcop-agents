# Q42 post-processing of out/q42/runs/*.npz (dynesty posteriors) -> out/q42/q42_summary.csv, q42_modes.txt
import glob,json,numpy as np,pandas as pd
from scipy.stats import gaussian_kde
def wq(x,w,q): o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return np.interp(q,c,x[o])
def modes(x,w,lo,hi):
    """weighted KDE (bw 0.06 dex*) in log10 T; peaks with prominence>10% of max; split at minima; returns [(Tpeak,weight)]"""
    g=np.linspace(lo,hi,400); k=gaussian_kde(x,weights=w,bw_method=0.25)(g)
    pk=[i for i in range(1,399) if k[i]>=k[i-1] and k[i]>k[i+1]]; pk=[i for i in pk if k[i]>0.1*k.max()]
    cuts=[]
    for a,b in zip(pk[:-1],pk[1:]):
        m=a+np.argmin(k[a:b+1])
        if k[m]<0.5*min(k[a],k[b]): cuts.append(m)   # real valley
    edges=[0]+cuts+[399]; res=[]
    for a,b in zip(edges[:-1],edges[1:]):
        sel=(x>=g[a])&(x<g[b]) if b<399 else (x>=g[a]); i=a+np.argmax(k[a:b+1]); res.append((10**g[i],float(w[sel].sum()/w.sum()),10**g[a],10**g[b]))
    return res
rows=[];txt=[]
for f in sorted(glob.glob('/workspace/out/q42/runs/*.npz')):
    d=np.load(f,allow_pickle=True); S,w,nm=d['S'],d['w'],list(d['names']); J=json.load(open(f.replace('.npz','.json')))
    P=lambda n:S[:,nm.index(n)]; TN,TS=10**P('lTc_N'),10**P('lTc_S')
    r=dict(J)
    for n,T in (('TcN',TN),('TcS',TS)): r[n+'_16'],r[n+'_50'],r[n+'_84']=wq(T,w,[.16,.5,.84])
    if 'lTe' in nm: Te=10**P('lTe'); r['Te_16'],r['Te_50'],r['Te_84']=wq(Te,w,[.16,.5,.84]); r['P50']=float(w[Te>=50].sum()); r['P_TeGtBoth']=float(w[(Te>TN)&(Te>TS)].sum())
    else:
        for z in 'NS': Te=10**P('lTe_'+z); r['P50_'+z]=float(w[Te>=50].sum())
    r['P_TcN_gt_TcS']=float(w[TN>TS].sum()); r['dlogT_med']=wq(P('lTc_N')-P('lTc_S'),w,.5); r['dlogT_16']=wq(P('lTc_N')-P('lTc_S'),w,.16); r['dlogT_84']=wq(P('lTc_N')-P('lTc_S'),w,.84)
    key=f.split('/')[-1][:-4]; txt.append(key)
    for n in ('lTc_N','lTc_S')+(('lTe',) if 'lTe' in nm else ()):
        txt.append(f'  {n}: '+'; '.join(f'T~{t:.0f} K [{a:.0f}-{b:.0f}] weight {wt:.2f}' for t,wt,a,b in modes(P(n),w,1.0,2.18)))
    rows.append(r)
D=pd.DataFrame(rows); D.to_csv('/workspace/out/q42/q42_summary.csv',index=False); open('/workspace/out/q42/q42_modes.txt','w').write('\n'.join(txt))
pd.set_option('display.width',280); pd.set_option('display.max_columns',40)
print(D.drop(columns=['sec','ncall','chi2min']).round(2).to_string()); print('\n'.join(txt))
