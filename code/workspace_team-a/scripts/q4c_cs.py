# Q4c-6: CS 7-6 as an abundance discriminator. Data: Q1 deblended component ratios CS 7-6/HCO+ 3-2 and CS 7-6/HCN 3-2 (Tb, K km/s, matched coarse beam; out/q1/q1_deblend_summary.csv, model 'nomissing',
# value + template-ensemble spread; + 14% cross-band calibration). Model: for every posterior draw (common Tkin, n(HCO+)) of the joint fit (classes A, B, F2; 'all' data; physT30) find the CS column
# that reproduces CS(7-6)/HCO+(3-2) = r_obs (RADEX CS, same Tkin and n as HCO+ , kappa=1) -> N(CS); then compare log[N(HCN)/N(CS)], log[N(HCO+)/N(CS)] and log[N(HCN)/N(HCO+)] between S and N
# with independent draws (posterior x measurement error). Also the 'raw' data ratios. Output out/q4c/cs_*.csv
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4c_joint as J
O=J.O; GCS=np.load(O+'grid4_cs.npz'); TRC=GCS['TR'][...,0]       # (T,n,N) CS 7-6
S=pd.read_csv('/workspace/out/q1/q1_deblend_summary.csv'); S=S[S.model=='nomissing']
def val(c,sp,q): r=S[(S.comp==c)&(S.species==sp)&(S.quantity==q)].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
rows=[]; rng=np.random.default_rng(5)
raw={}
for c in ('N','S'):
    for sp,q in (('HCO+','CS 7-6/HCO+ 3-2 (Tb)'),('HCN','CS 7-6/HCN 3-2 (Tb)')):
        v,s=val(c,sp,q); e=np.hypot(0.14*v,s); raw[(c,sp)]=(v,e); rows.append(dict(kind='data',region=c,quantity=q,value=v,err=e,spread_halfwidth=s))
# raw N/S contrasts
for sp in ('HCO+','HCN'):
    (vn,en),(vs,es)=raw[('N',sp)],raw[('S',sp)]; rows.append(dict(kind='data',region='S/N',quantity=f'CS7-6/{sp}3-2 ratio S over N',value=vs/vn,err=vs/vn*np.hypot(en/vn,es/vs),spread_halfwidth=np.nan))
rN,eN=raw[('N','HCO+')]; 
def draws(c,cl,nd=20000):
    P=J.posterior(c,'all',[cl]); w=P['w']/P['w'].sum(); idx=rng.choice(len(w),size=nd,p=w)
    T=P['T'][idx]; nO=P['logn_O'][idx]; NH=P['logNH'][idx]; NO=P['logNO'][idx]
    it=np.searchsorted(J.TK,T); inn=np.round((nO-J.LN[0])/0.25).astype(int)
    r,e=raw[(c,'HCO+')]; target=r*np.exp(rng.normal(0,e/r,nd))    # measured CS7-6/HCO+3-2 incl. error
    # HCO+ 3-2 TR at the draw
    iN=np.round((NO-J.LNC[0])/0.25).astype(int); trO=J.TO[it,inn,iN,2]
    trc=TRC[it,inn,:]                         # (nd,23)
    want=target*trO*(P['kappa'][idx]*0+1.0)
    # N_CS by interpolation on the monotonic rising part
    logNcs=np.full(nd,np.nan)
    for i in range(nd):
        t=trc[i]; ok=np.isfinite(t)
        if ok.sum()<3: continue
        m=np.maximum.accumulate(np.nan_to_num(t))   # monotone envelope
        if want[i]<m[0] or want[i]>m[-1]: continue
        # first crossing
        j=np.searchsorted(m,want[i]); j=max(j,1); x0,x1=m[j-1],m[j]; f=0 if x1==x0 else (want[i]-x0)/(x1-x0); logNcs[i]=J.LNC[j-1]+f*0.25
    return dict(logNH=NH,logNO=NO,logNcs=logNcs,T=T,n=nO,tau=P['tau32H'][idx])
res={}
for cl in ('A both thick, kappa=1','B both thick, same Tex (0.8-1.25), kappa free','F2 separate n, kappa=1, both thick'):
    D={c:draws(c,cl) for c in ('N','S')}
    for c in D: rows.append(dict(kind='model',cls=cl,region=c,quantity='fraction of draws where CS can be matched',value=float(np.isfinite(D[c]['logNcs']).mean())))
    def delta(a,b):
        # independent draws S vs N, using finite values
        s=np.array([x for x in a if np.isfinite(x)]); n=np.array([x for x in b if np.isfinite(x)])
        k=min(len(s),len(n)); d=rng.permutation(s)[:k]-rng.permutation(n)[:k]; return d
    for nm,fa in {'log[N(HCN)/N(CS)]':(lambda D:D['logNH']-D['logNcs']),'log[N(HCO+)/N(CS)]':(lambda D:D['logNO']-D['logNcs']),'log[N(HCN)/N(HCO+)]':(lambda D:D['logNH']-D['logNO']),'log N(CS)':(lambda D:D['logNcs'])}.items():
        d=delta(fa(D['S']),fa(D['N'])); p=np.percentile(d,[2.5,16,50,84,97.5])
        for c in ('N','S'):
            x=fa(D[c]); x=x[np.isfinite(x)]; rows.append(dict(kind='model',cls=cl,region=c,quantity=nm,value=np.median(x),lo16=np.percentile(x,16),hi84=np.percentile(x,84)))
        rows.append(dict(kind='model',cls=cl,region='S minus N',quantity=nm,value=p[2],lo16=p[1],hi84=p[3],lo95=p[0],hi95=p[4],P_negative=float((d<0).mean()),sigma=abs(np.median(d))/max(np.std(d),1e-9)))
R=pd.DataFrame(rows); R.round(4).to_csv(O+'cs_results.csv',index=False); pd.set_option('display.width',250,'display.max_colwidth',60); print(R.round(3).to_string())
