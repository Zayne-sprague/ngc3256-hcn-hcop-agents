# Q39 post-processing: per-configuration seed statistics from out/q39/runs/*.csv|npz  -> out/q39/q39_summary.csv, printed table
import glob,re,numpy as np,pandas as pd
rows=[]
for f in sorted(glob.glob('/workspace/out/q39/runs/*.csv')):
    d=pd.read_csv(f).iloc[0]; m=re.match(r'q39_([a-z0-9]+)_(noir|ir[a-z0-9]+)_(sh|sep)_q36([a-z]+)_s(\d+)\.csv',f.split('/')[-1])
    d['mode'],d['ir'],d['env'],d['scen'],d['seed']=m.group(1),m.group(2),m.group(3),m.group(4),int(m.group(5)); d['f']=f; rows.append(d)
D=pd.DataFrame(rows); out=[]
for (env,scen,mode,ir),g in D.groupby(['env','scen','mode','ir']):
    r=dict(env=env,scen=scen,mode=mode,ir=ir,nseed=len(g),logZ=g.logZ.mean(),logZ_sd=g.logZ.std(),chi2prof=g.chi2_prof.mean(),chi2prof_min=g.chi2_prof.min(),BIC=g.BIC.mean())
    X=[np.load(f.replace('.csv','.npz'),allow_pickle=True) for f in g.f]; nm=list(X[0]['names']); 
    # pool equal-weight posteriors of all seeds
    Xp=np.concatenate([x['X'] for x in X]); P=lambda n: Xp[:,nm.index(n)]
    for k in ('lTc_N','lTc_S'): r['Tc'+k[-2:]]=tuple(np.round(10**np.percentile(P(k),[16,50,84]),1))
    if 'lTe' in nm:
        ps=np.array([(10**x['X'][:,nm.index('lTe')]>=50).mean() for x in X]); r['P50_pooled']=ps.mean(); r['P50_seed_sd']=ps.std(ddof=1); r['P50_seed_sem']=ps.std(ddof=1)/np.sqrt(len(ps))
        r['Te']=tuple(np.round(10**np.percentile(P('lTe'),[16,50,84]),0))
    else:
        for n in ('N','S'):
            ps=np.array([(10**x['X'][:,nm.index('lTe_'+n)]>=50).mean() for x in X]); r['P50_'+n]=ps.mean(); r['P50_'+n+'_sd']=ps.std(ddof=1)
    out.append(r)
O=pd.DataFrame(out); O.to_csv('/workspace/out/q39/q39_summary.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(O.round(3).to_string())
