# Q36 step 5: aggregate the Q22-with-Tbg runs (6 seeds each): chi2, BIC, logZ, P(Te>=50) (log-T prior), T medians.
import glob,numpy as np,pandas as pd
rows=[]
CFG=[('sh',t,g,t+'_'+g) for t,g in ((  'q36','q36cmb'),('q36','q36nominal'),('q36','q36dusty'),('q36','q36high'),('q36','_tbg5'),('q36e5','q36cmb (core CMB, ENV Tbg=5 K)'.split()[0]),('q36_envmid','q36nominal'),('q36_envhi','q36nominal'),('q36_envhi','q36high'))]+[('sep','q36',g,'q36_'+g) for g in ('q36cmb','q36nominal','q36dusty','q36high','_tbg5')]
for env,tag,g,lab in CFG:
    fs=sorted(glob.glob(f'/workspace/out/q36/runs/{tag}_{env}_{g}_s*.npz'))
    if not fs: continue
    d=pd.concat([pd.read_csv(f.replace('.npz','.csv')) for f in fs])
    X=[];
    for f in fs: z=np.load(f); X.append(dict(zip([str(n) for n in z['names']],z['X'].T)))
    p={k:np.concatenate([x[k] for x in X]) for k in X[0]}
    q=lambda a: tuple(np.percentile(a,[16,50,84]))
    r=dict(env=env,grid=tag+':'+g,nseed=len(fs),chi2_mean=d.chi2_prof.mean(),chi2_min=d.chi2_prof.min(),BIC_mean=d.BIC.mean(),logZ_mean=d.logZ.mean(),logZ_sd=d.logZ.std(),
      Tc_N=q(10**p['lTc_N']),Tc_S=q(10**p['lTc_S']),nc_N=q(p['nc_N']),nc_S=q(p['nc_S']))
    if env=='sh':
        r['Te']=q(10**p['lTe']); r['P_Te>=50']=float((10**p['lTe']>=50).mean()); r['P_Te>=50_perseed']=[round(float((10**x['lTe']>=50).mean()),2) for x in X]; r['ne']=q(p['ne'])
    else:
        for s in 'NS': r[f'Te_{s}']=q(10**p[f'lTe_{s}']); r[f'P_Te{s}>=50']=float((10**p[f'lTe_{s}']>=50).mean())
    rows.append(r)
df=pd.DataFrame(rows); pd.set_option('display.width',250,'display.max_columns',40,'display.max_colwidth',60)
f=lambda t: '%.1f/%.1f/%.1f'%t
for r in rows:
    print(r['env'],r['grid'],'n=%d chi2 mean %.1f (min %.1f) BIC %.1f logZ %.2f+-%.2f'%(r['nseed'],r['chi2_mean'],r['chi2_min'],r['BIC_mean'],r['logZ_mean'],r['logZ_sd']),'| Tc_N',f(r['Tc_N']),'Tc_S',f(r['Tc_S']),
      ('| Te %s P(Te>=50)=%.2f seeds %s'%(f(r['Te']),r['P_Te>=50'],r['P_Te>=50_perseed']) if r['env']=='sh' else '| Te_N %s Te_S %s P(TeN>=50)=%.2f P(TeS>=50)=%.2f'%(f(r['Te_N']),f(r['Te_S']),r['P_TeN>=50'],r['P_TeS>=50'])))
df.to_csv('/workspace/out/q36/q36_q22_summary.csv',index=False)
