# Q22 round 2 post-processing: evidence/chi2 tables over seeds, prior-robust P(T) (flat-log-T vs flat-linear-T by importance reweighting), masses & X consistency.
import numpy as np, pandas as pd, glob, os
from scipy.ndimage import map_coordinates
R='/workspace/out/q22b/runs/'
rows=[]
for f in sorted(glob.glob(R+'*.csv')):
    d=pd.read_csv(f).iloc[0].to_dict(); d['file']=os.path.basename(f); rows.append(d)
D=pd.DataFrame(rows); D['drop']=D['drop'].fillna('').astype(str).str.replace('nan','')
g=D.groupby(['tag','env','grid']).agg(n=('seed','count'),k=('npar','first'),ndat=('ndat','first'),chi2=('chi2_prof','mean'),chi2N=('chi2_N','mean'),chi2S=('chi2_S','mean'),
   chi2_10N=('chi2_10_N','mean'),chi2_10S=('chi2_10_S','mean'),lnZ=('logZ','mean'),lnZsd=('logZ','std'),BIC=('BIC','mean'),PTe50=('P(Te>=50)','mean')).reset_index()
g.to_csv('/workspace/out/q22b/q22b_summary.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
print(g.round(2).to_string())
# ---- posteriors
LTN=['lTc_N','lTc_S','lTe','lTe_N','lTe_S']
def load(tag,env,grid):
    fs=sorted(glob.glob(R+f'{tag}_{env}_{grid}_s*.npz')); Z=[np.load(f,allow_pickle=True) for f in fs]
    if not Z: return None,None
    names=list(Z[0]['names']); return np.vstack([z['X'] for z in Z]),names
def wq(x,w,q): o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return np.interp(q,c,x[o])
out=[]
for grid in ('tbg2.73','_tbg5'):
  for env in ('sh','shf','sep','shH','same'):
    X,names=load('base',env,grid)
    if X is None: continue
    ix={n:i for i,n in enumerate(names)}; Ts=[n for n in names if n in LTN]
    for prior in ('logT','linT'):
        w=np.ones(len(X))
        if prior=='linT':
            for n in Ts:
                w=w*10**X[:,ix[n]]
        w/=w.sum()
        r=dict(grid=grid,env=env,prior=prior)
        for n in Ts:
            t=10**X[:,ix[n]]; r[n[2:]+'_med']=wq(t,w,.5); r[n[2:]+'_16']=wq(t,w,.16); r[n[2:]+'_84']=wq(t,w,.84)
            if n.startswith('lTe'): r['P('+n[2:]+'>=50)']=float(w[t>=50].sum())
        if 'lTe' in ix:
            te=10**X[:,ix['lTe']]; r['P(Te>=50)']=float(w[te>=50].sum()); r['P(Te>=30)']=float(w[te>=30].sum()); r['P(Te>=100)']=float(w[te>=100].sum())
            for c in ('N','S'): r['P(Te>Tc_'+c+')']=float(w[X[:,ix['lTe']]>X[:,ix['lTc_'+c]]].sum())
            r['P(Tc_N>=Tc_S)']=float(w[X[:,ix['lTc_N']]>=X[:,ix['lTc_S']]].sum())
            r['corr_lTe_lTcS']=float(np.corrcoef(X[:,ix['lTe']],X[:,ix['lTc_S']])[0,1])
        if 'lTe_S' in ix:
            for c in ('N','S'): r['P(Te_'+c+'>=50)']=float(w[10**X[:,ix['lTe_'+c]]>=50].sum())
        out.append(r)
P=pd.DataFrame(out); P.to_csv('/workspace/out/q22b/q22b_PT.csv',index=False); print(P.round(2).to_string())
# prior of P(Te>=50) under each T prior over [10,150]
print('prior P(T>=50): logT %.2f  linT %.2f'%((np.log(150)-np.log(50))/(np.log(150)-np.log(10)),(150-50)/(150-10)))
# ---- masses / X consistency (shared envelope, base, 2.73)
G={m:np.load(f'/workspace/out/q4c/grid4_{m}.npz') for m in ('hcn','hco+')}
mu=2.8*1.6605e-24; pc=3.086e18; Msun=1.989e33; PC=194.
SIZE={'S':0.455,'N':0.507}; MDYN={'S':0.58e8}
for env in ('sh','sep'):
    X,names=load('base',env,'tbg2.73'); ix={n:i for i,n in enumerate(names)}
    for Xab in (1e-9,3e-9,1e-8,3e-8):
        row=dict(env=env,X=Xab)
        for r in ('N','S'):
            NHc=10**X[:,ix['NHc_'+r]]/Xab; nc=10**X[:,ix['nc_'+r]]; fc=10**X[:,ix['lfc_'+r]]
            Lc=NHc/nc/pc; Ac=np.pi*(SIZE[r]/2*PC*pc)**2
            Mc=fc*NHc*mu*Ac/Msun
            ne=10**X[:,ix['ne' if env=='sh' else 'ne_'+r]]; NHe=10**X[:,ix['NHe' if env=='sh' else 'NHe_'+r]]/Xab; we=10**X[:,ix['lwe' if env=='sh' else 'lwe_'+r]]
            Me=we*NHe*mu*np.pi*(5*PC*pc)**2/Msun; Le=NHe/ne/pc
            ok=(Lc<=66)&(Le<=500)
            row.update({f'Lc_{r}_med':np.median(Lc),f'fok_{r}':ok.mean(),f'Mc_{r}':np.median(Mc),f'Me10_{r}':np.median(Me),f'fcore_{r}':np.median(Mc/(Mc+Me)),f'fcore_{r}_16':np.percentile(Mc/(Mc+Me),16),f'fcore_{r}_84':np.percentile(Mc/(Mc+Me),84),f'P(Mc<=Mdyn)_{r}':float((Mc<=MDYN['S']).mean()) if r=='S' else np.nan})
        print({k:(float('%.3g'%v) if isinstance(v,float) else v) for k,v in row.items()})
        pd.DataFrame([row]).to_csv('/workspace/out/q22b/q22b_mass.csv',mode='a',header=not os.path.exists('/workspace/out/q22b/q22b_mass.csv'),index=False)
