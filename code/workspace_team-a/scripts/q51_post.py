# Q51 post-processing of out/q51/runs: Tc, opacity/density branches, H13CN survival (Z ratio and importance), core mass vs Mdyn, X_crit, clumpiness (path length / core size).
# conventions: mu=2.8 m_H, D=40 Mpc (194 pc/")), X=N(HCN)/N(H2), M_core=fc*N_H2*mu*pi*(FWHM_HCN3-2/2)^2 (as Q22e/Q49), M_env = we*N_H2(env)*mu*A_beam with A=1.133*5.73*3.54 arcsec^2 (Xb4f beam)
# Mdyn(48 pc, +pressure not included) N 1.7(1.2-2.8)e8, S 0.58(0.41-0.78)e8 Msun (Q5f, report item 160); P evaluated at Mdyn and x0.5/x2 (S has a factor-2 systematic)
import sys,os,glob,json,warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
mu=2.8*1.6605e-24; pc=3.086e18; Msun=1.989e33; PC=194.; MDYN={'N':1.7e8,'S':0.58e8}; ABEAM=1.133*5.73*3.54
R='/workspace/out/q51/runs/'
def wq(x,w,q):
    o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return np.interp(q,c,x[o])
J=[json.load(open(f)) for f in sorted(glob.glob(R+'q51_*.json'))]
ctx={}
def head(model,reg,geom,dset):
    """exec the likelihood module head to get the geometry-dependent compact size S['hcn'] and TAU32/unpack"""
    k=(model,reg,geom,dset)
    if k in ctx: return ctx[k]
    env=dict(os.environ)
    if geom=='geo': os.environ.update(Q12_GEOM='1',Q12_U='1',Q12_SK='.15')
    else:
        for e in('Q12_GEOM','Q12_U','Q12_SK'): os.environ.pop(e,None)
    src=open('/workspace/scripts/q12c_twocomp.py').read().split('def smc(')[0]
    g={'__name__':'q12c'}; sa=sys.argv; sys.argv=['q12c',model,reg,dset,'g273','1']; exec(compile(src,'q12c','exec'),g); sys.argv=sa
    os.environ.clear(); os.environ.update(env); ctx[k]=g; return g
rows=[]
for j in J:
    d=np.load(R+f"q51_{j['model']}_{j['reg']}_{j['geom']}_h{j['h13']}_{j['dat']}_{j['dset']}_s{j['seed']}.npz",allow_pickle=True)
    S,w,nm=d['S'],d['w']/d['w'].sum(),list(d['names']); g=head(j['model'],j['reg'],j['geom'],j['dset'])
    P={n:S[:,i] for i,n in enumerate(nm)}
    p={n:S[:,i] for i,n in enumerate(nm)}; p=g['unpack'](S)
    reg=j['reg']; sz=g['S']['hcn']*np.exp(P.get('lk3',0*S[:,0]))        # compact HCN 3-2 FWHM (arcsec)
    fc=10**P['lfc']; NHc=10**P['NHc']; nc=10**P['nc']
    Mc1=fc*NHc*mu*np.pi*(sz/2*PC*pc)**2/Msun       # M_core at X=1
    L1=NHc/nc/pc                                      # path at X=1 (pc)
    two=j['model']!='1z'
    r=dict(model=j['model'],reg=reg,geom=j['geom'],h13=j['h13'],dat=j['dat'],dset=j['dset'],seed=j['seed'],logZ=j['logZ'],ess=j['ess'],chi2min=j['chi2min'],sec=j['sec'])
    for n,q in (('Tc',10**P['lTc']),('lognc',P['nc']),('logNHc',P['NHc']),('lfc',P['lfc'])):
        for pp in (16,50,84): r[f'{n}_{pp}']=wq(q,w,pp/100)
    r['P_T50']=float(w[10**P['lTc']>=50].sum()); r['P_T30']=float(w[10**P['lTc']<=30].sum())
    tau=g['TAU32']('hcn',P['lTc'],P['nc'],P['NHc']); 
    for pp in (16,50,84): r[f'tau32_{pp}']=wq(tau,w,pp/100)
    r['P_tau10']=float(w[tau>10].sum()); r['P_n6']=float(w[P['nc']>=6].sum())
    lx=P['NHc']-P['NOc']; r['logX_50']=wq(lx,w,.5)
    r['size_50']=wq(sz,w,.5)
    for pp in (16,50,84): r[f'Mc1_{pp}']=wq(Mc1,w,pp/100)     # M_core*X (Msun)
    xc=Mc1/MDYN[reg]                                         # X needed for Mcore<=Mdyn
    for pp in (16,50,84): r[f'Xcrit_{pp}']=wq(xc,w,pp/100)
    for X in (1e-9,3e-9,1e-8,3e-8):
        for f in (0.5,1,2): r[f'P_Mdyn_X{X:g}_f{f:g}']=float(w[Mc1/X<=f*MDYN[reg]].sum())
        Ls=L1/X; r[f'Lpath_X{X:g}_50']=wq(Ls,w,.5); r[f'Lpath_X{X:g}_16']=wq(Ls,w,.16); r[f'Lpath_X{X:g}_84']=wq(Ls,w,.84)
        r[f'Lfrac_X{X:g}_50']=wq(Ls/(sz*PC),w,.5)       # path / core FWHM
    if two:
        we=10**P['lwe']; NHe=10**P['NHe']; Me1=we*NHe*mu*ABEAM*(PC*pc)**2/Msun
        ff=Mc1/(Mc1+Me1)
        for pp in (16,50,84): r[f'corefrac_{pp}']=wq(ff,w,pp/100); r[f'Menv1_{pp}']=wq(Me1,w,pp/100)
        r['Te_50']=wq(10**P['lTc']*0+ (10**p['lTe']),w,.5)
    rows.append(r)
D=pd.DataFrame(rows); D.to_csv('/workspace/out/q51/q51_runs.csv',index=False)
# per-config summary: seeds pooled (equal run weights via mean of per-seed values; also seed range)
key=['dset','model','reg','geom','h13','dat']
num=[c for c in D.columns if c not in key+['seed']]
A=D.groupby(key)[num].mean(); A['nseeds']=D.groupby(key).size()
for c in ('logZ','Tc_50','Xcrit_50','P_Mdyn_X3e-09_f1','corefrac_50','Lpath_X1e-08_50'):
    if c in D: A[c+'_min']=D.groupby(key)[c].min(); A[c+'_max']=D.groupby(key)[c].max()
A['logZ_sd']=D.groupby(key)['logZ'].std()
A.to_csv('/workspace/out/q51/q51_summary.csv'); pd.set_option('display.width',250,'display.max_columns',60,'display.max_rows',200)
print(A[['nseeds','logZ','logZ_sd','chi2min','Tc_16','Tc_50','Tc_84','P_T50','tau32_50','P_tau10','P_n6','logX_50','Xcrit_16','Xcrit_50','Xcrit_84']].to_string(float_format=lambda v:'%.3g'%v))
_=(0)#.to_string())
print(A[['P_Mdyn_X1e-09_f1','P_Mdyn_X3e-09_f1','P_Mdyn_X1e-08_f1','P_Mdyn_X3e-09_f0.5','P_Mdyn_X3e-09_f2','Lpath_X1e-08_50','Lpath_X3e-09_50','Lfrac_X1e-08_50','Lfrac_X3e-09_50','corefrac_50']].round(3).to_string())
# H13CN: dlnZ(with limit - without) and survival
print('H13CN survival Z(h1)/Z(h0):')
for k,gp in D.groupby(['dset','model','reg','geom','dat']):
    z0=gp[gp.h13==0].groupby('seed').logZ.mean(); z1=gp[gp.h13==1].groupby('seed').logZ.mean(); s=np.exp(z1-z0).dropna()
    if len(s): print(k,'survival per seed',np.round(s.values,2),'mean %.2f'%s.mean())
# dlnZ(2cT-1z)
print('dlnZ(2cT-1z):')
for k,gp in D.groupby(['dset','reg','geom','h13','dat']):
    a=gp[gp.model=='2cT'].groupby('seed').logZ.mean(); b=gp[gp.model=='1z'].groupby('seed').logZ.mean(); s=(a-b).dropna()
    if len(s): print(k,np.round(s.values,2),'mean %.2f sd %.2f'%(s.mean(),s.std()))
