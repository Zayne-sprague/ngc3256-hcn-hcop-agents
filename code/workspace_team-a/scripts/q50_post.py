# Q50 post-processing: (1) summaries of single-nucleus NS fits (out/q50/runs) vs joint fits; (2) profile chi2(Tc) from NS samples;
# (3) H13CN 4-3 opacity datum applied to samples by importance reweighting: model HCN/H13CN 4-3 intensity ratio (1-exp(-tau))/(1-exp(-tau/X)) >= Rlim,
#     tau = core-grid HCN 4-3 line-centre tau (out/q36 grid_hcn_<reg>_cmb), X = 12C/13C ~ logU(40,100); Rlim from out/h13cn_limits.csv (N: 2.3 [zero CO], 9.5 [log-quad wing], 12 [exp wing]); soft edge 10%.
import glob,json,sys,numpy as np,pandas as pd
from scipy.ndimage import map_coordinates
rs=np.random.default_rng(1)
def wq(x,w,q): o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return np.interp(q,c,x[o])
G={r:np.load(f'/workspace/out/q36/grid_hcn_{r}_cmb.npz') for r in 'NS'}
LTK=np.log10(G['N']['TK']); LN=G['N']['LN']; LNC=G['N']['LNC']
def tau43(r,lT,n,lN):
    c=np.array([np.interp(lT,LTK,np.arange(len(LTK))),np.interp(n,LN,np.arange(len(LN))),np.interp(lN,LNC,np.arange(len(LNC)))])
    t=np.array(G[r]['tau'][...,3],float); t[~(t>0)]=np.nan
    return 10**map_coordinates(np.log10(t),c,order=1,mode='nearest')
def hw(S,w,nm,r,Rlim,X):
    P=lambda k:S[:,nm.index(k)]; t=tau43(r,P(f'lTc_{r}'),P(f'nc_{r}'),P(f'NHc_{r}')); t=np.nan_to_num(t,nan=1e3)
    R=(1-np.exp(-t))/(1-np.exp(-t/X)); L=1/(1+np.exp((np.log(Rlim)-np.log(R))/0.1)); return w*L,t
rows=[]
files=[(f,'single') for f in sorted(glob.glob('/workspace/out/q50/runs/*.npz'))]+[(f,'joint') for f in sorted(glob.glob('/workspace/out/q42/runs2/q42_*_cmb_sh_s*.npz'))]
for f,kind in files:
    d=np.load(f,allow_pickle=True); S,w,nm=d['S'],d['w'],list(d['names']); J=json.load(open(f.replace('.npz','.json'))); key=f.split('/')[-1][:-4]
    for r in 'NS':
        if f'lTc_{r}' not in nm: continue
        T=10**S[:,nm.index(f'lTc_{r}')]; row=dict(file=key,kind=kind,reg=r,mode=J['mode'],seed=J['seed'],logZ=J['logZ'],chi2min=J['chi2min'])
        row['T16'],row['T50'],row['T84']=wq(T,w,[.16,.5,.84]); row['P_T_le30']=float(w[T<=30].sum()); row['P_T_ge40']=float(w[T>=40].sum())
        # profile chi2 vs T (min over samples in log-T bins)
        b=np.linspace(1,2.176,10); ib=np.digitize(np.log10(T),b); ll=d['logl']; prof=[(-2*ll[ib==k].max()) for k in range(1,10) if (ib==k).any()]
        row['dchi2_prof_max']=max(prof)-min(prof); row['dchi2_prof_10_150']=';'.join(f'{p-min(prof):.1f}' for p in prof)
        if f'NHc_{r}' in nm and kind=='single' or r=='N':
            for Rl in (2.3,9.5,12.0):
                ww=0; 
                for X in (40,60,100):
                    w2,t=hw(S,w,nm,r,Rl,X); ww=ww+w2/3
                ww=ww/ww.sum(); Tn=wq(T,ww,[.16,.5,.84]); row[f'H13CN_R{Rl}_T16_50_84']='%.0f/%.0f/%.0f'%tuple(Tn); row[f'H13CN_R{Rl}_Zfac']=float(np.mean([hw(S,w,nm,r,Rl,X)[0].sum() for X in (40,60,100)])); row[f'H13CN_R{Rl}_P_le30']=float(ww[T<=30].sum())
            row['tau43_med']=wq(t,w,[.5])[0]; row['tau43_16_84']='%.1f-%.1f'%tuple(wq(t,w,[.16,.84]))
        rows.append(row)
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q50/q50_summary.csv',index=False)
pd.set_option('display.width',300); pd.set_option('display.max_columns',50); pd.set_option('display.max_colwidth',60)
print(df.round(2).to_string())
