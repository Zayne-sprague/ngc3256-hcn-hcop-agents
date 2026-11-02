import numpy as np,pandas as pd
d=np.load('/workspace/out/q39/grid_fine.npz');LNC=d['LNC']
H=pd.read_csv('/workspace/lit/harada_table13_parsed.csv');H['sp']=H.species_transition.str.replace(' ','')
def get(p,s):
    r=H[(H.pos==p)&(H.sp==s)].iloc[0];return r.I_Kkms,r.eI_Kkms
keys={'hcn':['HCN1−0','HCN3−2'],'hcop':['HCO+1−0','HCO+3−2'],'cs':['CS2−1','CS5−4']}
NF=np.arange(11.5,19.001,0.05)
from scipy.interpolate import CubicSpline
def model(sp):
    a=d[sp][...,:2,0];a=np.nan_to_num(a,nan=1e-3);l=np.log(1.0645*150*np.maximum(a,1e-3));return CubicSpline(LNC,l,axis=2)(NF)  # T,n,N,2
M={s:model(s) for s in keys}
def run(p,fixed_eta=None):
    out={}
    best=1e9;res=[]
    chi_sp={}
    for s in keys:
        y=np.log([get(p,k)[0] for k in keys[s]]);sg=np.sqrt(np.array([get(p,k)[1]/get(p,k)[0] for k in keys[s]])**2+0.05**2)
        etas=np.arange(-3,0.001,0.25)
        # chi2[T,n,eta,N]
        c=(((M[s][:,:,None,:,:]+etas[None,None,:,None,None]-y)/sg)**2).sum(-1)
        chi_sp[s]=c
    tot=chi_sp['hcn'].min(-1)+chi_sp['hcop'].min(-1)+chi_sp['cs'].min(-1)  # T,n,eta
    i=np.unravel_index(tot.argmin(),tot.shape);print(p,'min chi2 (no band scale, 6 lines)',tot.min(),i)
    # profile of A=log N_hcn - N_hcop
    prof={}
    for ih in range(len(NF)):pass
    ch=chi_sp['hcn'];co=chi_sp['hcop'];cc=chi_sp['cs'].min(-1)
    base=cc  # T,n,eta
    lag=np.arange(-20,21)
    pr=np.full(len(lag),1e9)
    for k,L in enumerate(lag):
        # N_hcn index = j+L
        if L>=0: a=ch[...,L:];b=co[...,:len(NF)-L]
        else: a=ch[...,:L];b=co[...,-L:]
        t=(a+b).min(-1)+base;pr[k]=t.min()
    return lag*0.05,pr-pr.min()
r={}
for p in 'NSC':r[p]=run(p)
lag,pn=r['N'];_,ps=r['S']
for L,a,b in zip(lag,pn,ps):print(round(L,2),round(a,2),round(b,2))
# N/S contrast profile: min over Ls of pn(LN)+ps(LS) with LN-LS=D
D=np.arange(-20,21)
for dd in D:
    v=min(pn[i]+ps[j] for i in range(41) for j in range(41) if i-j==dd) if True else 0
    print('D',round(dd*0.05,2),round(v,2))
