# Q47(c): apply the 2.2" 1-0 13C/12C ratios (own extraction, out/q47/ratios10.json) to our N-only nested-sampling posteriors (out/q50/runs, shared-nothing envelope,
# CMB grids) by importance reweighting: model R = TR(13C 1-0)/TR(12C 1-0) of the core component at the sample's (Tk,n,N), N(13C)=N(12C)/X, X ~ LN(45, 0.2 dex) truncated <=100
# (team-b's prior) or logU(30,100); model error 10% in R (LAMDA/H13CO+ file convention) added in quadrature.  H13CO+ uses the HCO+ file at N/X (as team-b).
# Reports T quantiles, P(T<=30), P(T>=40), surviving mass Zfac, for HCN-only, HCO+-only, both, and +4-3 limit (Rlim 9.5 as in q50_post).
import json,glob,sys,numpy as np,pandas as pd
from scipy.ndimage import map_coordinates
sys.path.insert(0,'/workspace/scripts')
def wq(x,w,q): o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return np.interp(q,c,x[o])
G12={m:np.load(f'/workspace/out/q36/grid_{m}_N_cmb.npz') for m in('hcn','hco+')}; G13=np.load('/workspace/out/q47/grid_h13cn_10.npz')
LTK=np.log10(G12['hcn']['TK']); LN=G12['hcn']['LN']; LNC=G12['hcn']['LNC']; LNC13=G13['LNC']
def interp(arr,lT,n,lN,axN):
    c=np.array([np.interp(lT,LTK,np.arange(len(LTK))),np.interp(n,LN,np.arange(len(LN))),np.interp(lN,axN,np.arange(len(axN)))])
    a=np.array(arr,float); a[~((a>0)&(a<1e3))]=np.nan; a=np.log10(a); return 10**map_coordinates(a,c,order=1,mode='nearest')
R=json.load(open('/workspace/out/q47/ratios10.json'))
DATA={'own_A':dict(hcn=R['HCN|N|A'],hcop=R['HCO+|N|A']),'own_T':dict(hcn=R['HCN|N|T'],hcop=R['HCO+|N|T']),'own_P':dict(hcn=R['HCN|N|P'],hcop=R['HCO+|N|P']),
      'teamb_T':dict(hcn=dict(R=0.0207,err=0.0034),hcop=dict(R=0.0231,err=0.0031)),'teamb_A':dict(hcn=dict(R=0.0220,err=0.0046),hcop=dict(R=0.0271,err=0.0049))}
def tau43(S,nm):
    G=np.load('/workspace/out/q36/grid_hcn_N_cmb.npz'); t=np.array(G['tau'][...,3],float); t[~(t>0)]=np.nan
    c=np.array([np.interp(S[:,nm.index('lTc_N')],LTK,np.arange(len(LTK))),np.interp(S[:,nm.index('nc_N')],LN,np.arange(len(LN))),np.interp(S[:,nm.index('NHc_N')],LNC,np.arange(len(LNC)))])
    return np.nan_to_num(10**map_coordinates(np.log10(t),c,order=1,mode='nearest'),nan=1e3)
rows=[]; rs=np.random.default_rng(3)
for f in sorted(glob.glob('/workspace/out/q50/runs/q50_N_*.npz')):
    d=np.load(f,allow_pickle=True); S,w,nm=d['S'],d['w'],list(d['names']); key=f.split('/')[-1][:-4]
    lT,n,NH,NO=[S[:,nm.index(k)] for k in('lTc_N','nc_N','NHc_N','NOc_N')]; T=10**lT
    tr12h=interp(G12['hcn']['TR'][...,0],lT,n,NH,LNC); tr12o=interp(G12['hco+']['TR'][...,0],lT,n,NO,LNC)
    for xprior in ('LN45','logU30-100'):
        Xs=np.exp(rs.normal(np.log(45),0.2*np.log(10),400)); Xs=Xs[Xs<=100] if xprior=='LN45' else 10**rs.uniform(np.log10(30),2,400)
        mods=[]
        for X in Xs:
            rh=interp(G13['TR'],lT,n,NH-np.log10(X),LNC13)/tr12h
            ro=interp(G12['hco+']['TR'][...,0],lT,n,NO-np.log10(X),LNC)/tr12o
            mods.append((rh,ro))
        for dn,dd in DATA.items():
            for use in('hcn','hcop','both'):
                L=np.zeros(len(w))
                for rh,ro in mods:
                    ll=0
                    if use in('hcn','both'): ll=ll-0.5*((rh-dd['hcn']['R'])/np.hypot(dd['hcn'].get('err',0),0.10*rh))**2
                    if use in('hcop','both'): ll=ll-0.5*((ro-dd['hcop']['R'])/np.hypot(dd['hcop'].get('err',0),0.10*ro))**2
                    L+=np.nan_to_num(np.exp(ll),nan=0.0)/len(mods)
                w2=w*L; Z=w2.sum(); w2/=Z; q=wq(T,w2,[.16,.5,.84]); ess=1/np.sum(w2**2); bs=[(T[rs.choice(len(w2),len(w2),p=w2)]<=30).mean() for _ in range(20)]
                rows.append(dict(file=key,fbad=float(np.mean(~np.isfinite(tr12h))),xprior=xprior,data=dn,use=use,Zfac=Z,ESS=ess,P_le30_bootsd=float(np.std(bs)),flag_ESS_lt100=bool(ess<100),T16=q[0],T50=q[1],T84=q[2],P_le30=float(w2[T<=30].sum()),P_ge40=float(w2[T>=40].sum()),P_le30_before=float(w[T<=30].sum())))
df=pd.DataFrame(rows); df.round(3).to_csv('/workspace/out/q47/reweight_1_0.csv',index=False)
df.round(4).to_csv('/workspace/out/q47/reweight_1_0.csv',index=False)
g=df.groupby(['xprior','data','use'])[['Zfac','ESS','P_le30_before','P_le30','P_le30_bootsd','P_ge40','T50']].mean().round(3); pd.set_option('display.width',200); print(g.to_string())

g.to_csv('/workspace/out/q47/reweight_summary_by_script.csv')
print(df.groupby(['use'])[['ESS']].describe().round(1).to_string()); print(df.groupby(['use','data','xprior'])['P_le30'].agg(['min','max','std']).round(3).to_string())
