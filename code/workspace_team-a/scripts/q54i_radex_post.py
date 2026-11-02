"""Q54 r3 fix2: joint non-LTE (q36_rt escape-prob, CMB bg, 150 km/s) fit of HCN+HNC+CN line ratios at N and S (2.5" matched beam, Tb integrals).
Shared (Tkin,nH2), separate column per species, 6 observed ratios (HCN 3-2/1-0, 4-3/3-2, HNC 3-2/1-0, CN 3-2/1-0, HNC/HCN 1-0, CN/HCN 1-0), lognormal errors incl. systematics.
Prior variants: base (n 1e4-10^7.5, T 10-100 K); +tau10(HCN) in 3-12 (Q47 H13CN). Output: out/q54/radex_post.json, radex_post.csv.
usage: q54i_radex_post.py"""
import numpy as np, json
Q='/workspace/out/'
g={'hcn':np.load(Q+'q36/grid_hcn_N_cmb.npz'),'hnc':np.load(Q+'q54/grid_hnc_cmb.npz'),'cn':np.load(Q+'q54/grid_cn_cmb.npz')}
TK,LN,LNC=g['hnc']['TK'],g['hnc']['LN'],g['hnc']['LNC']
# line TR arrays [T,n,N,line]
H=g['hcn']['TR']; HN=g['hnc']['TR']; C=g['cn']['TR']; Ht=g['hcn']['tau']
# observables per (T,n,N) for each species
hcn10,hcn32,hcn43=H[...,0],H[...,2],H[...,3]
hnc10,hnc32=HN[...,0],HN[...,2]
cn10=C[...,0]+C[...,1]; cn32=C[...,3]+C[...,4]
OBS={'N':dict(h32_10=(0.285,.10),h43_32=(0.535,.12),n32_10=(0.31,.14),c32_10=(0.155,.10),hnc_hcn=(0.48,.15),cn_hcn=(1.72,.10)),
     'S':dict(h32_10=(0.235,.10),h43_32=(0.32,.12),n32_10=(0.29,.14),c32_10=(0.103,.10),hnc_hcn=(0.39,.15),cn_hcn=(1.37,.10))}
def lg(x): return np.log(np.clip(x,1e-12,None))
def chi(mod,o): return ((lg(mod)-np.log(o[0]))/o[1])**2
tt=np.ones((len(TK),len(LN))); 
out={}; rows=[]
nsel=(LN>=4)&(LN<=7.5); tsel=TK<=100
for reg in 'NS':
    o=OBS[reg]
    # chi2 own-species terms [T,n,N]
    cH=chi(hcn32/hcn10,o['h32_10'])+chi(hcn43/hcn32,o['h43_32'])
    cN=chi(hnc32/hnc10,o['n32_10']); cC=chi(cn32/cn10,o['c32_10'])
    for var,xe in (('base',.15),('tau',.15),('ff2x',.69),('ff2x_tau',.69),('ff1.4x',.35)):
        # full 5D: [T,n,NH,NN,NC]
        A=cH[:,:,:,None,None]+cN[:,:,None,:,None]+cC[:,:,None,None,:]
        A=A+chi(hnc10[:,:,None,:,None]/hcn10[:,:,:,None,None],(o['hnc_hcn'][0],np.hypot(o['hnc_hcn'][1],xe)))+chi(cn10[:,:,None,None,:]/hcn10[:,:,:,None,None],(o['cn_hcn'][0],np.hypot(o['cn_hcn'][1],xe)))
        A=np.where(np.isfinite(A),A,1e9)
        m=np.zeros(A.shape,bool); m[np.ix_(tsel,nsel)]=True
        if var.endswith('tau'):
            ok=np.isfinite(Ht[...,0])&(Ht[...,0]>=3)&(Ht[...,0]<=12); m&=ok[:,:,:,None,None]
        W=np.where(m&(A<1e8),np.exp(-0.5*(A-A[m].min())),0.); W/=W.sum()
        r=dict(chi2min=float(A[m].min()))
        def q(vals,w):
            i=np.argsort(vals); c=np.cumsum(w[i]); c/=c[-1]; return [float(np.interp(p,c,vals[i])) for p in (.16,.5,.84)]
        ax={'T':(0,TK),'logn':(1,LN),'logN_HCN':(2,LNC),'logN_HNC':(3,LNC),'logN_CN':(4,LNC)}
        for k,(a,v) in ax.items():
            w=W.sum(axis=tuple(i for i in range(5) if i!=a)); r[k]=q(v,w)
        dH=LNC[None,None,:,None,None]*np.ones(A.shape); 
        d1=(LNC[None,None,None,:,None]-LNC[None,None,:,None,None]+0*A).ravel(); d2=(LNC[None,None,None,None,:]-LNC[None,None,:,None,None]+0*A).ravel()
        r['dex_HNC/HCN']=q(d1,W.ravel()); r['dex_CN/HCN']=q(d2,W.ravel())
        bins=np.arange(-2.125,2.2,.25); r['pmf_HNC']=np.histogram(d1,bins,weights=W.ravel())[0].tolist(); r['pmf_CN']=np.histogram(d2,bins,weights=W.ravel())[0].tolist(); r['pmf_bins']=bins.tolist()
        tau10=np.broadcast_to(Ht[...,0][:,:,:,None,None],A.shape).ravel(); r['tau10_HCN']=q(tau10,W.ravel())
        r['P(tau10_HCN>1)']=float(W.ravel()[tau10>1].sum())
        # model-predicted ratios at the best fit
        i=np.unravel_index(np.argmin(np.where(m,A,1e9)),A.shape); r['best']=dict(T=float(TK[i[0]]),logn=float(LN[i[1]]),logNH=float(LNC[i[2]]),logNN=float(LNC[i[3]]),logNC=float(LNC[i[4]]))
        out[f'{reg}_{var}']=r
        print(reg,var,'dN',r['dex_HNC/HCN'],r['dex_CN/HCN'],json.dumps({k:(np.round(v,2).tolist() if isinstance(v,list) else v) for k,v in r.items() if k!='best'}),r['best'],flush=True)
cen=0.5*(np.array(out['N_base']['pmf_bins'][1:])+np.array(out['N_base']['pmf_bins'][:-1]))
for var in ('base','ff2x','ff1.4x'):
    for k in ('HNC','CN'):
        pN=np.array(out['N_'+var]['pmf_'+k]);pS=np.array(out['S_'+var]['pmf_'+k]); D=np.subtract.outer(cen,cen) # S-N : axis0 S
        dd=(cen[:,None]-cen[None,:]).ravel(); w=(pS[:,None]*pN[None,:]).ravel(); i=np.argsort(dd); c=np.cumsum(w[i]); c/=c[-1]
        out[var+'_delta_S-N_'+k]=[float(np.interp(p,c,dd[i])) for p in (.16,.5,.84)]; out[var+'_P(S<N)_'+k]=float(w[dd<-0.01].sum()/w.sum()+0.5*w[np.abs(dd)<=0.01].sum()/w.sum())
        print(var,k,'delta log(N_X/N_HCN) S-N',np.round(out[var+'_delta_S-N_'+k],2),'P(S<N)',round(out[var+'_P(S<N)_'+k],2))
json.dump(out,open(Q+'q54/radex_post.json','w'),indent=1)
