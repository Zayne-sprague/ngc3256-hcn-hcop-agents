# Q2 sensitivity runs (reviewer fixes 1,2,4): usage  q2_sens.py <grid suffix: ''|_tbg5|_tbg10|_tbg20> <Tmin K>  -> out/q2/sens/sens_<suf>_T<Tmin>.csv
# Per region (N,S) and per class (A, B, F2, pooled equal-weight) over the same 8-datum vector as q4c_joint ('all' mode; no wing/core factor, which added no constraint):
#   chi2_min and best-fit, 16/50/84 % of T, log n(HCN), log N(HCN), log X, tau32(HCN), Tex ratio; marginal posterior P(T) per grid T; profile chi2(T) (min over n,N,kappa within class).
# Tmin = hard lower prior edge on Tkin (grid goes to 10 K).  Grids from q4c_grid.py (Tbg 2.73 K) or q2_grid_tbg.py (Tbg 5,10,20 K).
import sys,os
suf=sys.argv[1] if sys.argv[1]!='none' else ''; tmin=float(sys.argv[2]); os.environ['Q2_GRIDSUF']=suf; os.environ['Q2_TMIN']=str(tmin)
sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q4c_joint as J
os.makedirs('/workspace/out/q2/sens',exist_ok=True)
POOL={'A':'A both thick, kappa=1','B':'B both thick, same Tex (0.8-1.25), kappa free','F2':'F2 separate n, kappa=1, both thick'}
rows=[];prof=[]
for c in ('N','S'):
    d,use,chi,cls,aux=J.cube(c,'all','env'); sh=aux['sh']; B=lambda a: np.broadcast_to(a,sh)
    q={'T':B(J.TK[:,None,None,None,None]),'logn_H':B(J.LN[None,:,None,None,None]),'logNH':B(J.LNC[None,None,None,:,None]),'logNO':B(J.LNC[None,None,None,None,:]),
       'tau32H':B(aux['tH']),'TexRatio':B(aux['xH'])/B(aux['xO'])}
    q['logX']=q['logNH']-q['logNO']
    gbest=min(float(np.where(cls[n][0],chi[ik],1e9).min()) for n in POOL.values() for ik in cls[n][1])
    allw={};
    for key,nm in POOL.items():
        mask,iks=cls[nm]; W=[];V={k:[] for k in q}
        for ik in iks:
            cc=np.where(mask,chi[ik],1e9); w=(np.exp(-0.5*(cc-gbest))).ravel()/len(iks); s=w>1e-6; W.append(w[s])
            for k,a in q.items(): V[k].append(a.ravel()[s])
            # profile chi2(T)
            for iT,T in enumerate(J.TK): prof.append(dict(region=c,cls=key,T=T,chi2=float(cc[iT].min()),grid=suf,Tmin=tmin,_k=ik))
        allw[key]=(np.concatenate(W),{k:np.concatenate(v) for k,v in V.items()})
        w,v=allw[key]; cmin=min(float(np.where(mask,chi[ik],1e9).min()) for ik in iks)
        r=dict(region=c,cls=key,grid=suf,Tmin=tmin,chi2min=cmin,ndata=int(use.sum()))
        if w.sum()>0:
            for k in ('T','logn_H','logNH','logX','tau32H','TexRatio'):
                for p in (16,50,84): r[f'{k}_{p}']=J.wq(v[k],w,p/100)
            r['P(T<=40)']=float(w[v['T']<=40].sum()/w.sum()); r['P(T>=80)']=float(w[v['T']>=80].sum()/w.sum()); r['P(TexRatio<0.8)']=float(w[v['TexRatio']<0.8].sum()/w.sum())
            r['T_lo95']=J.wq(v['T'],w,.05); r['T_hi95']=J.wq(v['T'],w,.95)   # one-sided 95% limits
            r['wclass']=float(w.sum())   # relative evidence weight (flat priors, identical grid) for model comparison between classes
            for T in J.TK: r[f'pT{int(T)}']=float(w[v['T']==T].sum()/w.sum())
        rows.append(r)
    # pooled
    w=np.concatenate([allw[k][0] for k in POOL])/3; v={k:np.concatenate([allw[x][1][k] for x in POOL]) for k in q}
    r=dict(region=c,cls='pooled',grid=suf,Tmin=tmin,ndata=int(use.sum()))
    for k in ('T','logn_H','logNH','logX','tau32H','TexRatio'):
        for p in (16,50,84): r[f'{k}_{p}']=J.wq(v[k],w,p/100)
    r['P(T<=40)']=float(w[v['T']<=40].sum()/w.sum()); r['P(T>=80)']=float(w[v['T']>=80].sum()/w.sum()); r['P(TexRatio<0.8)']=float(w[v['TexRatio']<0.8].sum()/w.sum())
    r['T_lo95']=J.wq(v['T'],w,.05); r['T_hi95']=J.wq(v['T'],w,.95); r['chi2min']=min(x['chi2min'] for x in rows if x['region']==c and x['cls'] in POOL)
    for T in J.TK: r[f'pT{int(T)}']=float(w[v['T']==T].sum()/w.sum())
    rows.append(r)
    # branch best chi2 (Tex ratio <0.8 vs 0.8-1.25) for class A
    for key in ('A','F2'):
        mask,iks=cls[POOL[key]]; cc=np.where(mask,chi[iks[0]],1e9)
        for bn,bm in (('low','<0.8'),('eq','0.8-1.25')):
            xr=q['TexRatio']; m=(xr<0.8) if bn=='low' else ((xr>=0.8)&(xr<=1.25)); cb=np.where(m,cc,1e9); i=np.unravel_index(np.argmin(cb),sh)
            rows.append(dict(region=c,cls=key+' branch '+bm,grid=suf,Tmin=tmin,chi2min=float(cb.min()),T_50=float(q['T'][i]),logn_H_50=float(q['logn_H'][i]),logX_50=float(q['logX'][i]),TexRatio_50=float(q['TexRatio'][i])))
tag=f"{suf or 'tbg2.7'}_T{int(tmin)}"
pd.DataFrame(rows).round(4).to_csv(f'/workspace/out/q2/sens/sens_{tag}.csv',index=False)
pd.DataFrame(prof).drop(columns='_k').round(3).to_csv(f'/workspace/out/q2/sens/prof_{tag}.csv',index=False); print(tag,'done')
