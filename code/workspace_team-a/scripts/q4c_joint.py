# Q4c-3: joint non-LTE (RADEX, out/q4c/grid4_*.npz, dv=150 km/s; N/dv scaling for other dv) comparison of excitation / opacity / filling / abundance / separate-density interpretations of the
# low HCN/HCO+ at S, now including (i) the published J=1-0 data (Harada+2018 Table 13, 1.7"), (ii) an N-nucleus fit, (iii) a class with SEPARATE densities for HCN and HCO+ (common Tkin),
# (iv) model selection with proper parameter counts (BIC) and marginal likelihood (flat priors in log T,n,N on the RADEX grid; kappa log-uniform 0.03-1), (v) variants without the Band-5 2-1 ratios.
# Data vector (per nucleus):  d0 HCN 3-2/2-1, d1 HCN 4-3/3-2, d2 HCO+ 3-2/2-1, d3 HCO+ 4-3/3-2 (Q1b deblended component ladders, Tb; 14% scale + deblend spread),
#   d4 R32 = HCN/HCO+ 3-2 (1.7" integrated), d5 R10 = HCN/HCO+ 1-0 (Harada 1.7"), d6 HCN 3-2/1-0, d7 HCO+ 3-2/1-0 (mean of Harada's own 3-2 and ours at 1.7", error = half difference (+) 5%),
#   nuisance g = 3-2 vs 1-0 flux-scale factor (common to d6,d7; lognormal sigma 0.10).  Penalty: eta=Tint(3-2, beam-dilution fit)/T_R(3-2) must be <=1 (20% soft edge).
# Modes: 'pred' = fit d0-d4 only, then predict d5-d7 (independent test); 'all' = d0-d7; 'all_no21' = all without d0,d2.   Output out/q4c/joint_*.csv, out/q4c/joint_post_*.npz
import sys,warnings,itertools; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
O='/workspace/out/q4c/'
import os
GSUF=os.environ.get('Q2_GRIDSUF','')       # '' = Tbg 2.73 K grids; '_tbg10','_tbg20' = Q2 background-sensitivity grids (q2_grid_tbg.py)
G={m:np.load(O+f'grid4_{m}{GSUF}.npz') for m in ('hcn','hco+')}
TK,LN,LNC=G['hcn']['TK'],G['hcn']['LN'],G['hcn']['LNC']
KAP=np.round(10**np.linspace(-1.5,0,13),4); GG=np.exp(np.linspace(-0.3,0.3,9)); SG=0.10
L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv'); DIL=pd.read_csv('/workspace/out/q4b/dilution.csv'); LIT=pd.read_csv(O+'lit_vs_ours.csv')
def lad(c,sp,q): r=L[(L.comp==c)&(L.species==sp)&(L.quantity==q)].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
def lit(reg,q): return LIT[(LIT.region==reg)&(LIT.quantity==q)].iloc[0].value
def data(c):
    d=[]
    for sp in ('HCN','HCO+'):
        for q in ('3-2 [Xb4f]/2-1 (Tb)','4-3/3-2 [Xb4f 3-2] (Tb)'):
            v,s=lad(c,sp,q); d.append((v,np.hypot(0.14*v,s)))
    d=[d[0],d[1],d[2],d[3]]   # HCN32/21, HCN43/32, HCO32/21, HCO43/32
    r32=0.5*(lit(c,'H R32 (K km/s)')+lit(c,'ours R32 (K km/s, 1.7")')); d.append((r32,0.07*r32))
    r10=lit(c,'H R10 (K km/s)'); d.append((r10,0.05*r10))
    for sp in ('HCN','HCO+'):
        a=lit(c,f'H {sp} 3-2/1-0 (K km/s)'); b=lit(c,f'ours3-2/H1-0 {sp} (K km/s)'); m=0.5*(a+b); d.append((m,np.hypot(0.5*abs(a-b),0.05*m)))
    return np.array(d)
def prep(rmask='physT30'):
    TH,TO=G['hcn']['TR'],G['hco+']['TR']          # (T,n,N,line 0..3 = 1-0,2-1,3-2,4-3)
    return TH,TO
TH,TO=G['hcn']['TR'],G['hco+']['TR']; AH,AO=G['hcn']['tau'],G['hco+']['tau']; XH,XO=G['hcn']['Tex'],G['hco+']['Tex']
nT,nn,nN=len(TK),len(LN),len(LNC)
# 5-D cube axes: T, nH, nO, NH, NO
def bH(a): return a[:,:,None,:,None]      # HCN arrays (T,n,N) -> (T,nH,1,NH,1)
def bO(a): return a[:,None,:,None,:]      # HCO+ arrays (T,n,N) -> (T,1,nO,1,NO)
def cube(c,mode,prior='physT30'):
    """returns dict of class -> (chi2eff array (kappa,T,nH,nO,NH,NO)), mask arrays, predictions function inputs"""
    d=data(c); use=np.ones(8,bool)
    if mode=='pred': use[5:]=False
    if mode=='all_no21': use[[0,2]]=False
    if mode=='no_ladder_21_and_10': use[[0,2,5,6,7]]=False
    P=[bH(TH[...,2]/TH[...,1]),bH(TH[...,3]/TH[...,2]),bO(TO[...,2]/TO[...,1]),bO(TO[...,3]/TO[...,2])]
    rR32=bH(TH[...,2])/bO(TO[...,2]); rR10=bH(TH[...,0])/bO(TO[...,0]); rH=bH(TH[...,2]/TH[...,0]); rO=bO(TO[...,2]/TO[...,0])
    sh=(nT,nn,nn,nN,nN)
    tT=np.arange(nT)[:,None,None,None,None]; iH=np.arange(nn)[None,:,None,None,None]; iO=np.arange(nn)[None,None,:,None,None]
    dl=DIL[DIL.region==c].set_index('line'); TiH=dl.loc['HCN 3-2','Tint_K']; TiO=dl.loc['HCO+ 3-2','Tint_K']
    etaH=TiH/bH(TH[...,2]); etaO=TiO/bO(TO[...,2])
    pen=(np.clip(etaO-1,0,None)/0.2)**2+(np.clip(etaH-1,0,None)/0.2)**2
    ok=np.isfinite(sum(np.broadcast_to(p,sh) for p in P+[rR32,rR10,rH,rO]))&np.isfinite(pen)
    # physical prior: T>=30 K and N<=n+12.5 for both
    NHg=LNC[None,None,None,:,None]; NOg=LNC[None,None,None,None,:]; nHg=LN[None,:,None,None,None]; nOg=LN[None,None,:,None,None]
    pm=(NHg<=nHg+12.5)&(NOg<=nOg+12.5)&(TK[:,None,None,None,None]>=(30 if prior=='physT30' else float(os.environ.get('Q2_TMIN',0)) if prior=='env' else 0))
    chi=np.zeros((len(KAP),)+sh,np.float32)
    for ik,k in enumerate(KAP):
        lg=[]
        for g in GG:
            pr=P+[k*rR32,k*rR10,g*k*0+g*rH,g*rO]
            # kappa multiplies ratios between HCN and HCO+ (d4,d5); d6,d7 are single-molecule ratios (kappa cancels)
            pr=P+[k*rR32,k*rR10,g*rH,g*rO]
            c2=sum(((np.broadcast_to(pr[i],sh)-d[i,0])/d[i,1])**2 for i in range(8) if use[i])+(np.log(g)/SG)**2+pen
            lg.append(-0.5*c2)
        lg=np.array(lg); m=lg.max(0); chi[ik]=(-2*(m+np.log(np.exp(lg-m).sum(0)/len(GG)))).astype(np.float32)
    chi=np.where(ok&pm,chi,np.float32(1e9))
    # class masks
    tH=bH(AH[...,2]); tO=bO(AO[...,2]); xH=bH(XH[...,2]); xO=bO(XO[...,2]); diag=(iH==iO)
    B=lambda a: np.broadcast_to(a,sh)
    cls={'A both thick, kappa=1':(B(diag&(tH>1)&(tO>1)),[len(KAP)-1]),
         'B both thick, same Tex (0.8-1.25), kappa free':(B(diag&(tH>1)&(tO>1)&(xH/xO>0.8)&(xH/xO<1.25)),list(range(len(KAP)))),
         'C HCN thin, HCO+ thick, kappa=1':(B(diag&(tH<1)&(tO>1)),[len(KAP)-1]),
         'D both thin, kappa=1':(B(diag&(tH<1)&(tO<1)),[len(KAP)-1]),
         'F separate n (common Tkin), both thick, kappa free':(B((tH>1)&(tO>1)),list(range(len(KAP)))),
         'F2 separate n, kappa=1, both thick':(B((tH>1)&(tO>1)),[len(KAP)-1])}
    return d,use,chi,cls,dict(P=P,rR32=rR32,rR10=rR10,rH=rH,rO=rO,tH=tH,tO=tO,xH=xH,xO=xO,etaH=etaH,etaO=etaO,sh=sh,pm=pm)
def wq(x,w,q): i=np.argsort(x); c=np.cumsum(w[i]); c/=c[-1]; return np.interp(q,c,x[i])
def run(c,mode,prior='physT30',save=True):
    d,use,chi,cls,aux=cube(c,mode,prior); sh=aux['sh']; B=lambda a: np.broadcast_to(a,sh)
    nd=int(use.sum()); rows=[]
    TT=np.broadcast_to(TK[:,None,None,None,None],sh); NHh=np.broadcast_to(LNC[None,None,None,:,None],sh); NOo=np.broadcast_to(LNC[None,None,None,None,:],sh)
    nHh=np.broadcast_to(LN[None,:,None,None,None],sh); nOo=np.broadcast_to(LN[None,None,:,None,None],sh)
    Zs={}
    for nm,(mask,iks) in cls.items():
        best=1e30; W=[]; 
        for ik in iks:
            cc=np.where(mask,chi[ik],np.float32(1e9)); best=min(best,float(cc.min()))
        npar=4+(1 if nm[0] in 'BF' and 'kappa=1' not in nm else 0)+(1 if nm[0]=='F' else 0)
        nprior=(nT*nn*nN*nN) if nm[0]!='F' else (nT*nn*nn*nN*nN)
        # marginal likelihood relative to exp(-best/2): sum over prior cells (and kappa grid / count)
        Z=sum(np.exp(-0.5*(np.where(mask,chi[ik],1e9)-best)).sum() for ik in iks)/len(iks)
        Zs[nm]=(Z*np.exp(-0.5*best),best,npar,nprior)
        # posterior samples summary
        w=np.concatenate([np.exp(-0.5*(np.where(mask,chi[ik],1e9)-best)).ravel() for ik in iks]); 
        def cat(a): return np.concatenate([np.broadcast_to(a,sh).ravel() for ik in iks])
        kk=np.concatenate([np.full(int(np.prod(sh)),KAP[ik]) for ik in iks])
        s=w>1e-4*w.max(); w=w[s]
        info=dict(mode=mode,region=c,prior=prior,cls=nm,chi2_min=best,ndata=nd,npar=npar+1,BIC=best+(npar+1)*np.log(nd),logZ=np.log(Z)-np.log(nprior)+(-0.5*best))
        vals={'T':TT,'logn_H':nHh,'logn_O':nOo,'logNH':NHh,'logNO':NOo,'tau32H':aux['tH'],'tau32O':aux['tO'],'TexH':aux['xH'],'TexO':aux['xO'],'etaH':aux['etaH'],'etaO':aux['etaO']}
        for key,a in vals.items():
            v=cat(a)[s]; 
            for q,qq in ((16,.16),(50,.5),(84,.84)): info[f'{key}_{q}']=wq(v,w,qq)
            info[key+'_lo95']=wq(v,w,.025); info[key+'_hi95']=wq(v,w,.975)
        v=(NHh-NOo); v=cat(v)[s]; info['logNratio_50']=wq(v,w,.5); info['logNratio_lo95']=wq(v,w,.025); info['logNratio_hi95']=wq(v,w,.975)
        info['kappa_50']=wq(kk[s],w,.5)
        # predictions of ALL 8 data (without g) and 1-0 quantities: with kappa
        def prd(a,kfac): return (cat(a)[s]*(kk[s] if kfac else 1.))
        pr={'R10':prd(aux['rR10'],1),'R32':prd(aux['rR32'],1),'HCN32/10':prd(aux['rH'],0),'HCO32/10':prd(aux['rO'],0),'HCN32/21':cat(aux['P'][0])[s],'HCN43/32':cat(aux['P'][1])[s],'HCO32/21':cat(aux['P'][2])[s],'HCO43/32':cat(aux['P'][3])[s]}
        for key,v in pr.items():
            for q,qq in ((16,.16),(50,.5),(84,.84)): info[f'pred_{key}_{q}']=wq(v,w,qq)
            info[f'pred_{key}_lo95']=wq(v,w,.025); info[f'pred_{key}_hi95']=wq(v,w,.975)
        # best-fit point
        ik=min(iks,key=lambda i:float(np.where(mask,chi[i],1e9).min())); cc=np.where(mask,chi[ik],1e9); ii=np.unravel_index(np.argmin(cc),sh)
        for key,a in vals.items(): info['best_'+key]=float(np.broadcast_to(a,sh)[ii])
        info['best_kappa']=KAP[ik]
        for key,a in (('R10',KAP[ik]*aux['rR10']),('R32',KAP[ik]*aux['rR32']),('HCN32/10',aux['rH']),('HCO32/10',aux['rO'])): info['bestpred_'+key]=float(np.broadcast_to(a,sh)[ii])
        rows.append(info)
    # posterior class probabilities (prior-weighted) among A,B,C,D (same-n families) and among A,B,F
    return rows,d,use

def posterior(c,mode,names,prior='physT30',thr=1e-4):
    """weighted posterior samples (flattened cells) of the given classes (list of names, pooled with equal weight per class)  -> dict of arrays incl. 'w'"""
    d,use,chi,cls,aux=cube(c,mode,prior); sh=aux['sh']; out={}
    keys={'T':np.broadcast_to(TK[:,None,None,None,None],sh),'logn_H':np.broadcast_to(LN[None,:,None,None,None],sh),'logn_O':np.broadcast_to(LN[None,None,:,None,None],sh),
          'logNH':np.broadcast_to(LNC[None,None,None,:,None],sh),'logNO':np.broadcast_to(LNC[None,None,None,None,:],sh),'tau32H':np.broadcast_to(aux['tH'],sh),'tau32O':np.broadcast_to(aux['tO'],sh),
          'TexH':np.broadcast_to(aux['xH'],sh),'TexO':np.broadcast_to(aux['xO'],sh),'etaH':np.broadcast_to(aux['etaH'],sh),'etaO':np.broadcast_to(aux['etaO'],sh)}
    W=[];V={k:[] for k in keys}; V['kappa']=[]
    gbest=min(float(np.where(cls[n][0],chi[ik],1e9).min()) for n in names for ik in cls[n][1])
    for n in names:
        mask,iks=cls[n]
        for ik in iks:
            w=np.exp(-0.5*(np.where(mask,chi[ik],1e9)-gbest)).ravel()/len(iks); s=w>thr
            W.append(w[s]); 
            for k,a in keys.items(): V[k].append(a.ravel()[s])
            V['kappa'].append(np.full(s.sum(),KAP[ik]))
    out={k:np.concatenate(v) for k,v in V.items()}; out['w']=np.concatenate(W); return out
if __name__=='__main__':
    allr=[]
    for c in ('S','N'):
        for mode in ('pred','all','all_no21'):
            rows,d,use=run(c,mode); allr+=rows
            print(c,mode,'data',np.round(d[:,0],3),np.round(d[:,1],3),flush=True)
            for r in rows: print('  %-52s chi2=%7.2f n=%d k=%d BIC=%7.2f lnZ=%7.2f tauH %.1f [%.1f,%.1f] T %.0f n_H %.2f n_O %.2f logNH %.2f | pred R10 %.2f [%.2f,%.2f] H32/10 %.2f O32/10 %.2f'%(r['cls'],r['chi2_min'],r['ndata'],r['npar'],r['BIC'],r['logZ'],r['tau32H_50'],r['tau32H_lo95'],r['tau32H_hi95'],r['T_50'],r['logn_H_50'],r['logn_O_50'],r['logNH_50'],r['pred_R10_50'],r['pred_R10_lo95'],r['pred_R10_hi95'],r['pred_HCN32/10_50'],r['pred_HCO32/10_50']),flush=True)
    pd.DataFrame(allr).round(4).to_csv(O+'joint_summary.csv',index=False)
