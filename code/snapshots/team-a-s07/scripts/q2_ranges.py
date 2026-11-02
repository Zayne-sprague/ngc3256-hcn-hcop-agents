# Q2: credible ranges of n(H2), Tkin, N(HCN), N(HCO+), X, tau, Tex, filling from the joint RADEX posterior (machinery of q4c_joint.py:
# LAMDA hcn/hco+ grids out/q4c/grid4_*.npz, dv=150 km/s, Tbg 2.73 K, H2 collider, uniform sphere, flat priors in log T,n,N on the grid, T>=30 K physical prior).
# Adds to q4c: (1) pooled "both thick" posterior (classes A,B,F2 equal weight) with derived quantities (X, pressure, Tex ratio, eta from native Tb peaks);
# (2) explicit test of the S low-n Tex(HCN)<Tex(HCO+) branch vs the equal-Tex/high-n branch, per-datum pulls; (3) velocity-resolved constraint: the HCN/HCO+ 3-2 ratio in the
# line wings relative to the core at the SAME 0.30" beam (N 1.07+-0.13, S 1.06+-0.12, out/q9b/windows.csv) must follow from a column 0.25-1 dex lower (f=0.56..0.1);
# (4) H13CN/CS consistency. Output out/q2/*.csv
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4c_joint as J
O='/workspace/out/q2/'; rng=np.random.default_rng(11)
WING={'N':(1.07,0.13),'S':(1.06,0.12)}          # R_wing/R_core, 0.30" Xb4b, 100-260 km/s stacked / |v|<60 (q9b windows; errors in quadrature)
POOL=['A both thick, kappa=1','B both thick, same Tex (0.8-1.25), kappa free','F2 separate n, kappa=1, both thick']
SHIFTS=[0,1,2,3,4]
def shifted(a,s):                       # a (T,n,N): value at column index i-s along N (lower column), NaN if off-grid
    b=np.full_like(a,np.nan); b[...,s:]=a[...,:a.shape[-1]-s]; return b
def wq(x,w,q): return J.wq(x,w,q)
def analyse(c,wing=True):
    d,use,chi,cls,aux=J.cube(c,'all'); sh=aux['sh']; B=lambda a: np.broadcast_to(a,sh)
    Rcore=J.bH(J.TH[...,2])/J.bO(J.TO[...,2])
    if wing:
        mu,sg=WING[c]; L=np.zeros(sh,np.float32); nv=np.zeros(sh,np.float32)
        for s in SHIFTS[1:]:
            Rs=J.bH(shifted(J.TH[...,2],s))/J.bO(shifted(J.TO[...,2],s)); q=Rs/Rcore
            ok=np.isfinite(B(q)); L+=np.where(ok,np.exp(-0.5*((B(q)-mu)/sg)**2),0).astype(np.float32); nv+=ok
        Lw=np.where(nv>0,L/np.maximum(nv,1),0.)
    else: Lw=np.ones(sh,np.float32)
    keys={'T':J.TK[:,None,None,None,None],'logn_H':J.LN[None,:,None,None,None],'logn_O':J.LN[None,None,:,None,None],'logNH':J.LNC[None,None,None,:,None],'logNO':J.LNC[None,None,None,None,:]}
    der={'logX':B(keys['logNH'])-B(keys['logNO']),'tau32H':B(aux['tH']),'tau32O':B(aux['tO']),'TexH':B(aux['xH']),'TexO':B(aux['xO']),'TexRatio':B(aux['xH'])/B(aux['xO']),
         'etaH':B(aux['etaH']),'etaO':B(aux['etaO']),'logP':B(keys['logn_H'])+np.log10(B(keys['T']))}
    allv={**{k:B(v) for k,v in keys.items()},**der}
    W=[];V={k:[] for k in allv}; V['kappa']=[]; V['cls']=[]; gbest=min(float(np.where(cls[n][0],chi[ik],1e9).min()) for n in POOL for ik in cls[n][1])
    for n in POOL:
        mask,iks=cls[n]
        for ik in iks:
            w=(np.exp(-0.5*(np.where(mask,chi[ik],1e9)-gbest))*Lw).ravel()/len(iks); s=w>1e-5
            W.append(w[s])
            for k,a in allv.items(): V[k].append(a.ravel()[s])
            V['kappa'].append(np.full(s.sum(),J.KAP[ik])); V['cls'].append(np.full(s.sum(),n[0]))
    P={k:np.concatenate(v) for k,v in V.items()}; P['w']=np.concatenate(W); P['w']/=P['w'].sum()
    return P,d
def summ(P,c,tag,rows):
    w=P['w']
    for k in ('T','logn_H','logn_O','logNH','logNO','logX','tau32H','tau32O','TexH','TexO','TexRatio','etaH','etaO','logP','kappa'):
        v=P[k]; rows.append(dict(region=c,variant=tag,quantity=k,p2p5=wq(v,w,.025),p16=wq(v,w,.16),p50=wq(v,w,.5),p84=wq(v,w,.84),p97p5=wq(v,w,.975)))
    # probabilities
    for nm,m in (('P(logX<0, N(HCN)<N(HCO+))',P['logX']<-0.01),('P(logX>=0.5)',P['logX']>=0.5),('P(Tex_HCN/Tex_HCO+<0.8)',P['TexRatio']<0.8),('P(log n_H<4.5)',P['logn_H']<4.5),('P(log n_H>=5)',P['logn_H']>=5),('P(T<=40K)',P['T']<=40),('P(T>=80K)',P['T']>=80),('P(tau32H<3)',P['tau32H']<3)):
        rows.append(dict(region=c,variant=tag,quantity=nm,p50=float(w[m].sum())))
rows=[]; Ps={}
for c in ('N','S'):
    for wing in (False,True):
        P,d=analyse(c,wing); tag='ladder+R32+1-0 (q4c all)' if not wing else '+ wing/core ratio constraint'; summ(P,c,tag,rows); Ps[(c,wing)]=P
        print(c,tag,'done',flush=True)
R=pd.DataFrame(rows); R.round(4).to_csv(O+'q2_ranges.csv',index=False)
pd.set_option('display.width',250); print(R[R.variant.str.startswith('+')].round(3).to_string())
# (2) branch test for S (and N): posterior mass and per-datum pulls of the branches (ladder only data d0..d4 + 1-0)
br=[]
for c in ('S','N'):
    d,use,chi,cls,aux=J.cube(c,'all'); sh=aux['sh']
    for nm in ('A both thick, kappa=1','F2 separate n, kappa=1, both thick'):
        mask,iks=cls[nm]; ik=iks[0]; cc=np.where(mask,chi[ik],1e9)
        xr=np.broadcast_to(aux['xH']/aux['xO'],sh)
        for bn,bm in (('Tex(HCN)/Tex(HCO+)<0.8 (low-n HCN branch)',xr<0.8),('0.8-1.25 (equal Tex)',(xr>=0.8)&(xr<=1.25)),('>1.25',xr>1.25)):
            cb=np.where(bm,cc,1e9); i=np.unravel_index(np.argmin(cb),sh); best=float(cb.min())
            if best>1e8: br.append(dict(region=c,cls=nm,branch=bn,chi2min=np.nan)); continue
            pr=[float(np.broadcast_to(a,sh)[i]) for a in aux['P']]+[float(J.KAP[ik]*np.broadcast_to(aux['rR32'],sh)[i]),float(J.KAP[ik]*np.broadcast_to(aux['rR10'],sh)[i])]
            row=dict(region=c,cls=nm,branch=bn,chi2min=best,Zrel=float(np.exp(-0.5*(cb-cc.min())).sum()/np.exp(-0.5*(cc-cc.min())).sum()),
                     T=float(np.broadcast_to(J.TK[:,None,None,None,None],sh)[i]),logn_H=float(np.broadcast_to(J.LN[None,:,None,None,None],sh)[i]),logn_O=float(np.broadcast_to(J.LN[None,None,:,None,None],sh)[i]),
                     logNH=float(np.broadcast_to(J.LNC[None,None,None,:,None],sh)[i]),logNO=float(np.broadcast_to(J.LNC[None,None,None,None,:],sh)[i]))
            for lab,p,j in zip(('HCN32/21','HCN43/32','HCO32/21','HCO43/32','R32','R10'),pr,range(6)): row['pred_'+lab]=p; row['pull_'+lab]=(p-d[j,0])/d[j,1]
            br.append(row)
B=pd.DataFrame(br); B.round(3).to_csv(O+'q2_branches.csv',index=False); print(B.round(2).to_string())
# (4) H13CN consistency with pooled N/S posteriors: predicted HCN/H13CN (4-3 Tb) vs limits
G13=np.load('/workspace/out/q4c/grid4_h13cn.npz'); G12=J.G['hcn']
# (5) derived physical scales. Assumptions: line width per cloud 50 km/s (RADEX depends on N/dv, so N(cloud)=N(dv=150)/3, q4c phys_dv_check), X(HCN)=N(HCN)/N(H2)=1e-9,1e-8,1e-7 (assumed, not measured),
# 1 pc=3.086e18 cm. Implied coherent path L=N(HCN)/(X n_H2) and N(H2); compared with the CO-based GMC mean density (Brunetti+2022, MNRAS 515, 2928: median Sigma=470 Msun/pc2, R=88 pc;
# 1 Msun/pc2 = 4.6e19 H2 cm-2 incl. He (mean mass 2.72 m_H per H2); n=N(H2)/(2R)).
x=[]
for c in ('N','S'):
    P=Ps[(c,True)]; w=P['w']
    for lx in (-9,-8,-7):
        NHc=P['logNH']-np.log10(3.); NH2=NHc-lx; L=NHc-lx-P['logn_H']-np.log10(3.086e18)
        for nm,v in (('log N(H2) [cm-2]',NH2),('log L [pc]',L),('log Sigma_H2 [Msun/pc2]',NH2-np.log10(4.6e19))):
            x.append(dict(region=c,logX_HCN=lx,quantity=nm,p16=wq(v,w,.16),p50=wq(v,w,.5),p84=wq(v,w,.84)))
n_co=470*4.6e19/(2*88*3.086e18); x.append(dict(region='CO GMC',quantity='n_H2 mean [cm-3] (Sigma=470, R=88 pc; linear, not log)',p50=n_co,logX_HCN=np.nan))
pd.DataFrame(x).round(3).to_csv(O+'q2_scales.csv',index=False); print(pd.DataFrame(x).round(2).to_string())
