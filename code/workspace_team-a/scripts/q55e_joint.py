# Q55e: matched-beam joint RADEX model comparison of N vs S (replaces the 3-datum chi2 of q55_ns_single_param.py as the Q55 test).
# Data per nucleus (Tb ratios HCN/HCO+ in K km/s):
#   set A: 2.5" common circular beam: R10 (2015.1.00993), R32 (Xb4d, fixed window), R43 (2018.1.00493/00223, fixed window, uv-mismatch -> 25% shared nuisance),
#          + 13C/12C 1-0 at 2.2" (Q47: H13CN/HCN, H13CO+/HCO+)
#   set B: core 1.0-1.1" : R32 (Xb4d), R43 (2018.1.00493/00223); no 1-0 at this beam exists (best 1-0 beam 2.2")
# Parameters per nucleus theta=(Tkin, n, N(HCN), X=N(HCN)/N(HCO+), X13=12C/13C). Flat priors on the grid (log n 3-7.5, log N 13-17.5, logX -1..1,
# T in {20,30,50,80,150}, log X13 in {1.5,1.75,2.0}), N(HCN)/n<=10^12.5 (path<=100pc at X(HCN)=1e-8). Calibration of each J ratio: lognormal nuisance shared by N and S,
# marginalised (Gauss-Hermite). Evidence = prior-weighted mean likelihood (marginal likelihood); ln Z differences are ln Bayes factors.
import sys,json,itertools; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
import q40h_radex as R
O='/workspace/out/q55/'
TH=np.concatenate([np.load(O+'radex_low_hcn.npy'),np.load('/workspace/out/q40b/radex_hcn_J4.npy')],axis=2)[...,0]   # T,n,N(10..17.5),J
LO=np.load(O+'radex_low_hco+.npy'); MO=np.load('/workspace/out/q40b/radex_hco+_J4.npy')
TO=np.concatenate([LO[:,:,:8],MO],axis=2)[...,0]  # low 10..11.75, main 12..18.5
LH=np.arange(10,17.51,0.25); LOo=np.arange(10,18.51,0.25); assert TH.shape[2]==len(LH) and TO.shape[2]==len(LOo)
TK,LN,LNH,LX=R.TK,R.LN,R.LNH,np.arange(-1,1.01,0.25); LX13=np.array([1.5,1.75,2.0])
shape=(len(TK),len(LN),len(LNH),len(LX),len(LX13)); AX=['T','n','NH','X','X13']
iT,iN,iH,iX,i13=np.indices(shape)
lnH=LNH[iH]; lnO=lnH-LX[iX]; l13=LX13[i13]
def ix(l,grid): return np.rint((l-grid[0])/0.25).astype(int)
def tb(Tab,lcol,grid,J): return Tab[iT,iN,ix(lcol,grid),J]
M={}
for J,nm in ((0,'10'),(2,'32'),(3,'43')): M['R'+nm]=tb(TH,lnH,LH,J)/tb(TO,lnO,LOo,J)
M['H13_10']=tb(TH,lnH-l13,LH,0)/tb(TH,lnH,LH,0)*1.054; M['O13_10']=tb(TO,lnO-l13,LOo,0)/tb(TO,lnO,LOo,0)*1.054
M['H13_43']=tb(TH,lnH-l13,LH,3)/tb(TH,lnH,LH,3)*1.054
TRH={J:tb(TH,lnH,LH,J) for J in (0,2,3)}; TRO={J:tb(TO,lnO,LOo,J) for J in (0,2,3)}
M['fH10']=M['H13_10']*10**l13; M['fO10']=M['O13_10']*10**l13; M['fH43']=M['H13_43']*10**l13
M['tau32']=tb(np.concatenate([np.load(O+'radex_low_hcn.npy'),np.load('/workspace/out/q40b/radex_hcn_J4.npy')],axis=2)[...,1],lnH,LH,2)
M['tau10H']=tb(np.concatenate([np.load(O+'radex_low_hcn.npy'),np.load('/workspace/out/q40b/radex_hcn_J4.npy')],axis=2)[...,1],lnH,LH,0)
M['Tk']=TK[iT]; M['logn']=LN[iN]
valid=(lnH-LN[iN]<=12.5)&np.isfinite(M['H13_10'])&np.isfinite(M['O13_10'])&np.isfinite(M['H13_43'])&np.isfinite(M['R10'])&np.isfinite(M['R32'])&np.isfinite(M['R43'])&(M['R10']>0)&(M['R32']>0)&(M['R43']>0)
import os
X13MAX=float(os.environ.get('X13MAX','2.0')); valid&=(l13<=X13MAX+1e-9); TAG=os.environ.get('TAG','')
# ---- data -------------------------------------------------------------------------------------------------------------------
C=1.0; fl=0.05   # ratios are already Tb (K km/s) ratios (q4_lib.kb applied in the source scripts): no frequency factor
   # independent 5% floor per ratio (line-profile/deblend/model)
def lnd(v,e): return (np.log(v*C), np.hypot(e/v,fl))
SETS={
 'A_full': dict(R10=('N',0.906,0.010,'S',0.575,0.013), R32=('N',0.698,0.005,'S',0.386,0.006), R43=('N',0.624,0.015,'S',0.174,0.012), sig={'R10':0.07,'R32':0.05,'R43':0.25}, c13=True, abs='2.5'),
 'A_no43':dict(abs='2.5',R10=('N',0.906,0.010,'S',0.575,0.013), R32=('N',0.698,0.005,'S',0.386,0.006), sig={'R10':0.07,'R32':0.05}, c13=True),
 'A_no13':dict(abs='2.5',R10=('N',0.906,0.010,'S',0.575,0.013), R32=('N',0.698,0.005,'S',0.386,0.006), R43=('N',0.624,0.015,'S',0.174,0.012), sig={'R10':0.07,'R32':0.05,'R43':0.25}, c13=False),
 'A_R10R32':dict(abs='2.5',R10=('N',0.906,0.010,'S',0.575,0.013), R32=('N',0.698,0.005,'S',0.386,0.006), sig={'R10':0.07,'R32':0.05}, c13=False),
 'A_43wide':dict(abs='2.5',R10=('N',0.906,0.010,'S',0.575,0.013), R32=('N',0.698,0.005,'S',0.386,0.006), R43=('N',0.624,0.015,'S',0.174,0.012), sig={'R10':0.07,'R32':0.05,'R43':0.6}, c13=True),
 'A_43Xb47':dict(abs='2.5',R10=('N',0.906,0.010,'S',0.575,0.013), R32=('N',0.698,0.005,'S',0.386,0.006), R43=('N',0.515,0.008,'S',0.317,0.011), sig={'R10':0.07,'R32':0.05,'R43':0.10}, c13=True),
 'B_core3':dict(abs='1.0',R10=('N',1.059/C,0.055/C,'S',0.731/C,0.106/C), R32=('N',0.835,0.007,'S',0.402,0.008), R43=('N',0.814,0.019,'S',0.222,0.015), sig={'R10':0.05,'R32':0.05,'R43':0.25}, c13=False),
 'B_core3_no43':dict(abs='1.0',R10=('N',1.059/C,0.055/C,'S',0.731/C,0.106/C), R32=('N',0.835,0.007,'S',0.402,0.008), sig={'R10':0.05,'R32':0.05}, c13=False),
 'B_core3_13C':dict(abs='1.0',R10=('N',1.059/C,0.055/C,'S',0.731/C,0.106/C), R32=('N',0.835,0.007,'S',0.402,0.008), R43=('N',0.814,0.019,'S',0.222,0.015), sig={'R10':0.05,'R32':0.05,'R43':0.25}, c13=True),
 'B_core3_no43_13C':dict(abs='1.0',R10=('N',1.059/C,0.055/C,'S',0.731/C,0.106/C), R32=('N',0.835,0.007,'S',0.402,0.008), sig={'R10':0.05,'R32':0.05}, c13=True),
 'B_core':dict(abs='1.0',R32=('N',0.835,0.007,'S',0.402,0.008), R43=('N',0.814,0.019,'S',0.222,0.015), sig={'R32':0.05,'R43':0.25}, c13=False),
 'B_core32':dict(abs='1.0',R32=('N',0.835,0.007,'S',0.402,0.008), sig={'R32':0.05}, c13=False),
}
# 13C/12C 1-0 at 2.2" (Q47, flux ratios; mean of aperture and template estimators; error = max stat + half spread + 10% model)
# + H13CN 4-3 non-detection (Q47a, 4.6x3.1" beam): alpha=0+-limit/3 (3 sigma <0.10 N, <0.14 S)
C13={'N':{'H13_43':(0.0,0.0333),'H13_10':(0.0211,np.hypot(0.0065,0.0007+0.002)),'O13_10':(0.0234,np.hypot(0.0044,0.004+0.0023))},
     'S':{'H13_43':(0.0,0.0467),'H13_10':(0.0065,np.hypot(0.008,0.0005+0.0007)),'O13_10':(0.0129,np.hypot(0.007,0.0015+0.0013))}}
# absolute integrated Tb (K km/s) in the same beam: eta=I_obs/I_model must be <=1 (beam filling cannot exceed 1); soft edge 25% in ln
ABS={'2.5':{0:{'N':(135.68,149.8),'S':(59.67,103.8)},2:{'N':(40.86,59.73),'S':(14.30,37.92)},3:{'N':(21.49,37.99),'S':(6.52,27.36)}},
     '1.0':{0:{'N':(242.2,231.7),'S':(62.0,85.9)},2:{'N':(110.3,135.1),'S':(33.4,87.3)},3:{'N':(72.2,102.7),'S':(18.0,82.3)}}}
LW=1.0645*150.
def abs_pen(reg,key,Js):
    p=np.zeros(shape)
    for J in Js:
        for obs,TR in zip(ABS[key][J][reg],(TRH[J],TRO[J])):
            x=np.log(obs/np.clip(TR*LW,1e-9,None)); p+=-0.5*(np.clip(x,0,None)/0.25)**2
    return p
nod,wts=hermegauss(7); wts=wts/np.sqrt(2*np.pi)
def region_lik(reg,S):
    Js=[j for j in ('R10','R32','R43') if j in S]
    d={j:lnd(*(S[j][1:3] if reg=='N' else S[j][4:6])) for j in Js}
    combos=list(itertools.product(range(len(nod)),repeat=len(Js)))
    L=np.zeros((len(combos),)+shape); W=np.zeros(len(combos))
    base=np.zeros(shape)
    if S.get('abs'): base+=abs_pen(reg,S['abs'],[{'R10':0,'R32':2,'R43':3}[j] for j in Js])
    if S['c13']:
        for k,(v,e) in C13[reg].items():
            base+=-0.5*((v*1.054-M[k])/(e*1.054))**2   # model Tb ratio vs observed flux ratio->Tb
    for ic,cc in enumerate(combos):
        chi=base.copy(); w=1.0
        for j,k in zip(Js,cc):
            c=S['sig'][j]*nod[k]; w*=wts[k]
            chi+=-0.5*((d[j][0]-np.log(M[j])-c)/d[j][1])**2
        L[ic]=np.where(valid,np.exp(np.where(valid,chi,0)),0.0); W[ic]=w
    return L,W,Js
HYP={'all shared':[], 'T differs':['T'], 'n differs':['n'], 'N(HCN) differs':['NH'], 'HCN/HCO+ column ratio X differs':['X'],
     'T,n differ (excitation)':['T','n'], 'N(HCN),X differ (column+abundance)':['NH','X'], 'T,n,N(HCN) differ':['T','n','NH'], 'all four differ':['T','n','NH','X']}
def run(name,S,share13=True,verbose=True):
    LNr,W,Js=region_lik('N',S); LSr,_,_=region_lik('S',S)
    out={}
    for h,D in HYP.items():
        D13=D+([] if share13 else ['X13'])
        di=[AX.index(a) for a in D13]; si=[i for i in range(5) if i not in di]
        perm=si+di; ns=int(np.prod([shape[i] for i in si])) if si else 1; nd=int(np.prod([shape[i] for i in di])) if di else 1
        Z=0.0; WN=np.zeros((ns,nd)); WS=np.zeros((ns,nd))
        for ic in range(len(W)):
            a=LNr[ic].transpose(perm).reshape(ns,nd); b=LSr[ic].transpose(perm).reshape(ns,nd)
            if nd>1: z=(a.sum(1)/nd)*(b.sum(1)/nd); WN+=W[ic]*a*b.sum(1,keepdims=True)/nd; WS+=W[ic]*b*a.sum(1,keepdims=True)/nd
            else: z=a[:,0]*b[:,0]; WN+=W[ic]*(a*b); WS=WN
            Z+=W[ic]*z.sum()/ns
        lnZ=np.log(Z) if Z>0 else -np.inf
        out[h]=dict(lnZ=float(lnZ),perm=perm,ns=ns,nd=nd,WN=WN,WS=WS,si=si,di=di)
    return out
def cellvals(arr,res,which):
    # map (shared,diff) grid onto flat model array
    perm=res['perm']; a=arr.transpose(perm).reshape(res['ns'],res['nd']); return a
def wq(v,w,qs=(0.16,0.5,0.84)):
    o=np.argsort(v); v=v[o]; c=np.cumsum(w[o]); c/=c[-1]; return [float(np.interp(q,c,v)) for q in qs]
def summarise(name,S,share13):
    res=run(name,S,share13); best=max(r['lnZ'] for r in res.values()); rows=[]; js={}
    for h,r in res.items():
        d=dict(hyp=h,lnZ=r['lnZ'],lnBF_vs_best=r['lnZ']-best)
        for tag,Wt in (('N',r['WN']),('S',r['WS'])):
            if r['nd']==1 and h!='all shared': pass
        # posterior predictive (marginal over everything) of observables, per nucleus
        if np.isfinite(r['lnZ']) and r['lnZ']>-1e300:
            for tag,Wt in (('N',r['WN']),('S',r['WS'])):
                w=Wt.ravel()
                if w.sum()<=0: continue
                for k in ('R10','R32','R43','H13_10','O13_10','H13_43','fH10','fO10','fH43','Tk','logn','tau32','tau10H'):
                    d[f'{k}_{tag}']=wq(cellvals(M[k],r,tag).ravel(),w)
                lx=LX[np.indices(shape)[3]]; d[f'logX_{tag}']=wq(cellvals(lx,r,tag).ravel(),w)
                d[f'logNH_{tag}']=wq(cellvals(lnH,r,tag).ravel(),w); d[f'logX13_{tag}']=wq(cellvals(l13,r,tag).ravel(),w); d[f'P(X13=100)_{tag}']=float(w[(cellvals(l13,r,tag).ravel()==2.0)].sum()/w.sum())
            # sign of T_N - T_S and offset in X when that axis is in D (use pair weights is expensive: use marginals only -> independent approx)
        if r['nd']>1 or h=='all four differ':
            for ax,lab in ((0,'T'),(1,'n'),(3,'X'),(2,'NH')):
                vals={0:M['Tk'],1:M['logn'],3:LX[np.indices(shape)[3]],2:lnH}[ax]
                a=cellvals(vals,r,'N')
                pN=r['WN'].ravel(); pS=r['WS'].ravel(); vN=a.ravel()
                if pN.sum()>0:
                    # independence approximation between N and S marginals (shared calibration ignored)
                    uN,iN_=np.unique(vN,return_inverse=True); wN=np.bincount(iN_,pN)/pN.sum(); wS=np.bincount(iN_,pS)/pS.sum()
                    d[f'P({lab}_N>{lab}_S)']=float(sum(wN[i]*wS[:i].sum() for i in range(len(uN)))); d[f'P({lab}_N=={lab}_S)']=float((wN*wS).sum())
        rows.append(d)
    return rows
if __name__=='__main__':
    ALL={}
    for nm,S in SETS.items():
        for sh in (True,False):
            if not S['c13'] and not sh: continue
            rows=summarise(nm,S,sh); ALL[f'{nm}|X13 {"shared" if sh else "free"}']=rows
            print(f'\n== {nm}  X13 {"shared" if sh else "free per nucleus"}')
            for d in rows: print(f"  {d['hyp']:38s} lnZ {d['lnZ']:8.2f}  lnBF vs best {d['lnBF_vs_best']:7.2f}")
    json.dump(ALL,open(O+f'joint_summary{TAG}.json','w'),indent=1)
