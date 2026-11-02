# Q4b-3: HCN/HCO+ opacity at S (and N) as a model comparison of four physical interpretations, with a joint RADEX fit (common Tkin, n; N_HCN, N_HCO+ free) to the full ratio set at matched components:
#   data: HCN 3-2/2-1, HCN 4-3/3-2, HCO+ 3-2/2-1, HCO+ 4-3/3-2 (Q1b deblended component ladders vs Xb4f 3-2, 14% calibration + deblend spread; out/q1/q1_deblend_ladders_xb4f.csv)
#         HCN/HCO+ 3-2 peak ratio (S 0.43+-0.06 [0.46 at 0.25", 0.40 at 1.05", 0.36 1" templates]; N 0.85+-0.10) ; filling prior eta<=1 from the beam-dilution intrinsic Tint (out/q4b/dilution.csv, 20% error).
#   classes (tau = RADEX line-centre tau of the 3-2 line at dv=150 km/s):
#     A both thick (tau>1), equal filling (kappa=eta_HCN/eta_HCO+=1): R32 comes from different Tex (sub-thermal HCN)
#     B both thick, Tex(HCN)/Tex(HCO+) within 0.8-1.25 (thermalised/same Tex), kappa free (<=1, different filling)
#     C HCN thin (tau32<1), HCO+ thick, kappa=1
#     D both thin, kappa=1 (ratio = abundance ratio)
#   Flat priors in log T, log n, log N_HCN, log N_HCO+ (and log kappa in B). Reports min chi2, posterior mass (Bayes-factor-like, prior dependent), allowed ranges and predictions of each class.
# Priors: "physT30" = phys with Tkin>=30 K (starburst/AGN nuclei; gas not colder than dust); "grid" (flat on the RADEX grid) and "phys" (log N <= log n + 12.5, i.e. path length <= 100 pc for X<=1e-8). Output out/q4b/scenarios_summary.csv. Also a 'no 2-1' variant (drop the 3-2/2-1 ratios, which depend on the Band-5 deblending).
import sys,warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
O='/workspace/out/q4b/'
H=np.load(O+'grid_hcn.npz'); Q=np.load(O+'grid_hco+.npz'); TK,LN,LNC=H['TK'],H['LN'],H['LNC']
L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv'); DIL=pd.read_csv(O+'dilution.csv')
def lad(c,sp,q): r=L[(L.comp==c)&(L.species==sp)&(L.quantity==q)].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
R32={'S':(0.43,0.06),'N':(0.85,0.10),'Bridge':(0.52,0.07)}
def data(c,drop21=False):
    d=[]
    for sp in ('HCN','HCO+'):
        for q in ('3-2 [Xb4f]/2-1 (Tb)','4-3/3-2 [Xb4f 3-2] (Tb)'):
            v,s=lad(c,sp,q); d.append((v,np.hypot(0.14*v,s)))
    d.append(R32[c]); d=np.array(d)       # order: H32/21, H43/32, O32/21, O43/32, R32
    return d
TRH,TRO=H['TR'],Q['TR']; TAH,TAO=H['tau'],Q['tau']; TXH,TXO=H['Tex'],Q['Tex']
# arrays over (T,n,NH,NO)
def b(a,i,mode): return a[...,i][:,:,:,None] if mode=='H' else a[...,i][:,:,None,:]
hA=b(TRH,1,'H')/b(TRH,0,'H'); hB=b(TRH,2,'H')/b(TRH,1,'H'); oA=b(TRO,1,'O')/b(TRO,0,'O'); oB=b(TRO,2,'O')/b(TRO,1,'O')
r32m=b(TRH,1,'H')/b(TRO,1,'O'); r21m=b(TRH,0,'H')/b(TRO,0,'O'); r43m=b(TRH,2,'H')/b(TRO,2,'O')
shape=r32m.shape; sh=lambda a: np.broadcast_to(a,shape)
tH=sh(b(TAH,1,'H')); tO=sh(b(TAO,1,'O')); tH43=sh(b(TAH,2,'H')); xH=sh(b(TXH,1,'H')); xO=sh(b(TXO,1,'O'))
ok=np.isfinite(sh(hA))&np.isfinite(sh(hB))&np.isfinite(sh(oA))&np.isfinite(sh(oB))&np.isfinite(sh(r32m))&np.isfinite(tH)&np.isfinite(tO)&np.isfinite(tH43)&np.isfinite(xH)&np.isfinite(xO)&np.isfinite(sh(b(TRH,1,'H')))&np.isfinite(sh(b(TRO,1,'O')))&(sh(b(TRH,0,'H'))>0)&(sh(b(TRO,0,'O'))>0)&(sh(b(TRO,1,'O'))>0)
TT=np.broadcast_to(TK[:,None,None,None],shape); NN=np.broadcast_to(LN[None,:,None,None],shape); NH=np.broadcast_to(LNC[None,None,:,None],shape); NOo=np.broadcast_to(LNC[None,None,None,:],shape)
LPHYS=12.5   # physical prior: N/n <= X*L with X<=1e-8 and L=100 pc (3.1e20 cm): log N <= log n + 12.5
KAP=np.round(10**np.linspace(-1.5,0,13),4)
def wq(x,w,q): 
    i=np.argsort(x); c=np.cumsum(w[i]); c/=c[-1]; return np.interp(q,c,x[i])
def analyse(c,drop21,prior):
    pm=np.ones(shape,bool) if prior=='grid' else ((NH<=NN+LPHYS)&(NOo<=NN+LPHYS)&((TT>=30) if prior=='physT30' else True))
    d=data(c); use=np.array([not drop21,True,not drop21,True,True])
    dl=DIL[DIL.region==c].set_index('line'); TiH=dl.loc['HCN 3-2','Tint_K']; TiO=dl.loc['HCO+ 3-2','Tint_K']
    TRh=sh(b(TRH,1,'H')); TRo=sh(b(TRO,1,'O'))
    rows=[];preds=[]
    classes={'A both thick, kappa=1 (HCN sub-thermal/different Tex)':(pm&ok&(tH>1)&(tO>1),[1.]),
             'B both thick, same Tex, kappa free (different filling)':(pm&ok&(tH>1)&(tO>1)&(xH/xO>0.8)&(xH/xO<1.25),KAP),
             'C HCN thin, HCO+ thick, kappa=1':(pm&ok&(tH<1)&(tO>1),[1.]),
             'D both thin, kappa=1 (abundance ratio)':(pm&ok&(tH<1)&(tO<1),[1.]),
             'E thermalised n>=1e7 (any tau), kappa free (not pooled)':(pm&ok&(NN>=7),KAP)}
    allw={}
    for nm,(m,kl) in classes.items():
        best=(1e30,None); W=[]; Dat=[]
        for k in kl:
            P=[sh(hA),sh(hB),sh(oA),sh(oB),k*sh(r32m)]
            chi=sum(((P[i]-d[i,0])/d[i,1])**2 for i in range(5) if use[i])
            etaO=TiO/TRo; etaH=TiH/TRh; pen=(np.clip(etaO-1,0,None)/0.2)**2+(np.clip(etaH-1,0,None)/0.2)**2
            # for kappa=1 classes eta_H and eta_O must agree: that is the R32 data term; kappa<1 (B) lets them differ
            chi=np.where(m,chi+pen,1e9)
            i=np.unravel_index(np.argmin(chi),chi.shape)
            if chi[i]<best[0]: best=(chi[i],(i,k,[p[i] for p in P]))
            W.append((np.exp(-0.5*(chi-chi.min()*0))*m,k,P,etaH,etaO))
        # posterior mass relative to global min later
        allw[nm]=(best,W,kl)
    gmin=min(v[0][0] for v in allw.values())
    def summarise(nm,W,best,kl,scale_kl=True):
        mass=sum(w.sum() for w,*_ in W)/len(kl)*np.exp(0.5*gmin)   # relative to the best class point, divided by kappa prior count
        wt=[];vals={k:[] for k in ('T','logn','logNH','logX','tau32H','tau43H','tau32O','TexH','TexO','etaH','etaO','r21','r43','R43pred','H32/21','H43/32','O32/21','O43/32','R32','kappa','logNbeam025H')}
        for w,k,P,eH,eO,sc in W:
            s=w>1e-6*w.max() if w.max()>0 else w>0; wt.append(w[s]*sc)
            for key,arr in (('T',TT),('logn',NN),('logNH',NH),('logX',NH-NOo),('tau32H',tH),('tau43H',tH43),('tau32O',tO),('TexH',xH),('TexO',xO),('etaH',eH),('etaO',eO),('r21',k*r21m),('r43',k*r43m),
                            ('H32/21',P[0]),('H43/32',P[1]),('O32/21',P[2]),('O43/32',P[3]),('R32',P[4]),('logNbeam025H',NH+np.log10(np.clip(eH,1e-3,None)))):
                vals[key].append(arr[s])
            vals['kappa'].append(np.full(s.sum(),k)); vals['R43pred'].append((k*r43m)[s])
        wt=np.concatenate(wt)
        row=dict(prior=prior,region=c,drop21=drop21,scenario=nm,chi2_min=best[0],ndata=int(use.sum()),npar=4+(nm[0] in 'BE'),mass_rel=mass)
        if wt.sum()>0:
            for key in vals:
                v=np.concatenate(vals[key]); wm=wt[np.isfinite(v)]; v=v[np.isfinite(v)]
                lo,md,hi=[wq(v,wm,q) for q in (0.025,0.5,0.975)]
                row[key+'_med']=md; row[key+'_lo95']=lo; row[key+'_hi95']=hi
        i,k,pb=best[1]; row.update(dict(best_T=TT[i],best_logn=NN[i],best_logNH=NH[i],best_logNO=NOo[i],best_tau32H=tH[i],best_tau32O=tO[i],best_tau43H=tH43[i],best_TexH=xH[i],best_TexO=xO[i],best_kappa=k,
            best_etaO=TiO/TRo[i],best_etaH=TiH/TRh[i],pred_H32_21=pb[0],pred_H43_32=pb[1],pred_O32_21=pb[2],pred_O43_32=pb[3],pred_R32=pb[4],pred_R43=k*r43m[i],pred_R21=k*r21m[i]))
        return row
    pool=[]; pb=None
    for nm,(best,W,kl) in allw.items():
        W2=[(w,k,P,eH,eO,1./len(kl)) for (w,k,P,eH,eO) in W]
        rows.append(summarise(nm,W2,best,kl))
        if nm[0]!='E': pool+=W2
        if nm[0]!='E' and (pb is None or best[0]<pb[0]): pb=best
    rows.append(summarise('ALL pooled (A-D, prior-weighted)',pool,pb,[1.]))
    return rows
allrows=[]
for c in ('S','N','Bridge'):
    for d21 in ((True,) if c=='Bridge' else (False,True)):
        for prior in ('grid','phys','physT30'): allrows+=analyse(c,d21,prior)
R=pd.DataFrame(allrows)
for (c,d,pr),g in R.groupby(['region','drop21','prior']): R.loc[g.index,'mass_norm']=g.mass_rel/g.mass_rel[~g.scenario.str.startswith(('ALL','E'))].sum()
R.loc[R.scenario.str.startswith(('ALL','E')),'mass_norm']=np.nan
R.round(4).to_csv(O+'scenarios_summary.csv',index=False)
pd.set_option('display.width',300,'display.max_columns',80,'display.max_colwidth',40)
cols=['prior','region','drop21','scenario','chi2_min','ndata','npar','mass_norm','best_T','best_logn','best_logNH','best_tau32H','best_tau32O','best_tau43H','best_TexH','best_TexO','best_kappa','best_etaH','best_etaO','pred_H32_21','pred_H43_32','pred_O32_21','pred_O43_32','pred_R32','pred_R43','pred_R21']
print(R[cols].round(3).to_string())
print(R[['prior','region','drop21','scenario']+[c for c in R.columns if c.endswith(('_med','_lo95','_hi95')) and c.split('_')[0] in ('T','logn','logNH','logX','tau32H','tau43H','tau32O','etaH','etaO','R43pred','r21','logNbeam025H')]].round(2).to_string())
