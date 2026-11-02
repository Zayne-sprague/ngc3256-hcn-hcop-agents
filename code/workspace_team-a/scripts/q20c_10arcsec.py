# Q20c: our own grid-posterior reproduction of team-b p0028 (10" N/S ladders, HCN/HCO+ 2-1,3-2,4-3) with free beam filling.
# data: team-b p0028 10" integrated intensities (K km/s). Grids: out/q4c/grid4_*.npz (RADEX, LVG sphere, dv=150, TR in K; index 0..3 = 1-0..4-3).
# cases: eta=1 ; shared eta (flat ln eta [1e-3,1]) ; separate eta per molecule (|ln q|<=ln10) ; each at Tbg 2.73 and 5 K. Absolute Tb = I/(1.0645 w), w uniform 100-300 km/s.
import numpy as np, itertools, sys
O='/workspace/out/q4c/'
I={'N':{'hcn':[12.38,4.81,1.62],'hco+':[20.27,10.04,4.51]},'S':{'hcn':[11.66,4.26,1.35],'hco+':[20.89,9.95,4.28]}}
ERR=float(sys.argv[1]) if len(sys.argv)>1 else 0.15
def load(suf):
    G={m:np.load(O+f'grid4_{m}{suf}.npz') for m in ('hcn','hco+')}
    return G
LE=np.linspace(np.log(1e-3),0,41); W=np.linspace(100,300,5)
out=[]
for suf in ('','_tbg5'):
    G=load(suf); TK=G['hcn']['TK']; LN=G['hcn']['LN']; LNC=G['hcn']['LNC']
    TRh=G['hcn']['TR'][...,1:]; TRo=G['hco+']['TR'][...,1:]    # (T,n,N,3): 2-1,3-2,4-3
    TAU=G['hcn']['tau'][...,2]
    # joint grid axes: T, n, NH, NO -> (T,n,NH,NO,3) compact: tr_h[T,n,NH], tr_o[T,n,NO]
    ok_h=np.isfinite(TRh).all(-1)&(TRh>0).all(-1); ok_o=np.isfinite(TRo).all(-1)&(TRo>0).all(-1)
    wT=np.gradient(np.log(TK))
    for reg in ('N','S'):
        for case in ('eta1','shared','sep'):
            # evaluate chi2 over (T,n,NH,NO,eta params,w) with log-sum
            ih=np.array(I[reg]['hcn']); io=np.array(I[reg]['hco+'])
            Th=TRh[:,:,:,None,:]; To=TRo[:,:,None,:,:]
            valid=ok_h[:,:,:,None]&ok_o[:,:,None,:]
            # shape ratios (independent of eta)
            def lr(a,b): return np.log(a/b)
            rh=np.stack([lr(Th[...,1],Th[...,0]),lr(Th[...,2],Th[...,1])],-1); ro=np.stack([lr(To[...,1],To[...,0]),lr(To[...,2],To[...,1])],-1)
            dh=np.log([ih[1]/ih[0],ih[2]/ih[1]]); do=np.log([io[1]/io[0],io[2]/io[1]])
            c2=(((rh-dh)/ERR)**2).sum(-1)+(((ro-do)/ERR)**2).sum(-1)
            c2=np.where(valid,c2,np.inf)
            # cross-molecule + absolute terms need eta; margin over eta grids
            lrc=np.log(Th[...,1]/To[...,1])    # model HCN/HCO+ 3-2 TR ratio
            dc=np.log(ih[1]/io[1])
            # eta_h: ln eta from absolute 3-2: ln(Tb)=ln eta+ln TR
            if case=='eta1': etas=[(0.,0.)]
            elif case=='shared': etas=[(e,0.) for e in LE]
            else: etas=[(e,q) for e in LE for q in np.linspace(-np.log(10),np.log(10),9)]
            Lacc=np.full(c2.shape,-np.inf); acc=[]
            lw=np.log(1/ (W.size))
            # store marginal weights for eta
            etaw=np.zeros(len(etas)); Tpost=np.zeros(len(TK)); Pt=np.zeros(c2.shape)
            best_by=np.full(len(TK),np.inf); tau_acc=0; tw=0
            res=[]
            for k,(e,q) in enumerate(etas):
                # eta_hco = e ; eta_hcn = e+q?  q = ln(eta_hcn/eta_hco)
                eh,eo=e+q,e
                c2k=c2+((lrc+(eh-eo)-dc)/ERR)**2
                # absolute: Tb_obs(3-2) = I/(1.0645 w); chi2 on ln
                ab=np.zeros_like(c2k)-np.inf
                # marginalise w
                terms=[]
                for w in W:
                    tbh=ih[1]/(1.0645*w); tbo=io[1]/(1.0645*w)
                    ca=((np.log(Th[...,1])+eh-np.log(tbh))/ERR)**2+((np.log(To[...,1])+eo-np.log(tbo))/ERR)**2
                    terms.append(-0.5*(c2k+ca))
                L=np.logaddexp.reduce(np.stack(terms),0)-np.log(len(W))
                # priors: flat log T (weights wT), flat n and N grid; eta prior flat in ln eta (and flat in q)
                L=L+np.log(wT)[:,None,None,None]
                res.append(L)
            R=np.stack(res)  # (etas,T,n,NH,NO)
            R=np.where(np.isfinite(R),R,-np.inf)
            mx=R.max(); P=np.exp(R-mx); Z=P.sum()
            pT=P.sum((0,2,3,4))/Z; pe=P.sum((1,2,3,4))/Z
            el=np.array([ (e+q)/1 for e,q in etas]); eta_h=np.exp(np.array([e+q for e,q in etas])); eta_o=np.exp(np.array([e for e,q in etas]))
            Pk=P.sum((1,2,3,4))/Z
            # chi2 profile by T (best over all else)
            prof=-2*(np.logaddexp.reduce(R.reshape(R.shape[0],len(TK),-1),axis=(0,2)) if False else 0)
            minc=(-2*(R-np.log(wT)[None,:,None,None,None])).reshape(len(etas),len(TK),-1)
            profT=minc.min((0,2))
            tauw=(P*TAU[None,:,:,:,None]).sum()/Z
            cdf=np.cumsum(pT)
            Tmed=np.interp(0.5,cdf,TK); T16=np.interp(0.16,cdf,TK); T84=np.interp(0.84,cdf,TK)
            etam=np.exp((Pk*np.log(eta_o)).sum()) if case!='eta1' else 1.
            print(f"{reg} tbg{suf or '2.73'} {case:6s} lnZ={mx+np.log(Z):8.2f} T {Tmed:5.0f} [{T16:.0f},{T84:.0f}] P(T>50)={pT[TK>50].sum():.2f} P(T<=30)={pT[TK<=30].sum():.2f} eta_med {etam:.3f} profchi2(T=25)={profT[3]-profT.min():.1f} (T=80)={profT[8]-profT.min():.1f} (min {profT.min():.1f})  <tau32H>={tauw:.1f}")
            # profile at eta in restricted ranges
            if case=='shared':
                for lo,hi in ((1e-3,0.03),(0.05,1),(0.3,1)):
                    sel=(eta_o>=lo)&(eta_o<=hi); print(f"      eta in [{lo},{hi}]: min chi2 by T:", ' '.join(f"{int(t)}K:{minc[sel][:,i].min():.1f}" for i,t in enumerate(TK) if t in (20,25,30,50,80,150)))
