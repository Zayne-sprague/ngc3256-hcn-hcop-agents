# Q47 fix 3: posterior of X_S/X_N (12C/13C of the two nuclei, assumed common to HCN and HCO+) from the 2.2" 1-0 ratios, with the RADEX opacity correction.
# For each nucleus the core (Tk,n,N_HCN,N_HCO+) is drawn from the N-only / S-only nested-sampling posteriors (out/q50/runs, base, 2 seeds pooled, 12C ladders only), the model
# ratio R13/12 = TR(13C 1-0)/TR(12C 1-0) at N(13C)=N/X is computed from RADEX grids (out/q36 12C, out/q47 H13CN) for a grid of X, and L(X)=<exp(-chi2/2)> over the posterior samples (chi2 summed over H13CN/HCN and H13CO+/HCO+;
# error = empirical(max of formal, null scatter) + 10% model + systematic |R_A-R_T|).  Priors on X: logU(20,150) (primary) and LN(45,0.2 dex)<=100 (team-b).  Also the thin limit R=1/X.
import json,glob,numpy as np
from scipy.ndimage import map_coordinates
rs=np.random.default_rng(5); R=json.load(open('/workspace/out/q47/ratios10.json')); G13=np.load('/workspace/out/q47/grid_h13cn_10.npz')
XG=np.exp(np.linspace(np.log(15),np.log(200),120))
def lg(a): a=np.array(a,float); a[~((a>0)&(a<1e3))]=np.nan; return np.log10(a)
def ip(arr,G,lT,n,lN,axN):
    LTK=np.log10(G['TK']); c=np.array([np.interp(lT,LTK,np.arange(len(LTK)))*np.ones_like(lN),np.interp(n,G['LN'],np.arange(len(G['LN'])))*np.ones_like(lN),np.interp(lN,axN,np.arange(len(axN)))*np.ones_like(lT)]); return 10**map_coordinates(arr,c,order=1,mode='nearest')
def samples(reg,N=1500):
    S=[];W=[]
    for f in sorted(glob.glob(f'/workspace/out/q50/runs/q50_{reg}_base_s*.npz')):
        d=np.load(f,allow_pickle=True); nm=list(d['names']); S.append(np.array([d['S'][:,nm.index(k)] for k in(f'lTc_{reg}',f'nc_{reg}',f'NHc_{reg}',f'NOc_{reg}')]).T); W.append(d['w'])
    S=np.vstack(S); W=np.concatenate(W); W/=W.sum(); i=rs.choice(len(W),N,p=W); return S[i]
def lik(reg,est,sysmode='AT'):
    G12={m:np.load(f'/workspace/out/q36/grid_{m}_{reg}_cmb.npz') for m in('hcn','hco+')}; P=samples(reg); lT,n,NH,NO=P.T
    h12=ip(lg(G12['hcn']['TR'][...,0]),G12['hcn'],lT,n,NH,G12['hcn']['LNC']); o12=ip(lg(G12['hco+']['TR'][...,0]),G12['hco+'],lT,n,NO,G12['hco+']['LNC'])
    th=ip(lg(G12['hcn']['tau'][...,0]),G12['hcn'],lT,n,NH,G12['hcn']['LNC']); to=ip(lg(G12['hco+']['tau'][...,0]),G12['hco+'],lT,n,NO,G12['hco+']['LNC'])
    ll=np.zeros((len(XG),len(lT))); thin=[]
    for k,X in enumerate(XG):
        rh=ip(lg(G13['TR']),G13,lT,n,NH-np.log10(X),G13['LNC'])/h12; ro=ip(lg(G12['hco+']['TR'][...,0]),G12['hco+'],lT,n,NO-np.log10(X),G12['hco+']['LNC'])/o12
        c=0
        for m,rm in (('HCN',rh),('HCO+',ro)):
            d=R[f'{m}|{reg}|{est}']; sy=abs(R[f'{m}|{reg}|A']['R']-R[f'{m}|{reg}|T']['R']) if sysmode=='AT' else 0.
            c=c+((rm-d['R'])/np.sqrt(d['err']**2+(0.10*rm)**2+sy**2))**2
        ll[k]=-0.5*np.nan_to_num(c,nan=1e4)
    L=np.exp(ll-ll.max()).mean(1); return L/np.trapezoid(L,XG),dict(tau10_HCN=np.nanpercentile(th,[16,50,84]).round(2).tolist(),tau10_HCOp=np.nanpercentile(to,[16,50,84]).round(2).tolist())
def q(x,w,qs): c=np.cumsum(w); c/=c[-1]; return np.interp(qs,c,x)
out={}
for est in ('A','T','P'):
  for sysm in ('AT','none'):
    LN_,tN=lik('N',est,sysm); LS_,tS=lik('S',est,sysm)
    for pn,pr in (('logU20-150',lambda X:1/X*((X>=20)&(X<=150))),('LN45<=100',lambda X:np.exp(-0.5*(np.log10(X/45)/0.2)**2)/X*(X<=100))):
        pN=LN_*pr(XG); pN/=pN.sum(); pS=LS_*pr(XG); pS/=pS.sum()
        # ratio posterior by Monte Carlo
        xn=rs.choice(XG,200000,p=pN); xs=rs.choice(XG,200000,p=pS); r=xs/xn
        out[f'{est}|{sysm}|{pn}']=dict(XN=q(XG,pN,[.16,.5,.84]).round(1).tolist(),XS=q(XG,pS,[.16,.5,.84]).round(1).tolist(),ratio_q16_50_84=np.percentile(r,[16,50,84]).round(2).tolist(),ratio_q2p5_97p5=np.percentile(r,[2.5,97.5]).round(2).tolist(),
            P_ratio_gt1=float((r>1).mean()),P_ratio_le1p8=float((r<=1.8).mean()),P_ratio_gt1p8=float((r>1.8).mean()),tauN=tN,tauS=tS)
json.dump(out,open('/workspace/out/q47/xratio.json','w'),indent=1)
for k,v in out.items(): print(k,'XN',v['XN'],'XS',v['XS'],'rho',v['ratio_q16_50_84'],v['ratio_q2p5_97p5'],'P>1 %.2f P<=1.8 %.2f'%(v['P_ratio_gt1'],v['P_ratio_le1p8']))
print(out['A|AT|logU20-150']['tauN'],out['A|AT|logU20-150']['tauS'])
