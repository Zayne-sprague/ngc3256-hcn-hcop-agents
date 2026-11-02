# Q40 r2: N(HCN)/N(HCO+) at N marginalised over X (lognormal prior), and N/S contrast, for the three band-scale priors.
# S has no H13CO+ in Harada (variant H) so its fit is X-independent: computed once per prior width at logX=1.6 here.
# usage: python3 -I s237_q40r2_ratio_marg.py -> out/q40r2/ratio_marg.txt
import sys,os,numpy as np
sys.path.insert(0,'/workspace/scripts'); sys.argv=['x']
import s235_q40r2_Xmarg as S
s=S.s; O=S.O; ST=0.05; lag=np.arange(-20,29)
def quant(h,qs=(0.16,0.5,0.84),lg=lag):
    x=np.concatenate([[lg[0]-0.5],lg+0.5])*ST; c=np.concatenate([[0],np.cumsum(h)/h.sum()]); return [10**float(np.interp(q,c,x)) for q in qs]
L=[]
def P(x): print(x); L.append(x)
for k,sg in enumerate(S.SIGS):
    fS=O+f'fit_S_H_s{k}_lX160.npz'
    if not os.path.exists(fS):
        q=s.iso('H')['S']; s.X_iso={'H':q['H'],'O':q['O'],'C':[]}
        np.savez(fS,**s.run_one('S',1.6,'H',sg))
    rS=np.load(fS); hS=rS['h_A']/rS['h_A'].sum()
    for nm,(x0,sx) in {'LN(45,0.2dex)':(np.log10(45),0.2),'LN(60,0.2dex)':(np.log10(60),0.2),'flat':(None,None)}.items():
        lxs=[];lz=[];hs=[]
        for lx in S.LOGX:
            f=O+f'fit_N_H_s{k}_lX{int(round(lx*100))}.npz'
            r=np.load(f); lxs.append(lx); lz.append(-1.5*r['chibest']+np.log(r['Wc'].sum())+(0 if x0 is None else -0.5*((lx-x0)/sx)**2)); hs.append(r['h_A']/r['h_A'].sum())
        w=np.exp(np.array(lz)-max(lz)); w/=w.sum(); hN=(w[:,None]*np.array(hs)).sum(0)
        a=quant(hN); b=quant(hS)
        # N/S contrast distribution
        d={}
        for i in range(len(lag)):
            for j in range(len(lag)): d[lag[i]-lag[j]]=d.get(lag[i]-lag[j],0)+hN[i]*hS[j]
        kk=np.array(sorted(d)); hh=np.array([d[x] for x in kk]); c=quant(hh,lg=kk); pg=hh[kk>0].sum()+0.5*hh[kk==0].sum()
        P(f'band prior {sg}, X prior {nm}: N(HCN)/N(HCO+) at N = {a[1]:.2f} [{a[0]:.2f},{a[2]:.2f}] (marginal over X, Tk, n); S = {b[1]:.2f} [{b[0]:.2f},{b[2]:.2f}]; N/S = {c[1]:.2f} [{c[0]:.2f},{c[2]:.2f}], P(N/S>1)={pg:.2f}')
open(O+'ratio_marg.txt','w').write('\n'.join(L))
