# Q39: summarise s220d outputs: chi2 / evidence / BIC for one-gas vs two-gas, and column-ratio posteriors + profiles with a band scale shared between positions.
# usage: python3 -I s220e_q39_tg_summary.py > out/q39/tg_summary.txt
import numpy as np,json
from scipy.stats import norm
O='/workspace/out/q39/'; ST=0.05; NN=151; nT,nn,nE=30,51,73
from numpy.polynomial.hermite_e import hermegauss
x,wq=hermegauss(7); wq=wq/wq.sum(); SB=(0.10,0.15)
Z3,Z6=np.meshgrid(x*SB[0],x*SB[1],indexing='ij'); ZW=np.outer(wq,wq).ravel(); ZPEN=((Z3/SB[0])**2+(Z6/SB[1])**2).ravel()
LAGS={'A':np.arange(-20,29),'CO':np.arange(-40,41),'CH':np.arange(-40,41)}
def ld(pos,m): return np.load(O+f'tg_{pos}_{m}.npz')
def lse(a): m=np.max(a); return m+np.log(np.sum(np.exp(a-m)))
out={}
print('=== chi2 (data + band-scale prior), evidence, BIC (n_data=6) ===')
for pos in 'NSC':
    g3,g1,g2=ld(pos,'g3'),ld(pos,'g1'),ld(pos,'g2')
    # one gas
    chi1=float(g3['CMmin'].min()); lnZ1=np.log((ZW*g3['L'].sum((0,1))).sum())-0.5*(-g3['cb'])*0 - 0.5*float(g3['cb'])-np.log(nT*nn*nE*NN**3)
    # two-gas free (given z)
    L1=g1['L'].sum((0,1)); L2=g2['L'].sum((0,1))
    lnZ2=np.log((ZW*L1*L2).sum())-0.5*(float(g1['cb'])+float(g2['cb']))-np.log((nT*nn*nE)**2*NN**3)
    chi2f=float((g1['CMmin'].reshape(nT*nn,-1).min(0)+g2['CMmin'].reshape(nT*nn,-1).min(0)-ZPEN).min())
    # n2<=n1
    B1=g1['L'].sum(0); B2=g2['L'].sum(0); cB2=np.cumsum(B2,0)    # (n,z); constraint n2<=n1 has prior fraction 0.51
    lnZc=np.log((ZW*(B1*cB2).sum(0)).sum()/0.51)-0.5*(float(g1['cb'])+float(g2['cb']))-np.log((nT*nn*nE)**2*NN**3)
    m1=g1['CMmin'].min(0); m2=g2['CMmin'].min(0)    # (n,z) min over T
    cm2=np.minimum.accumulate(m2,0)
    chic=float((m1+cm2-ZPEN).min())
    bic=lambda c,k: c+k*np.log(6)
    print(f'{pos}: one-gas chi2min {chi1:.2f} (6 par) BIC {bic(chi1,6):.1f} lnZ {lnZ1:.2f} | two-gas free chi2min {chi2f:.2f} (9 par) BIC {bic(chi2f,9):.1f} lnZ {lnZ2:.2f} | two-gas n2<=n1 chi2min {chic:.2f} BIC {bic(chic,9):.1f} lnZ {lnZc:.2f} | lnBF(2g free/1g) {lnZ2-lnZ1:.2f}, lnBF(n2<=n1 /1g) {lnZc-lnZ1:.2f}')
    print('   best fits:',json.loads(str(g3['best'])),'|',json.loads(str(g1['best']))['N'],json.loads(str(g2['best']))['N'])
    # posterior of Tk,n for CS gas (two-gas) and gas1
    for nm,g in (('one-gas',g3),('gas1',g1),('gas2(CS)',g2)):
        w=(ZW*g['L']).sum(2); 
        TK=np.logspace(1,np.log10(300),30); LN=np.round(np.arange(3,8.001,0.1),2)
        pT=w.sum(1)/w.sum(); pn=w.sum(0)/w.sum(); q=lambda p,v:[float(v[np.searchsorted(np.cumsum(p),a)]) for a in (0.16,0.5,0.84)]
        print(f'   {nm}: Tk {q(pT,TK)}  log n {q(pn,LN)}')
# ----- stats with shared z
def corr(Mx,My,Px,Py,L):
    n=Mx.shape[-1]; h=[];p=[]
    for l in L:
        sx,sy=(slice(l,n),slice(0,n-l)) if l>=0 else (slice(0,n+l),slice(-l,n))
        h.append((Mx[:,sx]*My[:,sy]).sum(1)); p.append((Px[:,sx]+Py[:,sy]).min(1)-ZPEN)   # both P contain ZPEN once
    return np.array(h).T,np.array(p).T
def get(pos,model,key,sign):
    g3,g1,g2=ld(pos,'g3'),ld(pos,'g1'),ld(pos,'g2')
    if model=='one': h,p=g3['h_'+key],g3['p_'+key]
    elif model=='g1': h,p=g1['h_A'],g1['p_A']
    elif model=='two':
        if key=='A': h=g1['h_A']*(g2['L'].sum((0,1)))[:,None]; p=g1['p_A']+g2['CMmin'].min((0,1))[:,None]-ZPEN[:,None]
        elif key=='CO': h,p=corr(g2['M_C'],g1['M_O'],g2['P_C'],g1['P_O'],LAGS['CO'])
        elif key=='CH': h,p=corr(g2['M_C'],g1['M_H'],g2['P_C'],g1['P_H'],LAGS['CH'])
    lag=LAGS[key]
    if sign<0: lag=-lag[::-1]; h=h[:,::-1]; p=p[:,::-1]
    return lag,h,p
def quant(lag,h,qs=(0.16,0.5,0.84)):
    x=np.concatenate([[lag[0]-0.5],lag+0.5])*ST; c=np.concatenate([[0],np.cumsum(h)/h.sum()]); return [float(np.interp(q,c,x)) for q in qs]
def contrast(hA,hB,lag):   # shared z: sum_z ZW conv ; hA[z,lag]
    d={}
    for z in range(len(ZW)):
        for i in range(len(lag)):
            if hA[z,i]==0: continue
            for j in range(len(lag)):
                d[lag[i]-lag[j]]=d.get(lag[i]-lag[j],0)+ZW[z]*hA[z,i]*hB[z,j]
    k=np.array(sorted(d)); return k,np.array([d[a] for a in k])
def profc(pA,pB,lag):
    d={}
    for z in range(len(ZW)):
        for i in range(len(lag)):
            if not np.isfinite(pA[z,i]): continue
            for j in range(len(lag)):
                v=pA[z,i]+pB[z,j]-ZPEN[z]
                if v<d.get(lag[i]-lag[j],1e99): d[lag[i]-lag[j]]=v
    k=np.array(sorted(d)); return k,np.array([d[a] for a in k])
def marg(lag,h,p):
    hh=(ZW[:,None]*h).sum(0); pr=p.min(0); q=quant(lag,hh); kk=pr<pr.min()+1
    return [10**v for v in q],10**(lag[pr.argmin()]*ST),(10**(lag[kk].min()*ST),10**(lag[kk].max()*ST))
print('\n=== column ratios (shared band-scale marginalisation; med [16,84]; profile best, dchi2<1 interval) and N/S contrasts ===')
res={}
for model,stats in (('one',[('N(HCN)/N(HCO+)','A',1),('N(CS)/N(HCO+)','CO',1),('N(HCO+)/N(CS)','CO',-1),('N(CS)/N(HCN)','CH',1),('N(HCN)/N(CS)','CH',-1)]),('two',[('N(HCN)/N(HCO+)','A',1),('N(CS)/N(HCO+)','CO',1),('N(HCO+)/N(CS)','CO',-1),('N(CS)/N(HCN)','CH',1),('N(HCN)/N(CS)','CH',-1)])):
    for nm,key,sg in stats:
        D={pos:get(pos,model,key,sg) for pos in 'NSC'}
        for pos in 'NSC':
            lag,h,p=D[pos]; q,pb,pr=marg(lag,h,p); print(f'[{model}] {pos} {nm}: {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}] profile best {pb:.2f} [{pr[0]:.2f},{pr[1]:.2f}]')
            res[f'{model}|{pos}|{nm}']=dict(med=q[1],lo=q[0],hi=q[2],prof_best=pb,prof1=pr)
        for a,b in (('N','S'),('C','S')):
            lag,ha,pa=D[a]; _,hb,pb_=D[b]
            k,hh=contrast(ha,hb,lag); pgt=(hh[k>0].sum()+0.5*hh[k==0].sum())/hh.sum(); q=quant(k,hh)
            kp,pp=profc(pa,pb_,lag); dchi=pp[kp==0][0]-pp.min()
            print(f'   {a}/{b} {nm}: {10**q[1]:.2f} [{10**q[0]:.2f},{10**q[2]:.2f}] P(>1)={pgt:.3f} ({norm.ppf(max(pgt,1-pgt)):.1f} sig-equiv); profile best {10**(kp[pp.argmin()]*ST):.2f}, dchi2(=1) {dchi:.2f} -> {np.sqrt(dchi):.1f} sig')
            res[f'{model}|{a}/{b}|{nm}']=dict(med=10**q[1],lo=10**q[0],hi=10**q[2],Pgt1=float(pgt),prof_best=10**(kp[pp.argmin()]*ST),dchi2_at1=float(dchi))
json.dump(res,open(O+'tg_summary.json','w'),indent=1)
