# Q39 summary: posterior + profile-likelihood statistics of column ratios and N/S (and C/S) contrasts from s220b outputs (out/q39/fit_*.npz).
# usage: python3 -I s220c_q39_summary.py  [suffix list]   -> prints + out/q39/summary.json
import numpy as np,json,sys
from scipy.stats import norm
O='/workspace/out/q39/'; LAGS={'A':np.arange(-20,29),'CO':np.arange(-40,41),'CH':np.arange(-40,41)}; ST=0.05
def load(pos,v,sb): return np.load(O+f'fit_{pos}_{v}_{sb}.npz')
def quant(lag,h,qs=(0.16,0.5,0.84)):
    x=np.concatenate([[lag[0]-0.5],lag+0.5])*ST; c=np.concatenate([[0],np.cumsum(h)/h.sum()]); return [float(np.interp(q,c,x)) for q in qs]
def sub(hA,hB,lag):   # distribution of a-b for independent a,b on the same lag grid
    n=len(lag); out=np.zeros(2*n-1); 
    for i in range(n): out[i:i+n]+=hA[i]*hB[::-1]
    return np.arange(-(n-1),n)+2*lag[0]+ (n-1) - (n-1) , out   # lag index: a_i - b_j = lag_i - lag_j
def contrast(hA,hB,lag):
    n=len(lag); d={}
    for i in range(n):
        for j in range(n): d[lag[i]-lag[j]]=d.get(lag[i]-lag[j],0)+hA[i]*hB[j]
    k=np.array(sorted(d)); h=np.array([d[x] for x in k]); return k,h
def profc(pA,pB,lag):
    n=len(lag); d={}
    for i in range(n):
        for j in range(n): d[lag[i]-lag[j]]=min(d.get(lag[i]-lag[j],1e99),pA[i]+pB[j])
    k=np.array(sorted(d)); return k,np.array([d[x] for x in k])
res={}
def stat_row(name,sign,key,pos,v,sb):
    f=load(pos,v,sb); lag=LAGS[key]; h=f['h_'+key]; p=f['p_'+key]
    if sign<0: lag=-lag[::-1]; h=h[::-1]; p=p[::-1]
    q=quant(lag,h); pmin=p.min(); k=p<pmin+1
    return lag,h,p,q,(lag[k].min()*ST,lag[k].max()*ST),float(lag[p.argmin()]*ST)
def report(sb):
    for v,stats in (('one',[('N(HCN)/N(HCO+)',1,'A'),('N(CS)/N(HCO+)',1,'CO'),('N(CS)/N(HCN)',1,'CH'),('N(HCO+)/N(CS)',-1,'CO'),('N(HCN)/N(CS)',-1,'CH')]),('hh',[('N(HCN)/N(HCO+)',1,'A')])):
      for nm,sg,key in stats:
        rows={}
        for pos in 'NSC': rows[pos]=stat_row(nm,sg,key,pos,v,sb)
        for pos in 'NSC':
            lag,h,p,q,pr,pb=rows[pos]; print(f'[{v},{sb}] {pos} {nm}: median {10**q[1]:.2f} [{10**q[0]:.2f},{10**q[2]:.2f}] profile best {10**pb:.2f}, dchi2<1 range [{10**pr[0]:.2f},{10**pr[1]:.2f}]')
            res[f'{v}|{sb}|{pos}|{nm}']=dict(med=10**q[1],lo=10**q[0],hi=10**q[2],prof_best=10**pb,prof1=[10**pr[0],10**pr[1]])
        for a,b in (('N','S'),('C','S')):
            la,ha,pa,_,_,_=rows[a]; lb,hb,pb_,_,_,_=rows[b]
            assert (la==lb).all()
            k,hh=contrast(ha,hb,la); pgt=(hh[k>0].sum()+0.5*hh[k==0].sum())/hh.sum(); q=quant(k,hh)
            kp,pp=profc(pa,pb_,la); dchi=pp[kp==0][0]-pp.min(); 
            # one-sided profile significance toward ratio>1 (sign of best contrast)
            sgn=np.sign(kp[pp.argmin()]) 
            print(f'   {a}/{b} contrast {nm}: median {10**q[1]:.2f} [{10**q[0]:.2f},{10**q[2]:.2f}]  P(>1)={pgt:.3f} (gaussian-equivalent {norm.ppf(max(pgt,1-pgt)):.1f} sigma);  profile: best {10**(kp[pp.argmin()]*ST):.2f}, dchi2(ratio=1)={dchi:.2f} -> {np.sqrt(dchi):.1f} sigma')
            res[f'{v}|{sb}|{a}/{b}|{nm}']=dict(med=10**q[1],lo=10**q[0],hi=10**q[2],Pgt1=float(pgt),prof_best=10**(kp[pp.argmin()]*ST),dchi2_at1=float(dchi))
import os
TAG=os.environ.get('TAG','')
for sb in ('sb','nosb'): report(sb+TAG)
json.dump(res,open(O+f'summary{TAG}.json','w'),indent=1)
