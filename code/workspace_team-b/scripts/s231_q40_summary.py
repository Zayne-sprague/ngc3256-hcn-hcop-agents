# Q40 summary of s230 outputs: chi2min, cold/hot-branch dchi2, posterior Tk/n, A=N(HCN)/N(HCO+) per position and N/S contrast (marginal + profile). usage: python3 -I s231_q40_summary.py -> out/q40/summary.txt/json
import numpy as np,json,os
O='/workspace/out/q40/'; ST=0.05; lag=np.arange(-20,29)
TK=np.logspace(1,np.log10(300),30); LN=np.round(np.arange(3,8.001,0.1),2)
def quant(h,qs=(0.16,0.5,0.84),lg=lag):
    x=np.concatenate([[lg[0]-0.5],lg+0.5])*ST; c=np.concatenate([[0],np.cumsum(h)/h.sum()]); return [10**float(np.interp(q,c,x)) for q in qs]
def contrast(hA,hB,pA,pB):
    d={};p={}
    for i in range(len(lag)):
        for j in range(len(lag)):
            k=lag[i]-lag[j]; d[k]=d.get(k,0)+hA[i]*hB[j]; p[k]=min(p.get(k,1e9),pA[i]+pB[j])
    k=np.array(sorted(d)); return k,np.array([d[x] for x in k]),np.array([p[x] for x in k])
out={};L=[]
def P(s): print(s); L.append(s)
for v,Xs in (('N',[60]),('H',[25,30,40,60,90]),('A',[25,30,40,60,90]),('B',[40,60,90])):
  for X in Xs:
    R={}
    for pos in 'NSC':
        f=f'{O}fit_{pos}_{v}_X{X}.npz'
        if not os.path.exists(f): continue
        R[pos]=np.load(f)
    if len(R)<3: continue
    P(f'--- variant {v}, 12C/13C={X}')
    for pos,r in R.items():
        cm=r['chimin']; W=r['Wc']; W=W/W.sum()
        pcold=W[TK<=30].sum(); pT=W.sum(1); pn=W.sum(0)
        qT=[float(TK[np.searchsorted(np.cumsum(pT),q)]) for q in (.16,.5,.84)]; qn=[float(LN[np.searchsorted(np.cumsum(pn),q)]) for q in (.16,.5,.84)]
        A=quant(r['h_A']); pa=r['p_A']; pb=10**(lag[pa.argmin()]*ST)
        row=dict(lnZ=float(-0.5*r['chibest']+np.log(r['Wc'].sum())),chi2=float(r['chibest']),dchi2_cold=float(cm[TK<=30].min()-cm.min()),dchi2_mid=float(cm[(TK>30)&(TK<100)].min()-cm.min()),P_cold=float(pcold),Tk=qT,logn=qn,A=A,A_prof=float(pb))
        out[f'{v}|X{X}|{pos}']=row
        P(f'  {pos}: chi2min {row["chi2"]:.2f}, lnZ {row["lnZ"]:.2f}; dchi2 best(Tk<=30 K)-best = {row["dchi2_cold"]:.2f}, (30<Tk<100) {row["dchi2_mid"]:.2f}; posterior P(Tk<=30)={pcold:.3f}; Tk {qT[1]:.0f} [{qT[0]:.0f},{qT[2]:.0f}], log n {qn[1]:.1f} [{qn[0]:.1f},{qn[2]:.1f}]; N(HCN)/N(HCO+) {A[1]:.2f} [{A[0]:.2f},{A[2]:.2f}], profile best {pb:.2f}')
    for a,b in (('N','S'),('C','S')):
        k,h,p=contrast(R[a]['h_A'],R[b]['h_A'],R[a]['p_A'],R[b]['p_A'])
        q=quant(h,lg=k); pg=(h[k>0].sum()+0.5*h[k==0].sum())/h.sum(); d0=p[k==0][0]-p.min()
        out[f'{v}|X{X}|{a}/{b}']=dict(med=q[1],lo=q[0],hi=q[2],Pgt1=float(pg),prof_best=float(10**(k[p.argmin()]*ST)),dchi2_at1=float(d0))
        P(f'  {a}/{b} N(HCN)/N(HCO+): {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}] P(>1)={pg:.3f}; profile best {10**(k[p.argmin()]*ST):.2f}, dchi2(=1)={d0:.2f} -> {np.sqrt(d0):.1f} sigma')
open(O+'summary.txt','w').write('\n'.join(L)); json.dump(out,open(O+'summary.json','w'),indent=1)
