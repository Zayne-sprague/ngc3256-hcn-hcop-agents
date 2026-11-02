# Q4c-8: predicted HCN/H13CN 4-3 Tb ratio for the posterior draws of the joint fits (classes A, B, F2; 'all' data; physT30), N(H13CN)=N(HCN)/X with 12C/13C-type abundance ratio X=30,50,90 (assumed), same Tkin, n(HCN), filling;
# compared with the measured 3-sigma lower limits (ours >8.3 N, >8.3 S [wing-model conservative]; team-b p0007 /snapshots/team-b-s04 s41/s45/s49: >8.3 N (8-13), >4 S (4-7)). Output out/q4c/h13cn_pred.csv
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4c_joint as J
G=np.load(J.O+'grid4_h13cn.npz'); TRX=G['TR'][...,0]
rows=[]
for c,lim in (('N',8.3),('S',4.0)):
    for cl in ('A both thick, kappa=1','B both thick, same Tex (0.8-1.25), kappa free','F2 separate n, kappa=1, both thick'):
        P=J.posterior(c,'all',[cl]); w=P['w']/P['w'].sum()
        it=np.searchsorted(J.TK,P['T']); inn=np.round((P['logn_H']-J.LN[0])/0.25).astype(int); iN=np.round((P['logNH']-J.LNC[0])/0.25).astype(int)
        trH=J.TH[it,inn,iN,3]
        for X in (30,50,90):
            iNx=np.round((P['logNH']-np.log10(X)-J.LNC[0])/0.25).astype(int)    # nearest grid column (0.25 dex step) for N(HCN)/X
            ok=iNx>=0
            tr13=np.where(ok,TRX[it,inn,np.clip(iNx,0,None)],np.nan)
            # interpolate linearly in log N for the fractional offset
            f=((P['logNH']-np.log10(X)-J.LNC[0])/0.25); i0=np.floor(f).astype(int); fr=f-i0; i0c=np.clip(i0,0,len(J.LNC)-2); 
            tr13=TRX[it,inn,i0c]*(1-fr)+TRX[it,inn,i0c+1]*fr; tr13=np.where(i0>=0,tr13,np.nan)
            r=trH/tr13; m=np.isfinite(r)&(r>0); ww=w[m]/w[m].sum(); rr=r[m]
            o=np.argsort(rr); cw=np.cumsum(ww[o]); q=lambda p:np.interp(p,cw,rr[o])
            rows.append(dict(region=c,cls=cl,X=X,ratio_med=q(.5),ratio_lo95=q(.025),ratio_hi95=q(.975),limit_to_beat=lim,frac_consistent_with_limit=float(ww[rr>=lim].sum()),frac_valid=float(m.mean()),alpha_flux_med=0.949*1/q(.5)))
R=pd.DataFrame(rows); R.round(3).to_csv(J.O+'h13cn_pred.csv',index=False); pd.set_option('display.width',200); print(R.round(2).to_string())
