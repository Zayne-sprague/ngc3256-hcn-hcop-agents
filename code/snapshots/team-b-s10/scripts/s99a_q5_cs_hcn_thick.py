# Q5 round 2: N(CS)/N(HCN) implied by CS 7-6/HCN 4-3 (K ratio, Xb47 4.3x2.8") when HCN 4-3 is allowed to be optically thick.
# RADEX LVG grids (out/q3/grid_lvg.npz for HCN; out/q5/radex_cs_so.npz for CS; dv=150 km/s, Tbg 2.73 K, H2 collisions). (Tk, n) restricted to the ladder-allowed region from Q3 (Tk>=56 K, log n 4-6; plus a cool/dense corner 28 K),
# and N(HCN)/dv chosen to give tau(HCN 4-3) = 0.3, 1, 3, 10 (the log N grid is interpolated). Common Tk, n, velocity width and beam filling for CS and HCN assumed.
# Output out/q5/cs_hcn_thick.csv; usage: python3 -I s99a_q5_cs_hcn_thick.py
import numpy as np, pandas as pd
g3=np.load('/workspace/out/q3/grid_lvg.npz'); g5=np.load('/workspace/out/q5/radex_cs_so.npz'); TK,LN,LNC=g3['TK'],g3['LN'],g3['LNC']
R={'N':(0.1630,0.0084*3),'S':(0.1358,0.0077*3),'Bridge':(0.1369,0.0073*3)}   # R_K CS 7-6/HCN 4-3 ; error x3 = empirical-field + aperture/baseline systematic (~10-15%, see s98)
def logN_for(arr,target):       # arr over LNC, monotone rising part
    ok=np.isfinite(arr)&(arr>0); x,y=LNC[ok],np.log10(arr[ok]); i=np.argmax(y); x,y=x[:i+1],y[:i+1]; t=np.log10(target)
    return float(np.interp(t,y,x)) if len(y)>2 and y[0]<=t<=y[-1] else np.nan
rows=[]
for tk in (28,56,80,112,160):
  iT=int(np.argmin(abs(TK-tk)))
  for ln in (4.0,4.5,5.0,5.5,6.0):
    iN=int(np.argmin(abs(LN-ln)))
    tau=g3['hcn'][iT,iN,:,2,2]; TRh=g3['hcn'][iT,iN,:,2,0]; cs=g5['cs'][iT,iN,:,0,:]
    for t_target in (0.3,1.,3.,10.):
        ok=np.isfinite(tau); 
        if ok.sum()<3 or tau[ok].max()<t_target or tau[ok].min()>t_target: continue
        lN=float(np.interp(np.log10(t_target),np.log10(np.clip(tau[ok],1e-30,None)),LNC[ok]))   # tau monotone in N
        Th=float(10**np.interp(lN,LNC[ok],np.log10(np.clip(TRh[ok],1e-30,None))))
        for reg,(r,er) in R.items():
            for lab,rr in (('best',r),('-1s',r-er),('+1s',r+er)):
                lNcs=logN_for(cs[:,0],rr*Th); 
                rows.append(dict(region=reg,Tk=tk,logn=ln,tauHCN43=t_target,logN_HCN=lN,T_R_HCN43=Th,case=lab,logN_CS=lNcs,X_CS_HCN=10**(lNcs-lN) if np.isfinite(lNcs) else np.nan,tauCS=np.nan))
D=pd.DataFrame(rows); D.to_csv('/workspace/out/q5/cs_hcn_thick.csv',index=False)
b=D[(D.case=='best')]
pd.set_option('display.width',200,'display.max_rows',400)
print(b.pivot_table(index=['region','Tk','logn'],columns='tauHCN43',values='X_CS_HCN').round(3).to_string())
print('range of X (best) by tau:',b.groupby('tauHCN43').X_CS_HCN.agg(['min','median','max','count']).round(4))
print('region/tau(Tk>=56, n 1e4-1e5):'); bb=b[(b.Tk>=56)&(b.logn<=5.0)]; print(bb.groupby(['region','tauHCN43']).X_CS_HCN.agg(['min','median','max']).round(3))
