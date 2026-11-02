# Q4b-1: beam-dilution of the peak Tb of HCN/HCO+ 3-2 across the three beams (Xb4b 0.25", Xb4d 0.95", Xb4f 5.7x3.5") at N, S, Bridge.
# Model: round Gaussian source of FWHM s with intrinsic (cloud-ensemble averaged) brightness Tint; peak Tb in an elliptical beam (bmaj,bmin) = Tint * s^2/sqrt((s^2+bmaj^2)(s^2+bmin^2)).
# Two beams (0.25, 0.95) fix (s, Tint) exactly; Xb4f (5.7x3.5") is then a prediction (excess = extended/neighbour/low-level gas, not part of the compact source).
# Errors: Monte Carlo of the peak Tb (channel rms + 10% absolute flux scale per band-OUS, independent), 4000 draws; draws where no solution exists (ratio <= beam-area ratio) counted.
# Output out/q4b/dilution.csv
import numpy as np, pandas as pd
T=pd.read_csv('/workspace/out/q4/peak_tb.csv'); T=T[T.line.isin(['HCN 3-2','HCO+ 3-2'])]
rng=np.random.default_rng(1)
def F(s,b): return s**2/np.sqrt((s**2+b[0]**2)*(s**2+b[1]**2))
def solve(t1,t2,b1,b2):
    from scipy.optimize import brentq
    r=t1/t2; g=lambda s: F(s,b1)/F(s,b2)-r
    lo,hi=1e-3,5.
    if g(lo)*g(hi)>0: return np.nan,np.nan
    s=brentq(g,lo,hi); return s,t1/F(s,b1)
rows=[]
for reg in ('N','S','Bridge'):
    for line in ('HCN 3-2','HCO+ 3-2'):
        d={o:T[(T.region==reg)&(T.line==line)&(T.ous==o)].iloc[0] for o in ('Xb4b','Xb4d','Xb4f')}
        B={o:tuple(float(x) for x in d[o].beam.split('x')) for o in d}
        sig={o:np.hypot(d[o].eTb_K,0.10*d[o].Tb_peak_K) for o in d}
        s0,T0=solve(d['Xb4b'].Tb_peak_K,d['Xb4d'].Tb_peak_K,B['Xb4b'],B['Xb4d'])
        ss=[];tt=[];bad=0
        for i in range(4000):
            a=[rng.normal(d[o].Tb_peak_K,sig[o]) for o in ('Xb4b','Xb4d')]
            s,t=solve(a[0],a[1],B['Xb4b'],B['Xb4d'])
            if np.isfinite(s): ss.append(s);tt.append(t)
            else: bad+=1
        ss=np.array(ss);tt=np.array(tt); p=lambda x:np.percentile(x,[16,50,84])
        pred=T0*F(s0,B['Xb4f']) if np.isfinite(s0) else np.nan
        rows.append(dict(region=reg,line=line,Tb_025=d['Xb4b'].Tb_peak_K,Tb_095=d['Xb4d'].Tb_peak_K,Tb_57=d['Xb4f'].Tb_peak_K,s_fwhm_arcsec=s0,s_16=p(ss)[0],s_84=p(ss)[2],s_pc=s0*194, # 40 Mpc: 194 pc/arcsec
            Tint_K=T0,Tint_16=p(tt)[0],Tint_84=p(tt)[2],Tb57_pred_compact=pred,Tb57_excess_over_compact=d['Xb4f'].Tb_peak_K-pred,frac_no_solution=bad/4000,
            dilution_025=F(s0,B['Xb4b']),dilution_095=F(s0,B['Xb4d'])))
R=pd.DataFrame(rows); R.round(3).to_csv('/workspace/out/q4b/dilution.csv',index=False); pd.set_option('display.width',250,'display.max_columns',30); print(R.round(3).to_string())
