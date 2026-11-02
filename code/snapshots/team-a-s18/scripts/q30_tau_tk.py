# Q30: tau(12CO) from T12/T13 and conditional Tk table. Assumes same Tex, ff, beam for both isotopologues; R=12C/13C
import numpy as np, json
from scipy.optimize import brentq
H=6.62607e-27;KB=1.380649e-16;NU=230.538e9;TBG=2.73
J=lambda T: H*NU/KB/(np.exp(H*NU/KB/T)-1); Jinv=lambda j: H*NU/KB/np.log1p(H*NU/KB/j)
d=json.load(open('out/q30/co13_ratio.json')); tk=json.load(open('out/q23b/tk_ff.json'))
f=lambda t,R:(1-np.exp(-t))/(1-np.exp(-t/R))
def tau(r,R):
    if r>=R: return 0.0
    return brentq(lambda t:f(t,R)-r,1e-3,500)
out={}
for s in 'SN':
    # ratio inputs at 1.0": peak-channel and HCN-FWHM mean; conservative +-: 10% flux-cal/missing-flux differential (13CO more resolved-out) -> +-15% in ratio, plus noise
    for q in ['ratio_peak','ratio_FWHMmean']:
        r=d[f'1.0_{s}_core'][q]; e13=d[f'1.0_{s}_core']['e13_FWHMmean']/d[f'1.0_{s}_core']['T13_FWHMmean']
        row={'ratio':r}
        for R in [40,60,90]:
            row[f'R{R}']={'tau':[tau(r*u,R) for u in (0.85,1,1.15)]}   # 15% systematic range (low ratio -> high tau)
            # reverse: low ratio => larger tau
            row[f'R{R}']['tau']=[tau(r*u,R) for u in (0.85,1,1.15)][::-1][::-1]
        # Tk_min at 0.26" for peak/HCNmean with measured tau (ff=1,0.5,0.3), Q23 T_RJ(0.26)
        T=tk[s]['T_RJ_0.26']['peak' if q=='ratio_peak' else 'HCNmean']
        for R in [40,60,90]:
            tt=row[f'R{R}']['tau']
            row[f'R{R}']['Tk_min']={f'ff{ff}':[float(Jinv(T/(ff*(1-np.exp(-t)))+J(TBG))) for t in tt] for ff in [1.0,0.5,0.3]}
        out[f'{s}_{q}']=row
        print(s,q,round(r,1)); 
        for R in [40,60,90]: print('  R',R,'tau(ratio*0.85,1,1.15)',np.round(row[f'R{R}']['tau'],2),'Tk_min ff1',np.round(row[f'R{R}']['Tk_min']['ff1.0'],1),'ff.5',np.round(row[f'R{R}']['Tk_min']['ff0.5'],1),'ff.3',np.round(row[f'R{R}']['Tk_min']['ff0.3'],1))
json.dump(out,open('out/q30/tau_tk.json','w'),indent=1)
