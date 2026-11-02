# Q23: pooled direct test (observed red-minus-blue deficit vs ramp-model prediction, in-sample and leave-one-region-out) and HCN 4-3/3-2/HCO+ 2-1 correction summary.
import json, numpy as np
R=json.load(open('/workspace/out/q23/s145.json')); REG=['N','S','Bridge']
print('pooled N,S,Bridge: observed (blue - red) level [mJy/ch] (+ = red side depressed) vs predicted; chi2 of (obs-pred) over the 3 regions')
for k,d in R['direct'].items():
    o=np.array([d['rows'][r]['diff'] for r in REG]); e=np.array([d['rows'][r]['err'] for r in REG])
    # predicted blue-red = -pred_red (blue level modelled as 0)
    p=-np.array([d['rows'][r]['pred_red'] for r in REG]); pl=-np.array([d['rows'][r]['pred_red_loo'] for r in REG])
    m=o.mean(); em=np.sqrt(np.sum(e**2))/3
    print(f'  {k:20s} ({d["best"]:6s}) observed {m:+.2f}+-{em:.2f}; predicted {p.mean():+.2f} (LOO {pl.mean():+.2f}); chi2(obs-pred)={np.sum(((o-p)/e)**2):.1f} / 3 ; chi2(obs-0)={np.sum((o/e)**2):.1f} / 3 ; chi2 LOO={np.sum(((o-pl)/e)**2):.1f}/3')
print('\nHCN 4-3, HCN 3-2, HCO+ 2-1 integral change when the fitted ramp is removed (E0 = s30 estimator), % of raw: shapes tri100 / tri400 / emp')
for key in ('Xb47.spw16.HCN4','Xb4f.spw18.HCN3','Xb55.spw18.HCO+2'):
    for r in REG:
        row=R['est'][key][r]; e0=row['E0']
        s=' | '.join(f'{sh} {100*(row["templ"][sh]["E0"]/e0[0]-1):+.1f}%' for sh in ('emp','tri100','tri200','tri400'))
        print(f'  {key:18s} {r:6s} stat {100*e0[1]/e0[0]:.1f}%  {s}')
