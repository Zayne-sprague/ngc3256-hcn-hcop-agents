"""Q38: append the disk-decline robustness result to /output/report.json and flag the stale wording in entries 96, 111, 133, 99. usage: python3 -I s241_update_report_q38.py"""
import json
R=json.load(open('/output/report.json'))
res=R['results']
def find(prefix):
    for x in res:
        if x['quantity'].startswith(prefix): return x
    raise KeyError(prefix)
NOTE=(' [Q38 REVISION: the "disk" 0.52/0.34/0.14 (2-1/3-2/4-3) are E1/E2/W clump POSITION values at a 10" beam (7-8" from the N-S midpoint), NOT 8-12" annulus ratios (those are NO-GO for 3-2 and 4-3, see the disk-annuli entry). '
      'Robust: the sign and size of the 2-1 -> 3-2 step at E1 (R32/R21 = 0.64 at the position, 0.65-0.69+-0.15 from de-contaminated clump totals); insensitive to PB model (<2% on the ratio) and to zero level for Xb55/Xb4f (<3%). '
      'NOT robust: the 4-3 value 0.14 (E1 0.155, W 0.05; PB>=0.7 for a 12-m Gaussian PB excludes E1/E2/W; zero-level offsets of +-1x the far-field mean move E1 to 0.07-0.22; the de-contaminated E1 HCN 4-3 total is 0.00+-0.16 Jy km/s, R43<0.06): treat as an upper limit R43<~0.2. '
      'Not distinct from the nuclei: in uv-complete matched 2.2" data (12m+7m 3-2, 12m 1-0) R32/R10 = 0.74+-0.07 (E1), 0.89+-0.06 (W), 0.63+-0.11 (E2) vs 0.75+-0.01 at N, so "disk declines faster than the nuclei" is NOT supported at 2.2"; the stronger 10" decline carries nuclear-wing leakage. See Q38 entry.]')
for p in ('Matched-J HCN/HCO+ at a common circular 10"','Observed HCN and HCO+ J=2-1, 3-2, 4-3 ladders at the common 10"','Disk (E1) vs centre excitation contrast'):
    x=find(p)
    if '[Q38 REVISION' not in x['value']: x['value']+=NOTE
x=find('Disk (E1) vs centre excitation contrast'); x['confidence']='low'
x=find('Disk (E1) vs centre excitation contrast')
res.append(dict(
 quantity='Q38: is the disk (E1/E2/W clump) HCN/HCO+ decline with J robust? PB mask, zero level, nuclear leakage, short spacings; 10" ladder, de-contaminated clump totals and uv-complete matched 2.2" 1-0/3-2',
 target='clumps E1 (9" E of N), E2, W; N for reference; team-a p0032/p0033/p0043 challenge',
 value=('PARTLY REVISED. (1) Agreement with team-a p0032: 8-12" annulus ratios at 3-2 and 4-3 are not measurable (already NO-GO in our annulus entry); the 0.52/0.34/0.14 are clump-position values at 10", not an annulus ladder. '
  '(2) 10" ratios (K km/s) at E1/E2/W: 2-1 0.518/0.536/0.508; 3-2 Xb4f 0.329/0.339/0.338 (Xb4d 0.319/0.277/0.255); 4-3 0.155/0.191/0.048. Zero-level offsets (+-1x far-field mean): 2-1 <=0.013, Xb4f <=0.004; Xb4d 0.22-0.39 (E1; its far-field level is at PB~0.2, unreliable); Xb47 E1 0.07-0.22. PB-model effect on the ratio <=1.6%. '
  '(3) De-contaminated clump totals (NNLS, free constant): R21 E1 0.35+-0.04, W 0.49+-0.05, E2 0.44+-0.14; R32 E1 0.24+-0.05 (Xb4f) / 0.23+-0.04 (Xb4d), W 0.33+-0.05 (Xb4f), E2 <0.12 (Xb4f); R43 E1 <0.06, W <0.33, E2 <0.17 (2 sigma): R32/R21 = 0.69+-0.15 (E1 Xb4f), 0.65+-0.14 (Xb4d), 0.68+-0.13 (W). '
  '(4) uv-complete matched 2.2" (X961 12m 1-0; Q34 12m+7m 3-2): R10 E1 0.49+-0.02, E2 0.46, W 0.49, N 0.84; R32 E1 0.365+-0.03, E2 0.29+-0.05, W 0.43+-0.03, N 0.63; R32/R10 E1 0.74+-0.07, E2 0.63+-0.11, W 0.89+-0.06, N 0.75+-0.01; r31 (HCN, HCO+) clump/N: E1 0.80+-0.07, 0.81+-0.03, E2 0.34, 0.40, W 0.55, 0.46 (12m-only 3-2 recovers much less at E2: HCO+ 1.9 vs 4.7 K km/s, so E2/W r31 are lower limits on recovery). '
  'CONCLUSION: a modest decline from 2-1 to 3-2 at the E1 and W clumps (R32/R21 ~0.65-0.7+-0.15) survives; it is consistent with, not in tension with, p0043 (extended core-removed R_ext ratio 1.0 +0.3/-0.4) at ~1-2 sigma, and with the 1-0->3-2 step of the nucleus (0.75). The 4-3 clump ratio is a limit (<~0.2), and "disk declines faster than the nuclei" is not established; downgrade to: HCN/HCO+ at the clumps falls with J by ~25-35% from 1-0/2-1 to 3-2, similar to N, with 4-3 unconstrained.'),
 uncertainty='10" positions: stat 0.003-0.02, zero level +-0.01 (Xb55/Xb4f), +-0.1 (Xb4d, Xb47), PB <2%, nuclear leakage biases the 10" 1-0 and 2-1 ratios high by ~0.1-0.15 (R10 0.63 at 10" vs 0.49 at 2.2" for E1); clump totals: sign-flip + variant spread; 2.2": stat only plus 12m-only 1-0 recovery unquantified; calibration 5-10% per band is common to same-band ratios but enters R32/R21 and R32/R10 (add ~10%)',
 unit='ratio',
 assumptions=[dict(name='positions',value='E1/E2/W nominal positions of regions.py; 2.2" apertures re-centred on the HCO+ peak within 1"'),
  dict(name='zero_level',value='additive offset = f x far-field mean of each line map (f=-1..1), per-pixel line-free baseline otherwise; clump totals with free-sign constant'),
  dict(name='PB',value='12-m Gaussian PB model and a 7-m-like (FWHM x1.7) model bracket the true PB; same-J ratios change <2%'),
  dict(name='conversion',value='RJ T_b; 2-1/3-2/4-3 Jy ratio x (nu_HCO+/nu_HCN)^2; 2.2" ratios computed in K km/s from pbcor cubes')],
 confidence='low',
 evidence='Q38: scripts/s240_q38_disk_clumps.py -> out/q38/s240.txt/json; scripts/s241_update_report_q38.py; inputs s30 maps, q12b/clumpfit, q34 12m+7m cubes, 2015.1.00993.S X961'))
json.dump(R,open('/output/report.json','w'),indent=1)
