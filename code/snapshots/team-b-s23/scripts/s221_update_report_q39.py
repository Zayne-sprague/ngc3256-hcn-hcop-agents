"""Q39: write Q39 results into /output/report.json (idempotent) and add a supersession note to the Q35 column-ratio entry. usage: python3 -I s221_update_report_q39.py"""
import json
R=json.load(open('/output/report.json')); R['results']=[x for x in R['results'] if not str(x.get('quantity','')).startswith('Q39')]
NOTE=' [Q39 UPDATE: on a fine continuous grid the P(>1)=1.00 / 2.5 sigma statements above are over-precise; N(HCN)/N(HCO+) N/S is ~2.1 sigma (profile likelihood) and N(HCO+)/N(CS) N/S 0.63 is only 1.6 sigma; see the Q39 entries.]'
for x in R['results']:
    if str(x['quantity']).startswith('Q35: column-density ratios') and 'Q39 UPDATE' not in x['value']: x['value']+=NOTE
EV=('/workspace/scripts/s220a_q39_fine_grid.py -> out/q39/grid_fine.npz; s220b_q39_fine_fit.py + s220c_q39_summary.py -> out/q39/summary.txt (+ _cap17.0/_cap18.0 sensitivity); '
    's220d_q39_twogas.py + s220e_q39_tg_summary.py -> out/q39/tg_summary.txt/json (shared band scale, one-gas vs two-gas chi2/BIC/evidence); logs out/q39/*.log')
COM=[dict(name='data',value='Harada+2018 Table 13 peak pixels at 1.7" (lit/harada_table13_parsed.csv): HCN 1-0, 3-2; HCO+ 1-0, 3-2; CS 2-1, 5-4 (K km/s), error = table rms + 5% per line; literature input, cubes not available. The 4.3" ratio constraints of Q35 are not included (they changed nothing in Q35).'),
     dict(name='model',value='RADEX LVG (LAMDA hcn, hco+, cs; H2 collider; dv 150 km/s; Tbg 2.73 K), fresh grid Tk 10-300 K (30 log steps, 0.05 dex) x log n 3-8 (0.1) x log N 11.5-19 (0.1, cubic-interpolated to 0.05 dex); columns of HCN, HCO+, CS free and independent (flat in log N); log n and log Tk flat; ln eta flat on [-3.6,0], one eta per gas; band scales 1-0/2-1 (10%) and 3-2/5-4 (15%) Gaussian. TR<=0 grid points (Tex<Tbg) floored at 1e-3 K. Tk<=300 K cap: the one-gas N best fit sits at the cap.'),
     dict(name='statistics',value='Marginal posterior (log-flat priors, exact sums on the grid, no jitter) and profile likelihood (min chi2 incl. band-scale prior) of log column ratios on a 0.05 dex lag grid. N/S contrasts: band scale shared between N and S (same observations) in s220e; independent (conservative) in s220c. Profile significance = sqrt(chi2(ratio=1)-chi2min), one parameter. With 6 data and ~6 parameters the likelihood is not Gaussian, so these are indications, not p-values.')]
ADD=[
dict(quantity='Q39: N/S contrast of N(HCN)/N(HCO+) from the RADEX fit on a fine continuous (0.05 dex) grid, one-gas (HCN, HCO+, CS common Tk, n, eta) and HCN+HCO+ only / two-gas',
 target='N vs S (and Bridge C vs S); Harada 1.7" lines',
 value=('Q35 value was 1.64 [1.31,2.09], P>1=1.00 on a 0.2 dex grid with +-0.1 dex jitter. Fine grid, one gas, band scale shared N-S (s220e): marginal posterior N/S = 1.5 [0.42,2.6], P(>1)=0.73; profile-likelihood best 1.78, chi2(N/S=1)-chi2min = 4.5 -> 2.1 sigma. '
        'The marginal is broad because the N posterior has a second branch (cold/dense, optically thick lines) in which N(HCN)/N(HCO+) is unconstrained between 0.1 and 2.5 at Delta chi2 ~4.3 above the optimum (A_N = 3.2-3.8, i.e. log 0.5-0.55; S: A=2.0-2.5); that branch is not removed by capping log N at 17 or 18. '
        'Band scales held fixed: 1.54 [1.27,1.83], P=0.993 (2.4 sigma), profile 2.1 sigma (independent-position calc., s220c, N cap 17/18 gives the same). HCN+HCO+ only (4 lines, independent positions): 2.3 [1.5,3.1], profile 1.4 sigma (band scales marginalised) / 2.3 [1.75,3.0], 1.7 sigma (fixed). Two-gas (HCN+HCO+ gas 1, CS gas 2): 2.4 [1.6,3.3], profile 1.6 sigma. '
        'Absolute (one gas, profile best): N(HCN)/N(HCO+) = 3.6 N, 2.0-2.2 S, 2.2 C (A = 2.0-2.5, log 0.3-0.4 at S, C). Bridge/S = 1.02 [0.85,1.23] (one gas), 1.24 [0.92,1.65] (two gas): no contrast. '
        'CONCLUSION: the N/S excess of N(HCN)/N(HCO+) of ~1.5-2.3 is consistent in sign and size with Q35 and with the line-ratio excess (1.5-1.9 at every J, model independent), but its RADEX significance is ~2 sigma (profile), not P=1.00; the 0.2 dex grid had made it look sharper and had biased the one-gas chi2 high.'),
 unit='dimensionless',uncertainty='marginal 16-84% and profile intervals as stated; non-Gaussian, ~2 sigma at best',
 method='fine-grid RADEX marginal posterior and profile likelihood (s220b-e)',evidence=EV,confidence='low',assumptions=COM),
dict(quantity='Q39: N/S contrast of N(HCO+)/N(CS) (and N(HCN)/N(CS), N(CS)/N(HCN)) on the fine RADEX grid: does 0.63 survive?',
 target='N vs S (Bridge C vs S)',
 value=('N(HCO+)/N(CS) N/S, one gas: 0.75 [0.64,0.87] (2.7 sigma-equiv., P(>1)=0.03) if band scales are held fixed (independent-position calc.), profile 1.6 sigma (chi2(=1)-chi2min = 2.6, best 0.79); with shared band-scale marginalisation 0.80 [0.36,6.0], P(>1)=0.38, profile 1.6 sigma. Q35 had 0.63 [0.47,0.84]. '
        'So it does NOT survive as >2 sigma. N(HCN)/N(CS) N/S = 1.44 [0.96,3.6] (profile best 1.41, 0.9 sigma); N(CS)/N(HCN) = 0.69 [0.28,1.04]. Column ratios (one gas, profile best): N(CS)/N(HCN) 0.25 N, 0.35 S, 0.32 C; N(HCN)/N(CS) 4.0 N, 2.8 S, 3.2 C; N(HCO+)/N(CS) 1.1 N, 1.4 S, 1.4 C. '
        'Bridge/S: N(HCO+)/N(CS) 0.88 [0.71,1.05], N(HCN)/N(CS) 0.91 [0.72,1.15] (<1 sigma). In the two-gas model the CS column is unconstrained (two lines, four parameters), so no CS-referenced ratio can be quoted there. '
        'CONCLUSION: with CS as reference HCO+ is lower at N by ~0.6-0.8 only as a 1.6 sigma indication; HCN vs CS does not depart (0.9 sigma). Q35 statement "HCO+-poor, 2.5 sigma" is downgraded to ~1.6 sigma.'),
 unit='dimensionless',uncertainty='see value; profile intervals dchi2<1; non-Gaussian',
 method='s220b-e',evidence=EV,confidence='low',assumptions=COM),
dict(quantity='Q39: one-gas vs two-gas (CS from a separate, lower-density gas) fit at N, S, Bridge: chi2, BIC, evidence; Tk, n',
 target='N, S, Bridge(C); Harada 1.7"',
 value=('One-gas best chi2 (data + band-scale prior, fine grid) N 5.1, S 0.57, C 0.94 (Q35, 0.2 dex grid: 10.2, 2.1, 7.9): the coarse grid inflated the N misfit by ~2x, the one-gas model fits N acceptably (6 lines, ~6 parameters, nominal dof ~0), best fit at Tk = 300 K (the cap), log n ~4.5. '
        'Two-gas free: chi2 0.03 (N), 0.04 (S), 0.08 (C), 9 parameters. BIC (6 data points): N one-gas 15.8 vs two-gas 16.2 (Delta=+0.4); S 11.3 vs 16.2; C 11.7 vs 16.2. Log Bayes factor two-gas(free)/one-gas = +1.7 (N), -1.3 (S), -0.3 (C); with n2<=n1 (CS in the lower-density gas, prior-volume corrected) +1.0 (N), -2.2 (S), -2.2 (C). '
        'No preference for two gases at any position (N weakly favoured, <2.5 in ln B); CS-gas CS column is unconstrained. Tk, n posteriors: one gas N Tk bimodal (14 K [11,237], log n 6.3 [4.6,6.6] cold/dense branch vs the hot branch at the cap), S Tk 117 [29,237] K, log n 4.6 [4.3,5.1]; Bridge Tk 237 [167,300] K, log n 4.3; HCN+HCO+ gas: N Tk 18 [10,104] K, log n 4.0 [3.2,5.0]; S 46 [18,167] K, 4.0 [3.2,4.8]; CS gas N Tk 32 [11,167], log n 4.7 [3.7,5.7]; S 65 [14,211], 4.8 [4.3,5.6]. Temperature is not constrained at N, S or Bridge by these lines.'),
 unit='chi2; K; cm-3',uncertainty='posterior 16-84%; BIC with n=6 is only indicative',
 method='s220d-e',evidence=EV,confidence='low',assumptions=COM)]
R['results']+=ADD
json.dump(R,open('/output/report.json','w'),indent=1)
