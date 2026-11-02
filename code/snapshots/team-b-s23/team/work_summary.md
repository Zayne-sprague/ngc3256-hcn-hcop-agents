# Q39 work summary (worker, team-b)

## Question
Redo the Q35 RADEX column-ratio fit without the 0.2 dex grid-step jitter; do N(HCN)/N(HCO+) N/S (1.64) and N(HCO+)/N(CS) N/S (0.63) survive at >2 sigma; add a two-gas (CS from a lower-density gas) variant.

## What I did
- `scripts/s220a_q39_fine_grid.py`: new RADEX LVG grid (HCN, HCO+, CS; LAMDA, dv 150 km/s, Tbg 2.73 K) Tk 10-300 K (30 log steps), log n 3-8 (0.1), log N 11.5-19 (0.1) -> `out/q39/grid_fine.npz`. Q35's grid dropped/NaN'd nothing but was 0.2-0.25 dex; TR<=0 points (Tex<Tbg) are floored (a first version dropped whole cells and gave wrong fits; fixed).
- `scripts/s220b_q39_fine_fit.py` (+ `s220c_q39_summary.py`): data = Harada+2018 Table 13 1.7" (six lines, as s210c). Columns of HCN, HCO+, CS are independent and summed exactly on a 0.05 dex grid (log-flat priors), Tk/n on the fine grid, ln eta and band scales (10%/15%) marginalised analytically. Gives marginal posteriors AND profile likelihoods of the column ratios (no jitter). Positions independent. `NCAP`/`TAG` env vars give N-cap sensitivity runs (17, 18): same results.
- `scripts/s220d_q39_twogas.py` (+ `s220e_q39_tg_summary.py`): explicit eta grid and 7-node band-scale quadrature so that N and S share the band scale, and so that one-gas (g3), gas 1 = HCN+HCO+ (g1) and gas 2 = CS (g2, own Tk, n, eta) are comparable by chi2, BIC and evidence; n2<=n1 variant.
- `scripts/s221_update_report_q39.py` wrote three Q39 entries to /output/report.json and a supersession note to the Q35 column-ratio entry (backup `report_backup_q39.json`). Report validates (0 errors).
- Outputs: `out/q39/summary*.txt` (independent positions, band scale marginalised "sb" and fixed "nosb"), `out/q39/tg_summary.txt/json` (shared band scale; one vs two gas), logs.

## Results
- One-gas best chi2 on the fine grid: N 5.1, S 0.57, C 0.94 (Q35: 10.2, 2.1, 7.9). The "chi2 10 at N" was mostly grid coarseness. Best fit at N sits at the Tk=300 K cap, log n ~4.5; N posterior is bimodal (cold/dense optically thick branch at Delta chi2 ~4.3, where N(HCN)/N(HCO+) is unconstrained 0.1-2.5).
- N(HCN)/N(HCO+) N/S: Q35 1.64 [1.31,2.09] P=1.00. Fine grid one gas, shared band scale: marginal 1.5 [0.42,2.6], P(>1)=0.73; profile best 1.78, 2.1 sigma. Band scale fixed: 1.54 [1.27,1.83], 2.4 sigma marginal / 2.1 profile. HCN+HCO+ only: 2.3 [1.5,3.1], 1.4-1.7 sigma profile; two-gas 2.4 [1.6,3.3], 1.6 sigma. So the value is stable (1.5-2.4) but significance is ~2 sigma, not P=1.00.
- N(HCO+)/N(CS) N/S: Q35 0.63. Fine grid: 0.75 [0.64,0.87] with band scales fixed (marginal 2.7 sigma-equiv.), but profile 1.6 sigma, and 0.80 [0.36,6.0] with a shared band-scale marginalisation. Does NOT survive as >2 sigma. N(HCN)/N(CS) N/S 1.4 (0.9 sigma): no departure. Bridge/S contrasts are all <1 sigma (~1.0-1.2).
- One gas vs two gas: two-gas chi2 0.03/0.04/0.08 (9 par) vs 5.1/0.57/0.94 (6 par); BIC (n=6) N 16.2 vs 15.8, S 16.2 vs 11.3, C 16.2 vs 11.7; ln B(two/one) +1.7 N, -1.3 S, -0.3 C (n2<=n1: +1.0, -2.2, -2.2). No significant preference; the CS column is unconstrained in the two-gas model, so no CS-referenced ratio exists there. Tk is unconstrained at all three positions.

## Literature
Harada+2018 Table 13 (parsed in lit/harada_table13_parsed.csv) as data input only; no published column ratios exist for NGC 3256 N/S.

## Uncertain / remaining
- With 6 data points and about as many parameters the likelihood is non-Gaussian, so sigmas are indications; the profile and marginal answers differ (N branch).
- Tk <= 300 K cap and log-flat priors in n, N; 4.3" ratios (Q35 +4.3 variant) not re-included; Harada 1.7" values not on a common beam with the CS 7-6 deblend.
- What would settle it: H13CN/H13CO+ (or deeper 3 mm isotopologues) to break the optically thick branch; CS 2-1/5-4 on a common beam; a lower-J HCN/HCO+ (1-0) tied to a Tk tracer (e.g. CH3CN/H2CO) to constrain temperature.
