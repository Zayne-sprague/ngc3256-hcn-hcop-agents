# Q53 round 2 work summary
Script: scripts/s332_q53_commongrid.py (re-runnable, `python3 -I`); output out/q53/s332.txt/.json. Old: s330, s331 (13CO).
Method: s320 recipe (pbcor, RJ Tb, circular 2.5"/1.4" kernel, fixed N/S pixels, linear baseline 600<|v|<1200, other lines masked) but integrals now come from cumulative-flux interpolation on a common velocity axis (fractional channels).
Results (details in report.json, 6 Q53 entries):
- core50 inversion (N/S 0.82) was an artefact of per-line channel selection; common grid gives HCN/HCO+ N/S 1.36 (|v|<50), 1.48 (<100), 1.55 (<150) at 1-0; 3-2: 2.0-2.1 (1.4"), 1.8-1.85 (2.5"). HNC 1.6-1.8, CN 2.0 (1-0).
- Per-bin profiles (>5 sigma both lines only): contrast in every bin, no sign change. Sensitivity: v_sys +-20 and baseline order/range <1% on core; bin edges dominate (HCN 1.36-1.83).
- Wings 100-350 (1-0): HCN/HCO+ N 1.23, S 0.545 (N/S 2.25+-0.22, syst. 2.0-2.4). Null test: wing integrals have 1-3x formal noise at N. N-rich wing not reproduced at 3-2 (N wing 0.63 = core); CN wing ratios baseline-dependent (N/S 2.8-5.7), dropped.
- 3-Gaussian joint fit: broad HCO+ is 22% of S flux (11 K km/s inside |v|<100) vs 4% at N. Removing it: N/S 1.30 (HCN-free outflow) to 1.65 (rho=1); observed S-wing HCN/HCO+ 0.43-0.54 gives 1.43. HCO+ boost <=~12% of the contrast; equalising would need 28 K km/s HCN-free HCO+ in the S core.
- Opacity/excitation: matched 2.5" R32/10 core HCN N 0.306, S 0.217; HCO+ 0.404/0.359. CN 1-0 G1/G2 N 0.40, S 0.49 (thin 0.50): CN thin at both.
- Wording softened: "N/S line-ratio deficit", abundance not separated from excitation.
Literature: Michiyama+2018 ApJ 868, 95 (abstract: HCO+/CO 0.20+-0.04 in S outflow vs 0.08+-0.01 in S nucleus). Did not verify their HCN/HCO+ per-velocity values (the 0.43 S-wing figure is ours/the reviewer's).
Uncertain: decomposition model dependence; CN G1 blends; 13CO(2-1) opacity; no matched H13CN/H13CO+ in this item; R32/10 differs from the earlier "equal at both nuclei" statement (not reconciled).
