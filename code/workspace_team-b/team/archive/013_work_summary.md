# Q40 work summary: do 13C isotopologue data exclude the cold/dense thick branch of the N RADEX posterior?

## What I did
- Extended the Q39 fine-grid one-gas RADEX fit (s220a grid, s220b machinery, imported and monkeypatched) with 13C lines: same Tk, n, filling, band scale; column N/X, X=12C/13C = 25,30,40,60,90 (assumption; NGC 3256 value unknown).
- Data: Harada+2018 Table 13 (lit/harada_table13_parsed.csv, 1.7" peak pixel) H13CO+ 1-0 and 3-2 at N, H13CO+ 1-0 at C (S not tabulated); H13CN 1-0 and S/C values from our Q29 short-spacing-complete 2.2" maps (out/q29/s185.json), applied to the Harada 12C intensities. H13CN 3-2 is not covered by Xb4f (259.0 GHz is in the gap); H13CN 4-3 (Xb47) is blended with the CO 3-2 wing and its limits (Q4/Q16) are far weaker than the 1-0 detections, so it was not used. Isotopologue data enter in linear K km/s (non-detections OK).
- Variants: N (none, = Q39), H (Harada H13CO+ only, matched beam, cleanest), A (H + Q29 H13CN, S, C), B (Q29 aperture ratios for all).

## Scripts / outputs
scripts/s230_q40_isotopologue_fit.py (fits), s231_q40_summary.py (summary), s232_q40_residuals.py (per-datum model vs data at branch best fits), s233_q40_lowX.py (X=25,30), s234_update_report_q40.py (report). Outputs: out/q40/summary.txt/json, residuals.txt, fit.log, lowX.log, fit_*.npz.

## Results
- Without 13C: cold/dense branch dchi2=4.5 above hot (P(Tk<=30 K)=0.72).
- Harada H13CO+ pair at N: R(1-0)=21.8, R(3-2)=35, 3-2/1-0 = 0.28+-0.06. dchi2(cold - hot) = 32 (X=25), 27 (30), 17 (40), 5.2 (60), -1.3 (90). Cold branch predicts H13CO+ 3-2/1-0 = 0.4-0.7, hot 0.25-0.35. Adding H13CN (A): 41/38/28/10/0.2.
- So: thick branch excluded for X<=60 (>=2.3 sigma, >=4 sigma for X<=40). At X=90 not excluded but then N is poorly fit (chi2 17-23 for 8-9 data vs 4.6 without) so X>~70 in a one-zone model is itself disfavoured. N prefers X~30-40.
- Hot branch: Tk 170-300 K (cap), log n 4.6-4.8; N(HCN)/N(HCO+) N 3.5 [3.1,3.8] (X=40), 3.8 [3.3,4.7] (X=60); S 2.2 [1.9,2.6]; C 2.3 [2.1,2.6].
- N/S: 1.6-1.8, marginal P(>1) rises 0.52 -> 0.99, but profile significance stays 2.1 sigma (limited by 12C line S/N). Variant A gives 2.8-4.4 sigma but S fits badly (chi2 12-45 at X<=60; S isotopologue data want X>~90), so not adopted. Bridge/S ~1.0-1.1.

## Literature
Harada+2018 Table 13 (no published NGC 3256 N/S RADEX columns or H13CN data in it). Own Q29 (p0007-era H13CN/H13CO+ at 2.2"), Q39 grid.

## Uncertain
12C/13C unknown (results depend strongly on it); single-zone model; beam mismatch 1.7" vs 2.2"; S has no matched-beam H13CO+ (would need a 1.7" Band 3 map; Harada cubes not available); Tk capped at 300 K; no 13C LAMDA files (HCO+/HCN files reused, ~10%). What would settle it: matched-beam H13CN/H13CO+ 1-0 and 3-2 at S (deeper), or 12C/13C from 13CO/C18O.
