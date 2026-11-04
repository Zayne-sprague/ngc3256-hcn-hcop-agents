# Q59 work summary (worker, team-b): H13CN/HCN and H13CO+/HCO+ 1-0 at N and S, common 2.2-2.5" beam

## What I did
- Checked prior work first: team-a p0059 (cites 2015.1.00993 13C data at 2.2", says H13CN/HCN at ~2" is decisive, measured N/S "3+-3.5"); our own Q29 s185 (template/box fits at 2.2"), Q42b, Q44c. Reproduced the Q29 numbers, then built a Q56-style systematic budget instead of duplicating.
- s350: 384 variants = beam {2.2, 2.5"} x aperture {pixel, 0.5 theta, theta} x aperture centre {HCO+ or HCN +-150 peak, S-centred} x window {+-100/150/250, HCO+ template} x baseline {lin, const}; weak line from a joint fit with all lines in the window (SiO, HN13C, CCH with H13CO+; H13CN alone), errors from line-free rms x null-aperture calibration (floor = null std), +5% cross-EB scale.
- s351: inference: double ratio D, 13C-based N(HCN)/N(HCO+) contrast, 12C opacity vs 12C/13C prior; (a) continuum cross-EB check is impossible (cubes continuum-subtracted); X95d/X961 are consecutive EBs of one night, X961 matches X145 and Harada Table 13 to <=5% (Q29 s183).
- s352: report.json updated (3 new Q59 entries; Q29/Q42b opacity entries annotated). 250 results; backup report_backup_q59.json.

## Scripts / outputs
scripts/s350_q59_isotopologue_budget.py, s351_q59_inference.py, s352_update_report_q59.py; out/q59/s350.json, s350.txt, s351.json, s351.txt. Rerun: `python3 -I scripts/s350...; python3 -I scripts/s351...` (needs /archive 2015.1.00993 cubes X961 spw0, X95d spw0/spw1).

## Results (mean over variants; error = sqrt(median stat^2 + variant spread^2))
| ratio | N | S |
|---|---|---|
| H13CN/HCN | 0.0221+-0.0075 (stat 0.0044) | 0.0014+-0.0098 (<0.031, 3 sigma; not detected) |
| H13CO+/HCO+ | 0.0284+-0.0058 | 0.0128+-0.0064 |
Peak pixel (2.2"): H13CN/HCN N 0.031+-0.006, S -0.001+-0.013; H13CO+/HCO+ N 0.038+-0.005, S 0.0135+-0.0074. Harada+18 Table 13 (N, 1.7" peak): H13CO+/HCO+ 0.046+-0.005: consistent with a compact N core.
- D = N/S of the 13/12 ratio: H13CN/HCN 5.8 [2.1-13.5], 95% lower bound 1.2, P(D<1.2)=0.05; H13CO+/HCO+ 3.1 [1.7-7.9]. BF (D 2-3 vs 0.6-1.0): 6.7 and 18, but resting on one N detection and 1-1.5 sigma S values: suggestive, not decisive. D>1 for both species points to a nuclear-wide 13C/12C difference (opacity and/or 12C/13C), not an HCN-specific one.
- Thin-13C N(HCN)/N(HCO+) contrast C13 (N/S) = 4.8 [1.2-16.7] vs 12C intensity contrast 1.45; P(C13>1.45)=0.80, P(C13>2.5)=0.67, P(C13<1)=0.12. At N, HCN/HCO+ is 0.83 in 12C and 0.78 in 13C: no evidence HCN is more opaque than HCO+ at N.
- Opacity (equal Tex, thin 13C), X=12C/13C log-flat 40-100: tau(HCO+ 1-0) N median 1.3 [0.45-2.3], S <1.2 (95%); tau(HCN 1-0) N median 0.65 [0-1.7], <2.4 (95%), S <0.35. N observed 12C/13C-equivalents 45 (HCN), 35 (HCO+); S HCO+ 78.
- Verdict on p0059: right test, but 2015.1.00993 cannot decide it (our result matches their "3+-3.5"). Need sigma(H13CN/HCN)_S ~0.003 (~10x integration; X961/X95d depth 3.5-4.3 mK per 27 km/s channel).

## Literature
Harada+2018 (ApJ 855, 49) Table 13 and sect. 5.1 (13C abundance caveat); Tang+2019 (A&A 629, A6) and Martin+ for the 12C/13C range 40-100; Jimenez-Donaire+2017 (MNRAS 466, 49) X=50 convention; team-a p0059, team-b p0044, p0050.

## Uncertain / not done
- S H13CN consistent with zero; D and C13 posteriors are heavy-tailed and prior-dependent in the upper tail.
- Single-zone, equal-Tex, thin-13C assumptions; 13C data at 2.2" vs a ~1" core (Q44c radial hint of a core-concentrated 13C/12C); 12C/13C may differ N vs S (Q40 prefers X_N 31-50, X_S 65-85), degenerate with opacity.
- Cross-EB scale (X95d vs X961) is carried as 5%, not measured. No Q57 grid refit done here: feed these ratios (table above) into Q57/Q44 with the variant spread included.
- X94a (2016.1.00965) not used: it loses 25-45% of the 1-0 flux (resolved out), unsuitable for 13C/12C ratios.
