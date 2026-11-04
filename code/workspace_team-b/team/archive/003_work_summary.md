# Q44 work summary (worker, team-b)

## What was done
Fully matched one-zone RADEX fit (s250b/s235 machinery) with every datum at ONE 2.2" taper, r<2.2" aperture mean, +-150 km/s (SiO 2-1 free), no Harada peak pixel:
12C 1-0 (Q29 s183 X961), 12C 3-2 (Q34 s192 12m+7m), H13CN/H13CO+ 1-0 (Q29 s185_box pm150). CS is only available as Harada 1.7" peak pixel, so it enters scaled (variants MC/M10C, +-25%) and is flagged as non-matched.
The background job (s261_run_all.sh) hit its time limit after M, M12, M10, MC and most of M10C; M32 and M4 were NOT run (M4 is superseded by the Xb47 4-3 fits below; Q42 already showed H13CN 4-3 does not discriminate). Analysis script s262 fixed (X-grid index bug in the X_S/X_N contrast) and rerun.

## Scripts / outputs
scripts/s260_q44_radial.py, s260f_q44_radial_summary.py (radial test); s260b/c/d/e (Xb47-beam 4-3 variants); s261_q44_matched_fit.py + s261_run_all.sh (fits); s262_q44_analysis.py.
out/q44/analysis_q44.txt/json (matched fit), analysis.txt (s260e: 4-3 variants), s260a.txt, s260f.txt (radial), fit_*.npz, s261.log. Report entries: Q44a, Q44b, Q44c in /output/report.json.

## Results
* Q44a (primary M, LN(45,0.2dex) X prior, 10/15 and 20/25% band priors): N P(Tk<=30 K) = 0.013-0.029 (flat log X<=100: 0.020-0.042; flat to X<=178: 0.20), BF hot:cold 13-60 (2 for flat-178); Tk(N) ~130-155 K [60-250]. X_N = 66-77 [51-92]; S X_S = 80-89 [65-97] (lower bound, S H13CN/HCN undetected; ~130 if X<=178). X_S/X_N = 1.07 [0.84,1.43] (<=100), 1.24-1.31 [0.8,2.0] (<=178): no significant contrast. N(HCN)/N(HCO+): N 3.6 [3.1,4.5], S 2.2 [2.0,2.6], Bridge 3.1; N/S 1.6 (~2.3 sigma).
  Versus Q42 (mixed apertures): X_N 33-58 -> 66-77; P(cold, 20/25%) 0.2-0.58 -> 0.03-0.05. The Q42 non-exclusion came from mixing 1.7" peak-pixel 12C with 2.2" aperture 13C/12C. The thick branch is disfavoured (not excluded: flat X to 178 gives 0.2).
* Q44b: P(cold|N) 0.40/0.49 (12C only), 0.33 (1-0 12C+13C only, = prior), 0.33 (+CS 2-1), 0.013/0.029 (adding 3-2), 0.0002-0.011 (+CS 2-1, 5-4 scaled). The decision comes from 3-2 12C together with the 13C 1-0 data, not from 1-0 alone. 12C 4-3 cannot be at 2.2" (beam 4.14x2.70"); as a mismatched-beam supplement it raises dchi2(cold-hot) from 0 to 8-10 at N (P(cold) 0.002-0.033 for 10/15-40/40% band priors at r<2.2"; 0.02-0.36 at the 0.6" peak).
* Q44c radial: N H13CO+/HCO+ 0.0372+-0.0066 (<0.3") -> 0.0339 (<1.1") -> 0.0271 (<2.2") -> 0.0218+-0.0051 (<3.3"); gradient 0.0154+-0.0075 (2.0 sigma). H13CN/HCN gradient 1.5 sigma. S flat (0.015). Harada 0.0458+-0.0049 is consistent with our <0.3" value (1.1 sigma). The factor-2 drop is aperture (a real radial decline of 13C/12C at N, 2 sigma), not calibration.

## Literature
Harada+2018 Table 13 (1.7" peak pixel intensities, H13CO+/HCO+ 0.0458); Tang+2019, Martin+2019 for starburst 12C/13C 21-67 (X_N here, 66-77, is at/above the top of that range, an aperture average); Jimenez-Donaire+2017 (tau13=tau12/X).

## Uncertain
One-zone model with 1-3 dof; X prior and X<=100 truncation drive P(cold) (0.01 to 0.2); CS scaling and Xb47 beam are not matched; S 13C is a non-detection; no core-only matched fit (no 3-2 at <1"); M32/M4 not run.
