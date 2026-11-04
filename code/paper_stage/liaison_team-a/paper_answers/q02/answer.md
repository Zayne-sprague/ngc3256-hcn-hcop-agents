# q02 – isotopologue budget, double ratio, 12C/13C and opacity at N

**Short answer.** Neither of the Appendix A.7 numbers is in our team's work, and I cannot reproduce them. I could not see Appendix A.7 or the draft, so I compared against our own outputs. Our measured double ratios are smaller and much less significant, and they do not exclude an abundance-only difference (double ratio = 1). Section 3.3 ("consistent with abundance-only and with different 12C/13C") is the statement our work supports. Please do not use "5.8 [2.1-13.5], 95% lower bound 1.2" for HCN, or "3.1 [1.7-7.9]" for HCO+, unless the team that wrote them can point to a script that produces them.

## What our work measured (J=1-0, 2.2" aperture; `scripts/q47_ratios10.py` -> `out/q47/ratios10.json`; Q47(b) in the report)
Three estimators, T = template fit, A = aperture, P = peak pixel. Values are 13C/12C ± 1σ, and the double ratio is (13/12)_N/(13/12)_S:

| line | est. | N | S | double ratio, nominal | MC (Gaussian errors) median [16-84%] | P(≤1) |
|---|---|---|---|---|---|---|
| H13CN/HCN | T | 0.0217±0.0079 | 0.0070±0.0068 | 3.1 | 2.6 [1.3-7.0] | 0.08 |
| H13CN/HCN | A | 0.0204±0.0058 | 0.0060±0.0091 | 3.4 | 2.3 [1.1-6.6] | 0.09 |
| H13CN/HCN | P | 0.0319±0.0044 | 0.0216±0.0061 | 1.5 | 1.5 [1.1-2.1] | 0.08 |
| H13CO+/HCO+ | T | 0.0192±0.0039 | 0.0114±0.0092 | 1.7 | 1.5 [0.9-3.7] | 0.22 |
| H13CO+/HCO+ | A | 0.0276±0.0049 | 0.0144±0.0054 | 1.9 | 1.9 [1.3-3.1] | 0.04 |
| H13CO+/HCO+ | P | 0.0306±0.0031 | 0.0327±0.0043 | 0.9 | 0.9 [0.8-1.1] | 0.65 |

(Script: `paper_answers/q02/q02_calc.py`. The MC ratio is ill-defined when S is consistent with 0. S ≤ 0 occurs in 11-26% of draws for the T/A estimators of HCN, and in 11% for HCO+ T.)

- The double ratio depends on the estimator. At S the peak-pixel values (0.022, 0.033) disagree with the T/A values, and we treated that as a systematic (Q47(b)).
- Combined with a RADEX opacity treatment in Q47(b), the posterior for X_S/X_N = (12C/13C)_S/(12C/13C)_N has median 0.98, 16-84% 0.65-1.58 and 95% interval 0.45-2.55. The same bound is 0.72-1.39 for the median and 1.9-3.4 for the 97.5% upper limit, depending on estimator and prior. So there is **no evidence that the 12C/13C ratios differ, and an abundance-only difference (ratio 1) is not excluded**. The earlier team-b bound X_S/X_N ≤ 1.8 held with only 70-96% probability (`scripts/q47_xratio.py` -> `out/q47/xratio.json`).
- H13CN 4-3 gives only limits, and they depend on the CO 3-2 wing model: H13CN/HCN < 0.09-0.11 at N and < 0.135-0.15 at S (3σ), so HCN/H13CN > 9 at N and > 6.7 at S (`scripts/h13cn.py` -> `out/h13cn_limits.csv`). They don't constrain a double ratio.

Because "which budget is correct" cannot be decided by A.7, I recommend citing the table above (ratios ≈1.5-3.4 for HCN and ≈1-2 for HCO+, all with 1σ errors comparable to the signal), and stating that the lower 95% bound is below 1 for every estimator.

## 12C/13C at N if the 13C lines are optically thin (assumption: R = 13C/12C line ratio; ratios given in the request)
- H13CN/HCN = 0.022 gives 1/R = **45**. H13CO+/HCO+ = 0.028 gives 1/R = **36**. For the HCN number the uncertainty is 0.0217±0.0079 (T) or 0.0204±0.0058 (A), i.e. about 30-35%. For HCO+, 0.028±0.005 (A) or 0.019±0.004 (T), about 18-20%. So 12C/13C ≈ 45 (+27/-14) and 36 (+8/-6) at 1σ.
- This applies only if the 12C line is also thin. If the 12C lines are thick, these numbers are **lower limits**: R = τ13/(1-e^-τ12) ≥ 1/X, so X = (τ12/(1-e^-τ12))/R ≥ 1/R.

## Implied opacity at N
Solving τ12/(1-e^-τ12) = X·R with τ13 = τ12/X (thin 13C), for an assumed intrinsic 12C/13C = X:

| X | HCN (R=0.022): τ12 | τ13 | HCO+ (R=0.028): τ12 | τ13 |
|---|---|---|---|---|
| 40 | 0 | 0 | 0.23 | 0.006 |
| 60 | 0.6 | 0.010 | 1.2 | 0.019 |
| 66 | 0.8 | 0.012 | 1.4 | 0.021 |
| 90 | 1.6 | 0.017 | 2.3 | 0.025 |
| 135 | 2.8 | 0.021 | 3.7 | 0.027 |

- The ratios alone cannot fix τ: they give τ12 ≈ 0-3 over X = 40-135. The opacity is degenerate with the unknown X, so the result is τ12(HCN) ≈ 1-3 and τ12(HCO+) ≈ 1-4 for X ≈ 60-135. Our own X_N posterior is 66-135 (16-84%), giving τ12 ≈ 0.8-2.8 (HCN) and 1.4-3.7 (HCO+). With X = 1/R both lines are thin.
- This is lower than our RADEX posterior (Q47(b)): τ10(HCN) = 4.6 [3.0-8.7] and τ10(HCO+) = 0.5 [0.2-4.7] at N, and τ10(HCN) ≈ 5-12 at both nuclei in the RADEX grid. The simple thin-13C inversion therefore tends to give lower opacity than the excitation modelling for HCN, and it is consistent for HCO+. The two treatments differ, and we would quote the RADEX range as the model-based one and the inversion as a consistency check.
- The 12C/13C needed for τ(HCN) ≈ 5 is larger than 135 (extrapolating the table), and 13C lines would then not be thin either (τ13 > 0.04), so this cannot be pinned down with the 1-0 ratios alone.

## Assumptions
2.2" aperture, J=1-0; 13C lines thin, and for the opacity table a single-zone, same-excitation 12C and 13C line (R = τ-ratio only); X is not measured here; ratio uncertainties are the 1σ in `ratios10.json` (formal/empirical, larger used).
