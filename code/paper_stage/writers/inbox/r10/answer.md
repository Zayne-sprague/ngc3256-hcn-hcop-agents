# r10 — 4-3 vs 3-2 N/S HCN/HCO+ contrast at matched 1.0" and 1.5"

**Short answer: the data cannot decide it. The rise of the N/S contrast with J is real only if little 4-3 flux is missing. Bound at 1.0": 4-3 contrast 2.3 to 3.7 against 2.1 at 3-2, so the 4-3/3-2 contrast ratio is 1.07 to 1.73. The central pedestal-corrected value is 3.3 vs 2.1, i.e. 1.6.** Fully proportional correction (the maximum the Xb47 comparison allows) removes the rise completely.

Scripts: `r10_ladder.py` (matched-beam extraction) and `r10_corr.py` (missing-flux bounds). Outputs: `ladder.csv` and `ladder_corr.csv`. Both are in this directory.

## Method
- **Estimator.** Same as p0053 (`/workspace/scripts/q53g_fixedwin_all.py`). Window VC±350 km/s (N 2785, S 2750), linear baseline fitted at 450–950 km/s, T_b in K km/s, every cube convolved to a circular beam. Errors are statistical only.
- **4-3.** HCN from 2018.1.00493 (X1df1 spw29), HCO+ from 2018.1.00223 (X858 spw27).
- **3-2.** 2017.1.00379 Xb4b at 1.0" and 1.5", plus Xb4d at 1.5" (Xb4d's smallest taper is 1.1").
- **Missing flux.** The only single-tuning 4-3 pair is Xb47, which has a 4.1×2.7" beam. It cannot be convolved down to 1.0" or 1.5", so it enters only as a flux reference at 4.5" (`/workspace/out/q60/imgplane.csv`).
  - Recovery of the high-resolution 4-3 data relative to Xb47, with a flux of K km/s at 4.5": HCN N 0.876, HCN S 0.607, HCO+ N 0.769, HCO+ S 0.866.
  - The README figure "S HCN is 0.61 of Xb47" reproduces.
  - Two bounding corrections were applied to the 4-3 values:
    - (P) pedestal: add the extended flux (Xb47 − high-res at 4.5") as a uniform T_b offset. The fractional correction shrinks at smaller beams.
    - (F) fractional: divide by the recovery above, as if all lost flux sat inside the core. This is the maximum correction.
  - The 3-2 values are not corrected. Xb4b vs Xb4d differ by up to 14% for S HCN at 1.5" (21.9 vs 25.0 K km/s), which is a 3-2 systematic of that size.

## T_b (K km/s) and 4-3/3-2 at 1.0" (Xb4b 3-2)
| | 3-2 | 4-3 raw | 4-3/3-2 raw | 4-3/3-2 (P) | 4-3/3-2 (F) |
|---|---|---|---|---|---|
| HCN N | 114.3±0.9 | 73.8±1.1 | 0.646±0.011 | 0.655 | 0.737 |
| HCN S | 36.6±1.0 | 17.2±1.1 | 0.471±0.034 | 0.515 | 0.775 |
| HCO+ N | 136.8±1.2 | 90.7±1.6 | 0.663±0.013 | 0.691 | 0.862 |
| HCO+ S | 92.7±1.2 | 77.7±1.6 | 0.838±0.020 | 0.857 | 0.968 |

## T_b (K km/s) and 4-3/3-2 at 1.5" (Xb4b 3-2; Xb4d in brackets)
| | 3-2 | 4-3 raw | 4-3/3-2 raw | 4-3/3-2 (P) | 4-3/3-2 (F) |
|---|---|---|---|---|---|
| HCN N | 75.0±0.6 (76.5) | 45.7±0.6 | 0.608±0.009 (0.597) | 0.623 | 0.695 |
| HCN S | 21.9±0.6 (25.0) | 9.06±0.56 | 0.413±0.028 (0.362) | 0.488 | 0.680 |
| HCO+ N | 98.8±0.7 (98.2) | 63.6±1.1 | 0.644±0.012 (0.648) | 0.683 | 0.837 |
| HCO+ S | 58.6±0.7 (63.7) | 48.6±1.1 | 0.830±0.021 (0.763) | 0.860 | 0.958 |

## HCN/HCO+ and N/S contrast
| beam | 3-2 N / S | contrast 3-2 | 4-3 raw N / S | contrast raw | contrast (P) | contrast (F) |
|---|---|---|---|---|---|---|
| 1.0" | 0.836 / 0.395 | 2.12 | 0.814 / 0.222 | 3.67 | 3.34 | 2.26 |
| 1.5" | 0.759 / 0.374 (Xb4d 0.779 / 0.392) | 2.03 (Xb4d 1.99) | 0.718 / 0.186 | 3.85 | 3.26 | 2.37 |

- **Contrast ratio 4-3/3-2.** At 1.0" it is 1.73 raw, 1.58 (P) and 1.07 (F). At 1.5" it is 1.90 raw, 1.61 (P) and 1.17 (F).
- **Where the 3.3 comes from.** The "3.3 after missing-flux correction" quoted in p0053 is the pedestal (P) correction, reproduced here as 3.34.

## Does the contrast rise with J?
- **Uncorrected.** Yes, by a factor 1.7 to 1.9 (2.1 to 3.7 at 1.0"). The rise comes from S HCN. Its 4-3/3-2 is 0.47 against 0.65 at N, while HCO+ is 0.84 at S against 0.66 at N.
- **Sensitivity to the correction.** The S HCN recovery of 0.61 and the HCO+ recoveries of 0.77 to 0.87 are exactly the size needed to erase the rise under (F). Under (P) about 0.9 of the rise (1.6 of 1.7 in ratio terms) survives.
- **Independent checks that do not settle it.**
  - Xb47 alone at 4.5" gives N/S 0.515 / 0.317, a contrast of 1.62. The same pixels in the high-resolution data give 0.587 / 0.223, a contrast of 2.64. The 4.5" beam also blends the extended S disc.
  - Two-component uv-plane fits (p0053/Q60, `/workspace/out/q60/summary.csv`, hw350) give 4-3 contrasts of 1.4 (00493+00223, 2 components) to 1.6 (Xb47), and up to 1.6 to 1.9 with 3 components. These are not core-aperture quantities. The model systematics were never closed, and the 1.5× absolute flux offset between Xb47 and 00493/00223 (Q60 round 2) is unexplained.
- **Conclusion.** The 4-3 contrast at ≤1" is bounded between 2.3 and 3.7, central about 3.3 (pedestal). The 3-2 contrast is 2.0 to 2.1 (about 2.0 using Xb4d). The ratio of the two is 1.1 to 1.8 (central 1.6). A J-rise of ≥1.5 is plausible but not established. J-independence (ratio 1.0 to 1.1) is allowed only by the maximal-correction case (F).
- **What would settle it.** A 4-3 pair at ≤1" with short spacings (a 12m+7m combination, or a one-tuning HCN+HCO+ 4-3 observation at ≤1"). The 3-2 side also needs the 12m+7m taper at 1" (open in p0036).

## Caveats
- Statistical errors only. The per-project flux scale (about 10% each, so 14% on the 4-3 ratio) is not included, and it would not cancel in the contrast between N and S of the same pair.
- The recoveries come from a 4.5" comparison in which Xb47 has a 1.5× absolute flux offset relative to 00493/00223 (Q60). If the offset is a Xb47 scale error, the (P) and (F) corrections are overestimated and the raw contrast (3.7) is closer to the truth.
