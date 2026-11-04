**Short answer.** On a matched ACA-only footing the 3-2 and 4-3 N/S contrasts are the same within errors. At the coarsest common circular beam (6.5") N/S(4-3)/N/S(3-2) = 1.02 ± 0.18. At 5" it is 1.07 ± 0.17. N+S blending does not explain the 4-3 vs 3-2 difference, and it is not needed to make the two agree.

Source: `scripts/s438_q66_32_vs_43.py` → `out/q66s/uvmatched_32_43.csv`. This is team-b's p0066-follow-up (q66r2) table. `r09_table.py` and `r09_table.txt` in this folder are a read-out of that CSV, with no new fits.

**Recipe (identical for both lines).**
- ACA only: Xb4f for 3-2 (12-m+7-m data, 7-m-only cubes via the c32 images) and Xb47 for 4-3.
- Same momap code (`s431_q66_maps.py`), same circular restoring beam, same N and S positions and apertures, same window and baseline.
- uv cuts are scaled by the frequency ratio.
- N/S = (HCN/HCO+)_N / (HCN/HCO+)_S, in Jy km/s moments with the same ν² factor at N and S, so it cancels.
- Errors are 16–84% intervals from a null-aperture Monte Carlo (≥30 empty positions, 1.4826·MAD). I quote the half-width.

| beam / uv cut (4-3 ↔ 3-2) | N/S 3-2 | N/S 4-3 | 4-3/3-2 |
|---|---|---|---|
| 6.5", native uv (coarsest common beam) | 1.23 ± 0.16 | 1.25 ± 0.14 | **1.02 ± 0.18** |
| 5.0", 4-3 all uv ↔ 3-2 >8 kl | 1.36 ± 0.18 | 1.45 ± 0.13 | 1.07 ± 0.17 |
| 5.0", >15 kl ↔ >11.3 kl | 1.40 ± 0.19 | 1.60 ± 0.17 | 1.14 ± 0.19 |
| 5.0", >20 kl ↔ >15 kl | 1.56 ± 0.33 | 1.84 ± 0.26 | 1.18 ± 0.30 |

The ratio errors combine the two Monte Carlo errors in quadrature. They treat the two lines as independent, which is approximate because both use the same N and S positions.

**Individual lines (raw N/S flux ratios), 6.5" / 5" all-uv.**
- HCN: 3-2 1.23 / 1.48, and 4-3 1.38 / 1.76.
- HCO+: 3-2 1.00 / 1.09, and 4-3 1.10 / 1.21.
- The HCN-driven contrast is slightly larger at 4-3, and the HCO+ ratio is also slightly larger.
- Per-line errors are not tabulated. The contrast ratios above are the quantities with errors.

**Consistent with 1?**
- Yes at the native beams (6.5" and 5"), within about 0.4σ.
- The ratio rises to 1.14 and 1.18 at the higher uv cuts. These are still within about 1σ, and the points are not independent because they share the same data.
- All four pairs lie above 1. The sign is the same as the 12-m resolution trend in p0066, but the ACA data cannot establish it. The data allow 1.0–1.2.

**Does N+S blending explain any difference?**
- No. I tested this with the deblend models A, B and C from `s437`. They are the same models for both lines (two Gaussians at N and S, plus a constant, plus a broad component in B or a bridge component in C), and the results are in `out/q66s/deblend.csv`.
- Blending dilutes both contrasts. At 6.5" the deblended 3-2 values are 1.34 / 1.72 / 1.36 for models A / B / C, up from 1.23. The deblended 4-3 values are 1.24 / 1.34 / 1.26, up from 1.25.
- The deblended 4-3/3-2 ratio at 6.5" is 0.93 / 0.78 / 0.93. At 5" it is 1.04 / 0.97 / 1.01.
- So blending changes the absolute N/S by 10–40% in 3-2, more than in 4-3 under model B. It does not create a 4-3 excess. If anything, deblending moves the ratio slightly below 1 at 6.5". The finding that the ratio is consistent with 1 holds with and without deblending.
- At the higher uv cuts the deblended ratio is 1.0–1.3, which again is not a clear departure from 1.

**Why the numbers in the request differ.**
- The request quotes ACA 4-3 N/S of 1.62 / 1.49 / 1.32 at 4.5 / 5.5 / 7" (p0066, `out/q60/imgplane.csv`) and ACA 3-2 N/S of 1.07 at 6.5" (p0031).
- Those came from different recipes: earlier q60 and p0031 apertures and baselines. They are not matched to each other, and I did not reproduce them here.
- The matched recipe gives 1.25 at 6.5" for 4-3, consistent with the q60 trend (1.32 at 7"). It gives 1.23 for 3-2, higher than 1.07.
- I do not have an explanation for the 3-2 difference (1.23 vs 1.07) without re-running p0031, and I did not do that. The 1.07 should not be paired with the 4-3 numbers; the matched 1.23 should be used.
- If the paper uses the unmatched 1.07 against 1.32, the apparent 4-3/3-2 ratio would be about 1.23 and would be a recipe artefact.

**Caveats.**
- These are ACA-visible fluxes at 5–6.5", not small-scale contrasts, so they say nothing about the 1" J dependence in R2 (D = 1.00 ± 0.05, 1-0 vs 3-2 at 2.2").
- The 3-2 data at 5" (Xb4f 7-m-only maps) have lower S/N than 4-3, which is why their errors are larger.
- The 2018.1.00493/00223 (ACA18) rows are not used. S HCN is low there, and the b5u20 and b6p5 rows are unreliable.
