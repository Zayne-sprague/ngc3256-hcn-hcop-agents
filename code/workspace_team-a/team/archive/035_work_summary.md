# Q60 work summary (worker, round 2)

## Question
N/S HCN and HCO+ 4-3 contrast: 3.7 at 1" (2018.1.00493 + 00223 pipeline images) vs 1.6 (single-tuning Xb47). Which is right?

## What was done
Raw EBs restored (00493: 3 ACA EBs, a 3-pointing mosaic; 00223: 2 EBs; Xb47: 8 EBs; 12m 01246: 1 EB).
Not done: the tclean tapered per-EB imaging. The reviewer's visibility-dump crash (u,v 1-D) was fixed by not using getdata for u,v. Instead I fit the line-integrated, continuum-subtracted visibilities directly in the uv plane.
- `scripts/q60a_imgplane.py` image-plane integrals at 4.5/5.5/7" (the earlier result, `out/q60/imgplane.csv`).
- `scripts/q60c_uvspec.py` and `q60f_mosaic.py` (CASA, `rawdata.py casa`): per-EB, per-channel linear least-squares fit of two compact Gaussian components (N, S) (+ optional extended 8" component) to the visibilities in uv sets all / 8-15 / 15-25 / 25-50 / >50 klambda.
  - q60c used one phase centre. This decorrelates the 00493 mosaic.
  - q60f redoes the mosaic EBs with a per-field phase centre and a Gaussian primary-beam gain per field.
- `q60d_analyse.py` and `q60e_tables.py`: baseline fit on line-free channels (470<|dv|<edge km/s), integrals for +-250/350/500 km/s, curve of growth (`out/q60/uvfit.csv`, `cog.csv`).
- `q60g_summary.py`: per-project medians and scatter (`out/q60/summary.csv`).
- `q60h_joint.py`: feeds the new R43 into `q55e_joint` (`out/q60/joint43.txt`, `.json`).

## Results (uv plane, all baselines, +-350 km/s, 2 compact components)
| project | HCN N/S | HCN/HCO+ (Tb) N | HCN/HCO+ (Tb) S | contrast |
|---|---|---|---|---|
| Xb47 (8 EBs) | 2.13+-0.28 | 0.64+-0.07 | 0.39+-0.08 | 1.63+-0.27 |
| 00493 HCN x 00223 HCO+ (3x2 EBs) | 1.97+-0.21 | 0.71+-0.10 | 0.52+-0.05 | 1.37 |

(+-: sd between EBs.) The 12m EB X3cae has HCN N/S 1.97.

- The contrast is 1.3-1.9: 1.6 compact-only vs 1.9 with an extended component (Xb47), over window half-widths 250-500 km/s and over the uv sets.
- It rises with uv radius, to ~2.0+-0.5 at 25-50 klambda.
- **Adopted: 1.5+-0.3. The 3.7 is not reproduced.**
- The HCN N/S from 00493 (2.0) agrees with Xb47 (2.1) in the same fit. The "S HCN 0.61 of Xb47" deficit is therefore an image-plane/deconvolution (ACA mosaic X1df1 vs Xb47) effect, not a per-EB flux-scale difference.
- A pure flux-scale error is excluded as the cause of the contrast: it cancels in ratios.
- Fluxes of 00493/00223 are ~1.5x those of Xb47 for both nuclei and both lines. This is a common factor (calibration/model/PB), not differential.
- Xb47 EB X8300 gives an inverted HCN N/S (0.5; N flux 2-5x low). It was treated as failed and excluded; HCO+ in the same EB is normal.

## Q55 effect (ln BF, from `joint43.txt`)
- "Identical physics" stays excluded (lnBF <= -45) for every input.
- 2.5" set: X-only vs column+abundance differ by <0.5 (old 3.7 input: column+abundance won by 14.5 over X-only; Xb47 image input: X-only by 1.1).
- 1" core set: excitation (T,n) is best by 0.8-3 over column+abundance. With the old input n-only was best, col+ab at -0.08, so this is not a Q60 effect.
- Net: the 4-3 input no longer favours column+abundance. Column/abundance vs excitation stays undecided (|dlnBF|<3).
- Xb47 is the more trustworthy set: 8 EBs, known scatter, and HCN agrees with 00493.

## Not tested / uncertain
- No flux-calibrator audit and no per-EB continuum N/S check. The split MS contain the target only, so the calibrator flux tables are not available. This was not done and is not excluded.
- Fixed-position compact-component model: it does not equal clean-image apertures. Model-dependent by +-0.3 in the contrast.
- Why the pipeline-image ratio is 3.7 is not proven. Likely causes: HCN and HCO+ from different OUS with different uv coverage and deconvolution, with the S emission resolved out in the 00493 mosaic.
- The cog.csv spw-edge check is only partly done. The 12m EB has no negative-side baseline.
- No literature used for this item beyond Harada et al. 2018 (1.7" nuclear ratios), already cited in earlier rounds.
