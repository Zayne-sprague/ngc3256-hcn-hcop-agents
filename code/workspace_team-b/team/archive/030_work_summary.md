# Q56 round 2 work summary (worker, team-b)

## Question
Which is right for HCN R32/10 at N vs S: Q53 (S/N 0.71), Q34/"equal" (S/N 1.007), and why do they differ? Reviewer asked for systematics, uv-matching, non-LTE bound, reconciliation.

## Scripts / outputs
- scripts/s344_q56_systematics.py -> out/q56/s344.json: common 0.1" grid moment maps (1-0 X961 convolved to 2.5/3.5"; 3-2 12m+7m comb cubes), variants in beam, window (100/150/350), baseline (lin/const), estimator (pixel, HCO+-peak pixel, aperture 0.5 theta, aperture theta), 3x3 jackknife, image-plane uv high-pass at 0/7.7/13 klambda.
- scripts/s345_q56_summary.py -> out/q56/s345.txt, s345.json, s345_budget.json
- scripts/s346_q56_radex_bound.py -> out/q56/s346.txt (project RADEX grid, excitation-only vs column pairs)
- scripts/s347_update_report_q56r2.py -> /output/report.json (246 results; 2 old Q56 entries replaced by 3 final ones; Q53 and Q29 entries annotated). Backup report_backup_q56r2.json.

## Results
1. Reproduction of the s342 numbers on the new grid (HCN 2.5" N 0.306, S 0.277; HCO+ 0.400/0.421).
2. Error budget. HCN S/N over 48 variants (full uv) 0.91+-0.04 (0.83-0.98); with uv filters 0.95+-0.05; all 144: 0.93+-0.05 -> 1.4 sigma below 1, NOT significantly different from 1. HCO+ S/N 1.00+-0.05. Position jackknife: S R32/10 varies +-10% (N +-1.5%), S/N +-0.06-0.09. Flux calibration is common to N,S and cancels. D=contrast32/contrast10 = 1.07+-0.05 (range 0.97-1.20).
3. Reconciliation: Q53's 0.71 = uv-incomplete Xb4d 12m image (S 3-2 low by 25% HCN / 16% HCO+, N 0-4%). Q34's 1.007+-0.020 was an aperture-mean (r=FWHM) estimator with statistical-only error; aperture means dilute the compact N core (pixel R32/10 0.31 vs aperture 0.27) while S is unchanged, so the estimator choice moves S/N between 0.83 and 1.06. Both "equal" and "S 5-15% lower" lie inside this spread; neither is demonstrated. The N/S HCN/HCO+ contrast itself depends on aperture (1-0: 1.53 pixel -> 1.29-1.40 at r=FWHM; 3-2 1.69 -> 1.32-1.48).
4. uv-matching: raw X961 1-0 visibilities could not be restored (earlier request failed), so I used an image-plane high-pass of all maps (>7.7 and >13 klambda). This raises S R32/10 (HCN S/N 0.91 -> 0.95-0.96) and leaves D at 1.07; the Q34 12m-only (c12) vs comb test gave 4-15% flux changes. A true visibility-level matched-uv 1-0 imaging is NOT done.
5. Non-LTE: R32/10 S/N values alone do not exclude excitation-only (Tk,n) origin of the 1-0 contrast (~7% of consistent same-A/N pairs reach contrast >=1.35); they are a weak discriminator, as Q34 ln BF and team-a p0059 found. The old "excitation explains at most ~10%" is replaced: D-1 = 0.07+-0.05 (<0.17, 95%) assuming equal opacity/abundance.
6. Forum: consistent with p0058-r01/p0059 (1-0 N/S 1.45+-0.22 at 1", 1.44-1.58 at 2.5"). H13CN/HCN at ~2" remains the decisive test (not attempted this item).

## Literature
Harada+2018 (Table 13 via p0058), Michiyama+2018, team-a p0059; no new papers read.

## Uncertain
Image-plane uv filter is approximate; S centroid sensitivity (3-2 vs 1-0 peak offsets) dominates the S/N error; RADEX grid subsampled with flat priors and no beam filling/absolute intensities; R43/32 and S non-dust continuum (reviewer's questions) not addressed.
