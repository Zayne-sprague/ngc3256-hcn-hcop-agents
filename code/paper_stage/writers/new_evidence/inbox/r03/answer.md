# r03 — adopted core N/S HCN/HCO+ 3-2 contrast (team-a)

## Adopted values (image plane, Xb4b 12-m; Xb4d 12-m+ACA used for the 1" cross-check)
| beam | N/S (adopted) | total uncertainty |
|---|---|---|
| 0.3" (r<0.25", ~55 pc) | **2.2** | **+0.16/-0.15 (quote ±0.2)** |
| 1.0" (~190 pc) | **2.1** | **+0.12/-0.15 (quote ±0.15)** |

Basis. At 0.3", team-a medians over baseline windows are 2.13 [1.92-2.25] (r~0.1") and 2.21 [2.06-2.31] at 0.4"; team-b r<0.25" is 2.26. At 1.0", team-a Xb4b is 2.14 [2.09-2.18] and Xb4d 2.15; team-b is 1.98 (r<1") and 2.05-2.06 (r<0.5"). The differences are aperture (r<1" includes more extended gas) and the window convention, and they are within the budget. (p0065; out/q64r2/summary.txt, profile.txt, img_ratios.csv; team-b p0064, /snapshots/team-b-s35/out/q64r2/final_budget.txt.) The 1" value is a rounded blend of the two teams. I did not run a new fit.

## (a) Error components (quadrature, `r03/r03_budget.py` -> `r03_budget.txt`)
| component | 0.3" | 1.0" | source |
|---|---|---|---|
| stat (null apertures, line-free integrals) | 0.09 | 0.06 | team-b final_budget (0.09/0.06); team-a null.txt, summary.txt noise range |
| baseline / window convention | 0.07 (±0.03 on N=0.97, i.e. 3% of the ratio) | 0.01 | p0065: ±0.03 on N at ≤0.4", ≤0.01 at ≥1"; team-b 0.05 |
| array / flux recovery (Xb4b vs Xb4d at 1.05", ratio of ratios 0.98-0.99) | 0.05 (≤2%, carried inward from the 1"->0.6" constancy; not measured below 1") | 0.04 | summary.txt; p0064 |
| residual/restored (tclean, N 2-4%, S 5-8%) | 0 / -0.08 | 0 / -0.08 | team-b s414 (beamcheck_resid.txt). Team-a did not measure this. |
| uv-model form spread | +0.10 (allowance; see (d)) | ±0.10 | uvfit3 A vs image: 2.04-2.20 vs 2.14-2.26; see (c) |
Flux-recovery term (array + residual) is therefore ±0.04-0.05 up, -0.09 down, and it is not the dominant term. The sum is 0.16 (0.3") and 0.15 (1.0"). The 0.3" total is dominated by stat and baseline. At 0.3" I quote ±0.2 rather than 0.16 because the window extremes reach 1.92 and the Xb4b visibilities needed for a true recovery test do not exist (4/4 restores OOM).

## (b) Largest N/S bias allowed at <1"
- Measured differential flux loss: ≤2% at 1.05" (N/S bias ≤ +0.04), from the Xb4b/Xb4d ratio of ratios 0.98/0.99 (both teams). The bias at <1" is inferred from the 1"->0.6" constancy, not measured.
- Our own uv-model test (Xb4d+ACA, 6-component, baseline A; out/q64r2/uvfit3_A.log) gives N/S 2.04 at 1.0" and 2.20 at 0.6", i.e. at most -0.1 below the image. The uv models can raise it only with the rejected baseline B (2.35/2.65, chi2 47.7 vs 5.5).
- So the largest N/S shift we allow is about -0.1/+0.2 relative to the image-plane value (a 5-10% effect), and it does not change the ≥2 contrast.
- Team-b's high-pass bound (2.6-3.1 at uvmin 60-100 kλ): I would not use it as an allowed bias. It removes all flux below uvmin by construction, whereas the real data (Xb4d + ACA) contain the short spacings, and the measured loss is ≤2%. Team-a did not reproduce that test. I accept it only as a pessimistic bound on what a purely long-baseline measurement would give, and as a qualitative sign that the ratio rises when extended gas is filtered out. It is not strictly excluded above ~100 kλ (team-b).

## (c) uv fit vs image plane at 0.5-1"
Yes, within ~0.1 for the adopted baseline (Xb4d+ACA 6-component uv model, baseline A, redchi2 5.48; uvfit3_A.log): N/S 2.20 (N 0.96, S 0.44) at 0.6" and 2.04 (0.85/0.42) at 1.0", against image 2.26 and 2.14 (Xb4b), 2.15 (Xb4d, 1"). At 1.5" the uv model gives 1.90 vs image 2.03 (Xb4b)/1.99 (Xb4d). Below the Xb4d resolution (0.94x0.74") the uv model is only an extrapolation: 2.44 at 0.4" and 3.16 at 0.25" (model-form spread N(0.25") 0.98-1.5). That disagrees with the Xb4b image (2.1-2.3), so the uv fit is not used below 0.6". The Xb4b visibilities could not be restored, so there is no uv confirmation at 0.3". The image-plane statement stays "not uv-confirmed" there.

## (d) The 1.3-1.5 long-baseline N amplitude excess
Not explained. The N HCN/HCO+ amplitude ratio in annuli of uv distance is 0.59-0.92 at 3-90 kλ but 1.29±0.08 at 90-130 kλ and 1.49±0.28 at 130-200 kλ (out/q64/uvannuli.log; the same with and without a Bridge component). Actions: (1) a sixth, compact N component with free per-line flux and size (q64i_uvfit3.py) was tried (fitted HCN-only 0.05" component, 0.35-0.51 mJy in uvfit3_A/B.log); I did not verify that it removes the excess, and it gives the N(0.25") model-form spread 0.98-1.5; (2) the Xb4b image (N native 1.00 [0.92,1.06]) excludes the high end; (3) we withdrew the earlier uv-based N excess claim and do not use uv-derived numbers for the adopted contrast; (4) it enters the budget only as the +0.10 uv-form allowance at 0.3". Likely causes (unproven): non-Gaussian N shape or annulus blending. If real, it would push the true N core ratio and N/S up, so the adopted values would be conservative on the high side only.
