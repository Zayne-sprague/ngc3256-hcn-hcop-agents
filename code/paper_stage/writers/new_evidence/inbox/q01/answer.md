# q01 — dilution test of the compact N/S contrast

Script: `q01_compute.py` (reuses `/workspace/scripts/s410_q64_imageplane_calib.py` loader/moment/convolution/annulus code verbatim; aperture estimator = `/workspace/scripts/s192_q34_measure.py` recipe), `q01_plot.py`; numbers in `q01_values.json`. Data: Xb4b Band 6 pbcor cubes spw27 (HCN) / spw25 (HCO+), native beam 0.246"x0.237". Assumes 194 pc/arcsec. The convolution is in the image plane (flux-conserving kernel), same baseline (med, line-free |v-2780| 600-1500 km/s) and window (±350 km/s) at every beam.

## N/S (HCN/HCO+) vs beam, Xb4b only (stat = null scatter of 5 off-line pseudo-moments)
| beam | core r<0.25" (Table C1) | r<1" (Table C1) | aperture r=beam FWHM (uvtaper estimator) |
|---|---|---|---|
| native 0.25 | 2.27±0.11 | 2.03±0.10 | 2.23±0.11 |
| 0.3 | 2.26±0.11 | 2.03±0.10 | 2.25±0.13 |
| 0.6 | 2.21±0.11 | 2.02±0.09 | 2.11±0.11 |
| 1.05 | 2.09±0.10 | 1.98±0.08 | 1.98±0.08 |
| 1.5 | 2.01±0.08 | 1.94±0.07 | 1.85±0.07 |
| 2.2 | 1.89±0.07 | 1.85±0.07 | 1.66±0.07 |

Check: reproduces p0064 (2.26 at 0.3" r<0.25"; 1.98 at 1.05" r<1"; my stat errors 0.11/0.08 vs 0.09/0.06 published because null scatter is here propagated per-beam). Reference values: 12m+7m uvtaper 1.532±0.026 (1.5") and 1.441±0.023 (2.2"); 12m-only 1.578, 1.451; 1-0 1.616, 1.467 (s192/p0036).

## 1.5" comparison: where does 2.03 (team-a) vs 1.53/1.58 come from?
Xb4b image plane at 1.5" with the uvtaper aperture: **1.85±0.07**. Steps (all at 1.5"):
- team-a 2.03 → Table C1 core annulus 2.01 and r<1" annulus 1.94: estimator (annulus vs median of windows) accounts for ≈0 to −0.1; consistent within errors.
- annulus → uvtaper aperture on same Xb4b image: 1.94 (2.01) → 1.85: aperture/centering, −0.09 (−0.16).
- Xb4b image-plane aperture 1.85 → 12m-only uvtaper 1.58: −0.27 (≈ 3.5σ). Same estimator and aperture, so this is the array/uv-matching term: image-plane convolution of the 0.25" data vs a tapered 1.5" image from the Xb4d-project 12-m EB (Xe07), different uv coverage/weighting for the two lines. The same EB set's Xb4d image convolved to 1.5" gives 1.81 (s192), so a 1.81→1.58 shift also occurs within the compact-config data: it is the convolution-vs-taper (uv weighting) effect, not the Xb4b array itself. I did not taper Xb4b in the uv plane, so I cannot separate Xb4b-vs-Xe07 EB from convolution-vs-taper.
- adding 7-m: 1.58 → 1.53, −0.05 (≈1.4σ; small).
So at 1.5": estimator ~0–0.1, aperture ~0.1, array/uv-matching ~0.27, 7-m ~0.05 of the 2.03→1.53 drop of 0.50.

## Fraction of the drop explained by dilution alone
Drop = 2.26 (0.3", Table C1 core) → 1.44 (2.2", 12m+7m uvtaper) = 0.82.
- Table C1 core estimator (r<0.25" convolved to 2.2"): 2.26→1.89, i.e. 0.37 = **~45%** (±15%).
- uvtaper aperture estimator on Xb4b convolved: 2.25→1.66, 0.59 of 0.81 = **~73%** (±16%).
- r<1" annulus: 2.03→1.85 (only 0.18).
Dilution alone therefore explains roughly half to three-quarters (45–73%, depends on the estimator; the aperture estimator is the like-for-like one for the 1.44) of the drop. The remainder, 1.66→1.44 (≈0.2, ≈3σ; 27–55% of the drop), is the uv-matched (taper) vs image-plane-convolution/flux-recovery term, not dilution. Contrast does not fall to ≈1 with dilution: at 2.2" the flux-complete value (1.44) coincides with 1-0 (1.47).
Caveats: errors are statistical (null scatter); the fraction errors are approximate; at ≥1.5" the N and S apertures overlap (1.45" separation), which dilutes the contrast in any estimator.
