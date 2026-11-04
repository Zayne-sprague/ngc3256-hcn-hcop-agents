# Work summary: Q4 (opacity, N-core HCN excess, HCN-vib, absorption), round 2 (reviewer fixes)

All scripts in /workspace/scripts, outputs in /workspace/out/q4r2 (round-1 outputs s41-s46 unchanged). Report updated by `scripts/s51_update_report_q4r2.py` (idempotent).

## Fix 1: width-ratio / R(v) null test  (`s47_width_null.py` -> out/q4r2/s47_width_null.txt/json)
- Xb4b HCN 3-2 and HCO+ 3-2 both convolved channel-wise to a common circular 0.30" beam; N r<0.25" aperture.
- FWHM HCN 149.9, HCO+ 138.2 km/s, ratio **1.085 +- 0.015** (400 noise realisations on the observed spectra) +- 0.013 (pixel bootstrap); native beams 1.091 +- 0.018. Centroid HCN-HCO+ -3.2 +- 0.8 km/s (noise). Outer annuli 0.25-0.5" and 0.5-1": 1.038 +- 0.010, 1.038 +- 0.009.
- Null from 190 beam-spaced r<0.15" apertures (S/N>=8 in both lines, r<1.2" of N): ratio median 1.036, 16-84% 0.98-1.09, 17% of apertures >= 1.084; centroid difference 16-84% -1.0..+1.5 km/s (inner 13 apertures -5..+7.5). The core value is ~1 sigma of the real (profile + spatial) scatter, which is ~3x the statistical error. **Not significant; the "hint of opacity" is withdrawn** (R(v) is the same information).
- Spatial-weights model (HCN := HCO+ cube x smoothed per-pixel R map, identical local profiles; velocity gradient 119 km/s/arcsec = 36 km/s over the beam): predicts ratio 1.003, centroid -0.1 km/s, flat R(v). The reviewer's explanation (more compact HCN + gradient) therefore goes the wrong way (HCN would be narrower) and does not reproduce 1.085 either; the observed excess is just within the field scatter.
- Opacity expectation (Gaussian-profile broadening): FWHM/thin = 1.02, 1.18, 1.49, 1.96, 2.33 at tau 0.1, 1, 3, 10, 30. Core lines are narrower (138/150 km/s) than the 0.25-1" annuli (153-164), opposite to saturation broadening.

## Fix 2: T_b discriminator  (`s48_tb_discriminator.py` -> out/q4r2/s48_tb_discriminator.txt/json, s48_allowed_N/S.csv)
- Observed T_b(HCN)/T_b(HCO+) at 0.25": N 1.461/1.556 = 0.94; S 0.792/1.869 = 0.42.
- RADEX grid (6925 cells; shared Tkin, n, filling; independent N(HCN), N(HCO+)), accept +-15% on the ratio and eta<=1. Single-Tex slab: ratio -> 1.00 for tau_O>=5 (1.05 at tau_O=3).
- N core (29% of grid allowed): both-thick (tau>=3): log N(HCN)/150 km/s 15.75-19 (median 17.75), eta 0.01-0.75 (median 0.07), Tex(HCO+) 7-150 K (median 30), A unconstrained (0.3-18). Both-thin: log N 14.25-16.5, eta 0.02-0.99, A = 1.0-5.6 (median 1.8, close to the LTE value 1.6), HCN sub-thermal (TexH/TexO 0.35-0.96). The thin branch needs no coincidence, only A~1-6. **Not separable** with these data.
- S core (4.4% allowed): in every branch Tex(HCN)/Tex(HCO+) ~0.5-0.6 or HCN thin with A<1 (69% of thin cells, 37% of thick cells); equal-Tex saturated lines are excluded. Which line is thicker at S is not decided by T_b alone (both HCN-thick/HCO+-not and HCO+-thick/HCN-thin cells exist). Thin-branch A: N 1.8 vs S 0.56 (factor ~3).

## Fix 3: H13CN injection-recovery  (`s49_h13cn_injection.py` -> out/q4r2/s49_h13cn_injection.txt, s49_injection.csv, s49_coverage.json)
- Geometry (correcting the reviewer's premise): the band edge 342.396 GHz is at v=-224 km/s, on the **blue** (high-frequency) side, so only 2.8% of an HCN-like profile is outside the band. The CO 3-2 core, not the band edge, limits the fit windows: 341.80-342.20 GHz (v -54..+293) holds 75% of the template flux (Fisher information 33% from v<0, 67% from v>0); 341.85-342.12: 42%; 341.80-342.28: 96% (CO contaminated, few families fit). The v>0 half alone (341.80-342.14) gives alpha_3s<0.12 for all families.
- Injection into the real spectrum: recovered (bias <=0.01 in the main window). Simulations (wing truth from real fit in 3 families + template + noise segments from the line-free spectrum; 144 per level and window, 4 recovery families): bias +0.013 +-0.020 (alpha_inj 0), +0.013 (0.03), -0.003 (0.1), +0.007 (0.2), +0.009 (0.3); **false-limit rate 0/144 at every level** (conservative); >=3 sigma detections 0, 26, 85, 95, 100%; median 3 sigma UL on pure noise 0.06 (p84 0.10, max 0.12). Narrow window is worse (bias up to +0.04+-0.05).
- Real N: main-window best alpha 0.02-0.04 (1-2.9 sigma), narrow window 0.08-0.10 (s41: 0.09, a1 0.035-0.095) = HCN/H13CN ~10-30 **marginal**, wing-degenerate, not a detection; limit alpha<0.075-0.12 (HCN/H13CN>8.3).
- Sensitivity: thin 12C/13C=30-100 means alpha=0.010-0.033, below the 3 sigma level (0.06). The data only detect tau(HCN 4-3)>~2-3 (X=30), 4-5.5 (X=50), 8-11 (X=100) for alpha=0.075-0.105, or 12C/13C<~10-17.

## Fix 4: core fraction band  (`s50_core_fraction.py` -> out/q4r2/s50_core_fraction.txt, s50_core_bound.csv)
- f3(r<0.5" / r<1.7") = 0.36 with Xb4b, 0.34 against the Xb4d total (17.0 Jy km/s); scanned 0.15-0.50.
- L_tot = Tb(4-3)/Tb(3-2) at 6.5": HCN 0.50 (+-15% scale), HCO+ 0.64. For a thick core (L_c=1) f4 = f3 L_c/L_tot = 0.40, 0.73, 1.0 for f3=0.2, 0.36, 0.5, i.e. larger than f3 (the round-1 assumption f4=f3 was too lenient).
- tau_core(HCN 4-3) upper limit (thick core, L_H13=9.5): X=30/50/100: 7/13/28 (f3=0.2), 4/7/15 (0.36), 3/5.5/11 (0.5). Over the whole band tau_max = 3-137 (median 10-19); a tau>=10 core is allowed in 52-79% of scenarios (22% at X=30, 63% at X=50, 100% at X=100). Thick core neither excluded nor required.

## Uncertain / remaining
- Core opacity vs abundance is still undecided (no H13CN/H13CO+/HNC at <0.5"); T_b, width and ladder constraints are satisfied by both branches.
- Single-zone RADEX with shared filling; grid-weighted fractions are not posteriors. Injection noise has no CO-wing structure; alpha grid step 0.02. L_tot comes from a 6.5" beam-matched pair (Band 6 vs 7 scale 15%).
- Report: /output/report.json updated (4 Q4 entries revised, 3 added; interpretation, caveats, recipe step 13).
