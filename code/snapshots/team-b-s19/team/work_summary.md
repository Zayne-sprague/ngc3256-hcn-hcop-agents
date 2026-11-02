# Q26 work summary (worker, round 1)

## Question
HCN/HCO+ 3-2 ratio at native 0.24" (50 pc) and matched 55/120/250 pc beams for N, S, Bridge, E1, E2 and the outflow lobes, with velocity-resolved maps. Is p0026's result (HCO+ T_b equal at N and S, HCN 2x lower at S) a significant ratio contrast, and does it depend on scale?

## What I did
- `scripts/s170_q26_ratios.py [Xb4b|Xb4d|Xb4f]`: pbcor cubes -> per-channel RJ T_b -> flux-conserving convolution to circular beams (0.26" common beam for both lines = 50 pc; 0.2835, 0.6186, 1.2887" = 55/120/250 pc; 6.5" for Xb4d/Xb4f) -> per-pixel baseline from line-free channels -> moment-0 and spectra. Regions as in Q6, plus the S blue lobe and red ridge (Q17 positions and windows). Ratio estimators: aperture (r = beam FWHM) at the HCO+ peak, same pixel, peak/peak, 3-channel peak T_b. Noise comes from the line-free channels of the aperture-mean spectrum. Also produced: >5 sigma masked ratio maps and 3 sigma blue/red wing maps (`out/q26/maps_*.npz`).
- `scripts/s171_q26_analysis.py`: tables, N/S contrast MC, power-law scale dependence, wing/core, 50 km/s bins, Xb4d/Xb4f consistency, blue-only baseline test, figure `out/q26/q26_ratio_figure.png`. Output: `out/q26/s171_summary.txt` (Tables 1-9) and `s170_*.json`.
- `scripts/s172_update_report_q26.py`: adds 3 results and 1 recipe step to /output/report.json (valid against the schema; backup `report_backup_q26.json`).
- **Bug found and fixed on the way:** the HCN (spw27) and HCO+ (spw25) cubes have velocity axes offset by ~117 km/s (and 4.40 vs 4.38 km/s channels). A first run applied the HCO+ axis to HCN. Moment-0 ratios barely changed, but the wing/core and bin ratios were garbage. Now HCN is interpolated onto the HCO+ grid. Any earlier team script that shares one velocity vector between the two lines should be checked (s160 and the Q17 library use per-line axes, so are fine).

## Results (T_b ratio HCN/HCO+, aperture; 14% flux scale not included in the stat errors)
| | 50 pc | 55 pc | 120 pc | 250 pc | Xb4d 250 pc | 6.5" Xb4d / Xb4f |
|---|---|---|---|---|---|---|
| N | 0.99 | 0.99 | 0.84 | 0.69 | 0.71 | 0.43 / 0.46 |
| S | 0.44 | 0.43 | 0.40 | 0.37 | 0.39 | 0.39 / 0.43 |
| Bridge | 0.57 | 0.56 | 0.50 | 0.48 | 0.46 | 0.42 / 0.45 |
| E1 | 0.21 | 0.24 | 0.24 | 0.14 | 0.26 | |
| E2 | 0.42 | 0.43 | 0.50 | 0.52 | 0.28 | |
| Sblue | 0.84 | 0.84 | 0.77 | 0.85 | | |

Same-pixel ratio at the HCO+ peak: N 1.01, S 0.48 at 50 pc. Stat errors are 0.01-0.03 for N, S and Bridge.
- **N/S contrast = 2.1 +- 0.2 at 50-120 pc** (aperture 2.25, 2.32, 2.11; same-pixel 2.09, 2.14, 2.24), 1.9 at 250 pc (Xb4d 1.8). It is significant (>40 sigma stat). The flux scale cancels in it, and the line-free baseline choice changes it by <1%. The estimators spread it over 1.9-2.3. Peak T_b N/S at 50 pc: HCO+ 0.95 (equal), HCN 2.03, so p0026 is confirmed: the contrast is entirely in HCN.
- **Scale dependence:** none between 50 and 120 pc (same-pixel flat to 5%). The aperture contrast falls ~15% by 250 pc, because the N ratio falls with beam size (alpha = -0.23+-0.03; S -0.10+-0.03; Bridge -0.10+-0.04) while S is flat. At 6.5" the contrast is gone (1.07-1.10, with blending). So the ~2x contrast of the Q1 1"-scale data is not a dilution artefact of the coarser beam; it holds at 50 pc.
- **Consistency:** Xb4b and Xb4d agree to <=0.03 for N, S and Bridge at 250 pc. Xb4d and Xb4f agree to <=0.04 at 6.5". The E1 and E2 ratios differ between arrays (E2 0.52 vs 0.28, E1 0.14 vs 0.26) because their HCN is faint and Xb4b resolves out extended flux, so quote them as array-dependent +-0.15.
- **Velocity-resolved (relative to local HCO+ centroid):** at 55 pc the N ratio is flat, 0.95-1.09 in all 50-km/s bins from -125 to +125 km/s (core 0.98+-0.01, wings 1.14+-0.20 and 1.09+-0.11). S is flat at 0.36-0.43 over -125..+75, 0.67 in the -175 bin, and wings 0.45+-0.04 (blue) / 0.61+-0.27 (red). At 250 pc N shows a red rise (core 0.68, red 1.20+-0.12), which is absent at 55 pc, so it is beam mixing (consistent with p0010). The HCN wings are not detected in E1 and E2.
- **Lobes:** the S blue lobe and red ridge have a ratio >= 1 at 250 pc in Xb4d (1.43+-0.16 and 1.27+-0.17), against 0.39 at the S core. Xb4b gives 0.96+-0.19 and 1.36+-0.32, so the sign agrees but the amplitude is array-dependent. Neither is detected at >5 sigma in HCO+ at 50 pc, so there is no sub-250 pc lobe ratio.
- **Literature:** Harada+2018 give HCN/HCO+ 1-0 of ~1.0 (N) and ~0.6 (S), high at the S outflow. The J-independent part of the N/S contrast is therefore real, but the 3-2 contrast (2.1) is larger than the 1-0 contrast (1.5-1.7; team-a p0007-r01 says 1.49 at 1.7"). Nishimura+2024 (CON-quest) report global 3-2 ratios 0.4-2.3 (mean 1.1) at 0.3"; I found no NGC 3256 entry in the text I read. Our N 1.0 and S 0.4-0.5 lie inside their range.

## Uncertain / not done
- The HCN zero level (team-a: factor ~4 in some wings) matters only for faint wings. The line-free baseline test is on the lf channels only. I did not test a continuum-proportional baseline, or HCN absorption at S.
- HCO+ peak-based positions: the HCO+ and HCN peaks can differ at S. The same-pixel and peak/peak estimators agree within 0.05, so this is a small effect.
- E1, E2, W and the lobes are limited by HCN or HCO+ S/N; the ridge ratios use a single r=1.29" aperture and one window.
- Missing short spacings: Xb4b at 250 pc agrees with Xb4d for N, S and Bridge, but the lobes and E clumps are not recovered consistently. The conclusion on N/S needs no extra short-spacing correction at the 0.05 level.
- No physical conversion (n, T) of the ratio contrast was done in this round. It is an input to Q27/Q-excitation work: the ratio is 1.0 vs 0.45 at the same HCO+ T_b, so it is a column/excitation/abundance contrast in HCN, not in HCO+ (cf. p0030).
