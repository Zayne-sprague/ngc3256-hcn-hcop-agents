# Worker summary - Q0, round 2 (team-b) - answers to the reviewer's 6 fixes

Re-run order (all `python3 -I`, from /workspace): s11_wing_tests -> s11b_wing_geometry -> s12_line_table -> s12b_spatial_scatter -> s13_region_velocities -> s14_continuum -> s15_crossous -> s16_b5_overlap -> s17_ch3oh_ladder -> build_report2 (reads out/report_base_s7.json = round-1 report, writes /output/report.json, 88 results). Round-1 scripts (s0-s10, build_report.py) unchanged.

## 1. Is the S/Sblue blue wing real?  YES (existence); fluxes are lower limits
Script s11_wing_tests.py (+ s11b_wing_geometry.py, figure out/s11_wing_maps.png, numbers out/s11_wing_tests.txt/.json).
- Empirical noise: beam-aperture blue-wing (-450..-150 km/s) vs line-free pseudo-wing maps (-900..-600, +450..+750), 14k-28k pixels. Xb4d Sblue peak: HCO+ 0.245 Jy/beam km/s = 13.8 sigma (0 of 13869 pseudo pixels exceed; max 0.083), HCN 0.346 = 19 sigma; Xb4b 8.1 / 12.6 sigma. S nucleus wing 8.8 sigma (HCO+), 4.0 sigma (HCN).
- Sidelobe model (Sblue = a*S + b*N fitted on core channels, same spectral shape as the dirty-beam response): predicts 0.009 (HCO+), 0.022 (HCN) vs observed 0.182+-0.022, 0.139+-0.021 in Xb4d; wing/core 0.22 (Sblue) vs 0.035-0.055 (S) vs <=0.01 (N). Even a 10% sidelobe of S gives 0.014.
- Independent array: Xb4b convolved to the Xb4d beam gives 0.154 / 0.215 at the Sblue peak (Xb4d 0.245/0.346 at its peak, 0.18/0.14 nominal) and 0.114/0.087 at S (Xb4d 0.156/0.074): same structure with a different dirty beam.
- The "10x brighter at 0.94" than 0.25"" argument is wrong as an artefact test: uniform emission is 12.1x brighter in Jy/beam in the larger beam (area ratio). The 0.25" red-wing negative is only -1.6 sigma. Caveat: Xb4b HCN alone cannot discriminate sidelobes (0.016+-0.010 vs 0.012).
- Geometry: blob is an elongated (~2.3x0.7", long axis PA~15 deg), centroid 1.8-2.0" S of the S nucleus (PA 183-189 deg), flux ~0.6-0.75 Jy km/s: attached to S, "Sblue" is its southern end. Consistent qualitatively with the literature (outflow from the S nucleus, Michiyama+2018 abstract); the published outflow PA was not verified.

## 2. Unified per-line / per-region table
s12_line_table.py -> **out/s12_line_table.csv** (494 rows = 21 lines x regions x OUS x {pix, ap}). One convention: template = region's own HCO+ profile of the same OUS (window = >10% of peak, region-specific centre); I_T = matched-filter amplitude x template integral (Jy/beam km/s), W_sum (window box sum), peak (3-ch smoothed), per-channel rms (line-free, lag-1 inflated), sig_stat, sig_emp (template slid over line-free channels = baseline structure), sig_base, sig_cal (5% B6 / 10% B5,B7 + pb term), 3-sigma UL, K km/s, pb, blend response. Lines: HCN/HCO+ 2-1/3-2/4-3, HOC+, H41b, CH3OH blend, CH3CN, CS 7-6, H2CS, HC15N, SO 8_8-7_7, SO 5_6-4_5, NS e/f, SiS, HC3N, HCN v2=1f (3-2, 4-3). Limits: Xb4d typically 0.04-0.07, Xb4b 0.02-0.04 Jy/beam km/s (3 sigma). Per-line 1G/2G (v, FWHM, chi2r) for HCN/HCO+: out/s7_line_fits.csv; peaks/Tb there too.

## 3. Tentative lines: budget and confidence (spatial scatter test: s12b_spatial_scatter.py -> out/s12b_spatial_scatter.json)
| line | formal | vs spatial field scatter | confidence |
|---|---|---|---|
| CH3OH J=1-3 E blend (Xb55) | N 0.93, S 1.33 Jy/beam km/s, 17/13 sigma | 14 / 18.5 sigma | high; J=4,5 absent => Trot<~11-12 K (LTE thin), cold/subthermal; spw20/22 overlap is line-free and consistent, but contains no CH3OH |
| HOC+ 2-1 (Xb55) | N 0.96+-0.09 (11 sigma), S 0.73+-0.12 | 9.5 / 6.8 sigma | medium; single line, frequency assignment, no Eu<300 K competitor; N+S blended in the beam; HCO+/HOC+ = 37+-3 (N), 51+-8 (S), same order as NGC 1068 (30-100) and Harada+2018 NGC 3256 (~0.04 HOC+/HCO+) - literature from search snippets only; Sato+2022 not found |
| SO 5_6-4_5 (Xb4f S) | 0.28+-0.06, 5.0 sigma | 2.8 sigma | low/unconfirmed |
| 253.9 GHz feature (Xb4f S) | 0.53+-0.05, 9.8 sigma | 4.2 sigma (field max 0.86) | low/unconfirmed (broad blend/baseline); the "e-partner absent => NS disfavoured" claim is WITHDRAWN (limit 0.58 > expected 0.5) |
| HC3N 28-27 (N) | Xb4f 0.24+-0.06 (3.8), Xb4d 0.04+-0.014 (2.9) | 2.0 / 5.5 sigma | low |
Non-detection at 1" and 0.25" is NON-diagnostic: a smooth Xb4f source would dilute to 0.010 (SO), 0.018 (253.9) Jy/beam km/s in Xb4d, below the 3-sigma limits (0.05, 0.046). So these are not shown to be artefacts, just not confirmed.

## 4. Continuum (s14_continuum.py -> out/s14_source_catalogue.csv, s14_common_beam.csv, s14_alpha.csv, s14_continuum.txt)
- Blind >5 sigma finder, local pb-aware rms: Xb4b 38 peaks (16 with reliable 2D Gaussian fits), Xb4d 8, Xb47 3, Xb4f 1, Xb55 1; 0 peaks < -5 sigma in any image. Positions with Condon errors (0.001-0.02") + 0.1-beam astrometric floor: N (156.96345,-43.90388), S (156.96341,-43.90533), sep 5.22+-0.04"; SE (156.96684,-43.90635) 25 sigma; E1, E2 as before. S is 0.22x0.20" (marginally resolved) and its 0.24" peak holds only ~43% of its Xb4d Gaussian flux. Remaining Xb4b peaks are clumps in the N complex.
- Common beams (5.84x3.52" and 8.95x5.98"): N/S peak ratio 1.4-1.6 at every matched beam => "S compact / N extended" in round 1 was beam dilution + resolved-out flux. r<8" B6 totals at the common beam: Xb4b 24, Xb4d 35, Xb4f 51 mJy (48% / 68% recovered), so indices using Xb4b/Xb4d anchors are biased (alpha 3-4.7 or negative). Defensible: alpha(B5-B6)=1.26+-0.27, alpha(B6-B7)=1.92+-0.38, alpha(B5-B7)=1.54+-0.20 (+~0.5 uv-coverage systematic): flatter than pure dust. N and S are merged at B5 (per-nucleus B5 index meaningless).

## 5. Per-region systemic velocities (s13_region_velocities.py -> out/s13_region_velocities.csv, s13_region_vsys.csv)
Half-maximum midpoints (HCO+/HCN 3-2, Xb4b/Xb4d), km/s radio-LSRK: N 2792+-19, S 2721+-21, Bridge 2792+-7, E1 2832+-5, E2 2860+-10, W 2776+-5, Sblue 2734+-9. N-S offset +71+-28; N-S midpoint 2756+-15 (~ the nominal 2750). Round-1 "2780+-5" is the blended low-resolution centre. 1G centroid/FWHM flagged INVALID wherever chi2r>2 (N, S at all beams) in every round-1 line row of the report; templates/windows of all weak lines are per region.

## 6. Flux scale / pb / cross-OUS (s15_crossous.py -> out/s15_crossous.txt)
sig_cal column + text: 5% B6, 10% B5/B7; +5% (pb<0.9) or +10% (pb<0.75) pb term from a 12 m Gaussian model (pb E1 0.86/0.71/0.54, E2 0.82 B6, W 0.72 B6; B7 E1/E2/W 0.54-0.70; assumed array type, not measured). Xb4b convolved to the Xb4d beam vs Xb4d: continuum peak ratios 0.92-1.02, HCO+/HCN core moment-0 0.91-1.01 for bright regions (rms ~6%), 25-35% for faint ones. All cross-OUS comparisons are image-plane only and need matched beams AND uv coverage.

## Still uncertain
- Wing fluxes are uncleaned (lower limits; systematic <=0.01-0.02 Jy/beam km/s); outflow PA vs literature not checked against the papers.
- pb accuracy/array type for E1/E2/W (esp. Band 7); flux scale values assumed nominal.
- Tentative lines (SO 5_6-4_5, 253.9 GHz, HC3N) need other data; HOC+ only one line; CS/H2CS/29SiO and CH3OH J=1-3 stay blended; H13CN not separable from CO wing.
- sig_emp has few independent windows (n_emp column) for Xb47 spw20/22 lines; sig_tot falls back to statistical there.
- Some round-1 entries kept as is (coverage, mfs/cont line contamination, mask check) - not re-derived.
