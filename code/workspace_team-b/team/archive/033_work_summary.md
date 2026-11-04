# Q20 work summary (worker, team-b)

**Question:** is the red ridge (+180..320 km/s, 3.3" N of S) seen in independent lines with the same profile, and do ridge-vs-lobe line ratios discriminate S counter-lobe from Bridge streaming?

## Answer
- **Weak independent confirmation in HCN/HCO+ 4-3 (Band 7), ~3.4 sigma stacked.** Ridge excess relative to the 3-2 excess (s): HCN 4-3 1.64±0.73 (2.3σ), HCO+ 4-3 3.6±1.2 (2.9σ), stack 2.14±0.63. Lobe: 1.25±0.41, 1.16±0.49, stack 1.21±0.32 (3.8σ). The same fit recovers s≈1.3±0.4 for the 3-2 ridge. Mirror controls are null. Robust to the template taper (1.5" vs 2.5").
- **Band 7 cannot do the 1" test** (beam 4.1x2.7", N 1.9" and S 3.3" from the ridge). So all comparisons are at one 4.5x3.0" beam, using a template fit that assumes the 3-2 spatial shape, with 35-45% errors.
- **Profile:** the 4-3 ridge also starts at ~200 km/s (fraction in 160-200 ≈0 for HCO+, 0.18±0.17 for HCN); the HCO+ lobe falls from the lowest |v| (centroid 179 +30/-55 vs ridge 282 +31/-22: ~2σ). In HCN the shape difference is not distinguishable. Mild support only for the p0017 shape difference.
- **Ratios do not discriminate:** r43/32 (Tb) HCN ridge 0.72 (+0.51/-0.34) vs lobe 0.80 (+0.38/-0.29); HCO+ ridge 1.31 (+1.06/-0.54) vs lobe 0.50 (+0.24/-0.22). HCN/HCO+ of the excess is ~1 in both places (3-2 ridge 1.26 +0.97/-0.52, lobe 0.96 +0.31/-0.25), well above the core ratio (0.4-0.5), so it is lobe-like but the Bridge-like alternative is excluded at only ~1.4σ.
- **CS 7-6, SO, HC3N, H13CN: all non-detections, far above any expected level.** CS 7-6 ridge s<2.96 (9-26x the core ratio). SO 6_5-5_4 at 1" r=1" apertures: ridge <0.082 Jy km/s (<0.32 of the HCN 3-2 excess); lobe <0.089. HC3N 28-27: ridge <0.131, lobe <0.117. H13CN 4-3: s<7.6, uninformative (CO 3-2 tail on the baseline, spw ends at -225 km/s so no blue mirror; H13CN 3-2 is not covered). Stacking cannot help.
- **Origin stays undecided** (S counter-lobe ~55-65%, unchanged).

## Scripts (in /workspace/scripts, all re-runnable with `python3 -I`)
- s110a_q20_spectra.py: native peak-pixel spectra and per-channel noise (out/q20/s110a*.csv/png/txt). Noise mJy/beam: HCN43 5.2, HCO+43 6.0, CS 4.6, H13CN 4.3, SO 1.3.
- s110b_q20_maps.py: matched-beam window and 40 km/s bin maps (out/q20/s110b_maps.npz).
- s110c_q20_templatefit.py [taper]: template fit, shifted-template null, bins (out/q20/s110c*.txt/json).
- s110d_q20_so_hc3n_1arcsec.py: SO and HC3N at 1.05" with the s106a-style aperture null (out/q20/s110d*).
- s110e_q20_summary.py: ratios, stacks, ULs, shape statistics (out/q20/s110e.txt).
- s111_update_report_q20.py: adds 3 results and 1 caveat to /output/report.json (now 177 results).

## Uncertain / not done
- The template assumption (4-3 excess shaped like the 3-2 one); the HCN 4-3 core-wing coefficient is negative (nuclear continuum or absorption) and is not physically modelled.
- The HCO+ 4-3 ridge excess is above the thermal limit (1σ), pointing to unmodelled systematics.
- The 3-2 and 4-3 fit errors are treated as independent in the ratio errors (overstated if anything).
- The line-independent test is limited to HCN/HCO+ 4-3. The weak-line limits (CS, SO, HC3N, H13CN) are too shallow to test anything. A JWST/CO velocity-sign map beyond 250 pc is still needed for the origin.
