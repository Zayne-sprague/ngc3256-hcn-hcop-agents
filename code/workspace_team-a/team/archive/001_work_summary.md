# Work summary: Q15 (red-shifted dense-gas ridge ~3.3" N of S; is the S outflow bipolar in dense gas?)

## What I did
Re-measured the ridge independently in Xb4d and Xb4b (both channel-wise convolved to a 1.05" circular beam, per-pixel linear baseline), in the Xb4d tapered cubes (baseline-free check) and in Xb47 4-3 (4.5" beam). Scripts (re-runnable, `python3 -I scripts/<name>.py`, outputs in out/q15/):
- scripts/q15a_cache.py: caches baseline-subtracted HCN/HCO+ 3-2 cubes (Xb4d, Xb4b).
- scripts/q15b_ridge.py (-> q15b.log): red-wing maps, blind maxima, sliding line-free-window null, mirror/core-scaled controls, 2D Gaussian centroid.
- scripts/q15c_pv.py (-> q15c.log): aperture spectra, per-velocity-bin centroids, p-v cut along the S-N axis.
- scripts/q15d_baseline.py (-> q15d.log): red-side baseline systematics (all Q3 free channels are blue of the line).
- scripts/q15e_flux.py (-> q15e.log, q15e_table.csv, q15e_fig.png): ridge/lobe/N fluxes, Xb4d vs Xb4b, random-aperture noise, baseline variants.
- scripts/q15f_tapered.py (-> q15f.log): tapered cubes (30 km/s channels).
- scripts/q15g_43.py (-> q15g.log): HCN/HCO+ 4-3 vs 3-2 at 4.5".
- scripts/q15h_kin.py (-> q15h.log): core centroids, wing profile along S-N, dynamical times.
- scripts/build_report_q15.py: adds Q15a-h results, recipe step, caveat, interpretation to /output/report.json (valid against schema; also fixed the invalid 'medium-low' confidence enums to 'low').

Windows: red +230..+470 km/s (2960-3200), blue mirror -470..-230 (2260-2500), v0=2730; r=1" aperture mean; 1"=194 pc (D=40 Mpc assumed).

## Results
- **Real?** HCN: ridge 0.126+-0.022 Jy/beam km/s (5.9 sigma pointwise, noise = 400 random apertures in the same window; baseline systematic +-0.05, i.e. +-40%); tapered, baseline-free cube 0.093 vs control scatter ~0.01; no blue wing at the ridge (-0.002); wing/core 17% at the ridge vs 2% at N, and the ridge is 1.7" from N (beam leakage negligible). Recovered in Xb4b: 0.106+-0.023 (ratio 0.84+-0.22 to Xb4d). HCO+: only 0.038+-0.018 (2.1 sigma; tapered 0.051, ~4 sigma; team-b's 5.7 sigma uses a +180..+320 window that includes the core wing). Detected in HCN (medium), tentative in HCO+.
- **Position/velocity:** red-wing centroid 3.2-3.5" (window-dependent to 4.0") N of S = 1.7" S of N, constant with velocity; 0.6-0.7" from Michiyama+2018 OF-Red. Red HCN wing is a continuous tongue from S+1.75" to N (peak 0.19 at S+3.25", N itself 0.07-0.10, S..S+1" <0.02). Wing <v> +290 (HCN) / +221 (HCO+) vs lobe -262/-265 km/s. Core centroids: S 2753, N 2782, Bridge 2786, ridge 2793, lobe 2722-2736.
- **Ratios/mass:** HCN/HCO+ ridge 3.3+-1.8, lobe 2.0+-0.7 (cannot separate). Ridge/lobe HCN flux 1.23+-0.26 (Xb4d; 3.5+-2.2 in Xb4b because the lobe is resolved out: 0.30+-0.19 recovery). If the ridge is the counter-lobe the dense outflow mass is 1.7-2.2x the blue-lobe-only value.
- **4-3:** no usable constraint (blended with N/S at 4.5", HCN 4-3 baseline-limited: -0.33+-0.12 or +0.23+-2.0).
- **Origin: undecided.** For counter-lobe: OF-Red coincidence, symmetric |v| and flux, both wings detached ~2" from S. For Bridge/N: ridge core velocity matches N/Bridge (+40 km/s vs S), continuous tongue to N, N's own red-only wing. Against Bridge/N: wing +260 km/s vs 29 km/s N-S core difference. Bipolarity of the dense S outflow is permitted, not established; do not add the ridge to the S mass budget yet.

## Uncertain / remaining
Red-side baseline extrapolation (±40% HCN flux; decides HCO+); pointwise (a-priori) significance only, non-Gaussian red-window noise; no CO or low-J wings at <=1" in the archive for kinematic discrimination; 4-3 unusable; fluxes are 1.05" aperture values from uncleaned channels; dynamical times unprojected. Q13 stays ranked above Q15.
