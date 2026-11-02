# Q25 work summary (round 2)
Question: does HCN/HCO+ 1-0 at a 10" beam show the same 8-12" disk decline (~11%, Q11/p0032) as 2-1?

## Answer
Yes within the systematics, but the decline is an upper bound and cannot be pinned below ~5%.
8-12"/4-8" ratio factor: Michiyama 1-0 **0.90-0.94** (decline 6-10%), Xb55 2-1 **0.89-0.91** (decline 9-11%). R21/R10 (8-12 over 4-8) is 0.95-1.01 across all scenarios (flat to +-5%). Harada 1-0 (0.73-0.78) is recovery-limited and not trusted at 8-12".

## What was done (this round)
- `scripts/q25_systematics.py` -> `out/q25/q25_systematics.log` (earlier: q25_ratio10/growth/uvtest/radex*.py, `build_report_q25.py`).
- Bowl: mean of the 10" moment-0 in r=22-35" is -0.59 (HCN), -0.73 (HCO+) Jy/beam km/s for Michiyama; -0.54/-1.0 Harada (std 1.8-2.5, noisy); Xb55 -0.44/-0.81 but that ring lies outside its primary-beam half-power (poorly constrained). Line-free windows have |mean| <= 0.12, rms 0.05: the bowl is present only at line velocities (annulus spectra: -2 to -3 mJy/beam at line velocities vs 14-16 / 25-28 mJy peaks in the 8-12" ring), as expected for missing short spacings. Add-back scenarios (0, 0.5, 1 x ring offset): Michiyama factor 0.902/0.921/0.937, Harada 0.733/0.759/0.779, Xb55 0.889/0.892/0.895 (0.914 with the reviewer cumulative-flux 16-30" bound, `team/checks/q25_bowl_estimate.py`). Maximal add-back is an upper bound on the correction only if there is no real emission at 22-35".
- Window test (+-250/350/500 km/s): Michiyama 0.902/0.902/0.901, Xb55 0.887/0.889/0.898, Harada 0.757/0.733/0.663 (Harada unstable). 8-12" HCN and HCO+ spectra have the same width and no wings beyond +-350 km/s at the noise.
- uv/MRS: only image products were used; no visibilities restored, so MRS was **not measured** (headers carry no config history). Both 1-0 sets are 12m Cycle-3; min baseline ~15 m gives MRS ~10-30" depending on configuration (estimate). The choice of Michiyama rests on the image-plane Fourier test (Harada/Michiyama amplitude 0.6-0.7 at 1-4 kl, cumulative flux ratio 0.75 at 12" and 0.5 at 20") and the high-pass stability of the factor (Michiyama 0.90-0.92, Xb55 0.89-0.86 for 0.5-3 kl cuts), not the flux level; neither set is complete (cumulative flux falls beyond 16" in both). A common-uv re-imaging of all three (the 2015 visibilities and Xb55) would settle it; not done.
- Report updated: Q25a/b/c softened, Q25d added with the scenario table, Q11 entries carry an upper-bound note for the 2-1 11% (`scripts/build_report_q25_fix.py`, idempotent after `build_report_q25.py`).
- RADEX not re-run: a 5% R21 change is ~0.02 dex in log X vs a 0.25 dex grid step, so log X = +0.50 constant (delta 0.00 +-0.25) and an allowed 0-0.5 dex lower outer density stand; no abundance gradient or excitation change is required by the 1-0/2-1 data.

## Literature
Harada+2018 (ApJ 855,49) per-nucleus 1-0 ratios only (~1.0 N, ~0.6 S at 1.2"); Michiyama+2018 (ApJ 868,95) outflow paper; no published 8-12" disk ring ratio. I could not extract their uv/config details from the lit tool text.

## Uncertain
Bowl amplitude in 2-1 (ring outside PB); true emission at 22-35"; MRS of all three sets; Harada's 8-12" flux; the RADEX fit was not repeated per scenario; the 4-point per-ring statistics are small.
