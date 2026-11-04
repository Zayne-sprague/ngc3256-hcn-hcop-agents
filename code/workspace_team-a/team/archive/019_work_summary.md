# Q32 work summary (round 2, after reviewer "fix")

## What I did
- `scripts/q32_maps.py` (loader), `scripts/q32_models.py`, `scripts/q32_tk2.py`, `scripts/build_report_q32r2.py`
  -> `out/q32/models.json|models_full.log|tk2.json|tk2.log|maps.png`.
- Reviewer fix 1 (Q32b restated): rigorous limits = ff=1 values; ff/smooth-core values labelled conditional; core+pedestal/ff(core) numbers (44/29-30 K) withdrawn;
  "replaces" and "ff=0.5 not supported" dropped (ff=0.5 is not excluded: sub-0.18" clumping is unconstrained).
  One point I did NOT follow: T0/ff on top of the deconvolved T0 would double-correct (T0 is already deconvolved for a smooth Gaussian); the sub-clumping grid (ff_sub 0.7/0.5/0.3) is given instead.
- Fix 2: 5 profile families (Gauss+const, Gauss+free-size extended, single Gauss, exponential+const, Moffat+const) -> s, ff(0.26"), T0, Tk range. Confidence set to low.
- Fix 3: direct 2D forward fit (beam-convolved core + extended + const + 0-3 clumps, R=1.0/1.5") to the native map. (uv-plane fit on restored visibilities NOT done.)
- Fix 4: all fits repeated for peak channel, 3-ch mean, HCN-window mean; explicit flux ratio is the first sentence of Q32a.

## Results (peak channel, K RJ)
- Flux in beam 0.26"/1.0": 0.18 (S), 0.12 (N) [window 0.17/0.11]. Resolved at 0.18" regardless of model.
- S: s=0.35-0.71" (best 0.38-0.52), ff(0.26")=0.60-0.89 (0.68-0.80), T0=29-47 K; 2D fits s=0.38-0.57, ff 0.68-0.83.
- N: s=0.14-0.27", ff(0.26")=0.32-0.52; single Gaussian rejected; map is filamentary, N ff indicative only; 2D fits unstable.
- Tk: rigorous ff=1 native: S>=40, N>=21 K (0.26": 36, 18-19; window mean S>=30, N>=18). Smooth-core T0-based: S 34-53, N 18-33 K. Sub-clumping raises these (grid in report).
## Literature
Brunetti+2021 (NGC 3256 CO smooth to ~55 pc) cited as support for smoothness assumption only; not re-read this turn.
## Uncertain
Profile family (dominant), pedestal real vs missing flux (12m-only; not tested with 7m/SD cubes), sub-beam clumping, window-vs-peak, 13CO core size (0.85" beam), uv-plane fit not done, N core non-Gaussian.
