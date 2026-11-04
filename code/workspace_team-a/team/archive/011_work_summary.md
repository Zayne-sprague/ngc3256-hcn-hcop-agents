# Work summary: Q7 (N-core HCN/HCO+ 3-2 radial profile, disk 2-1 ratio; check of team-b p0005/p0007)

## What I did
Independent reproduction with our own code (no team-b scripts). Central ±7" crops of the Xb4b and Xb4d pbcor cubes (HCN spw27, HCO+ spw25) were cached, baseline-subtracted per pixel (linear fit, 520-900 km/s blue and >520 km/s red, outside HCN v2=1f at +401 km/s and HOC+ at -1000 km/s), and integrated to moment-0 maps. Both lines were convolved to a common circular beam: Xb4b to 0.30", and to 1.05" for a matched-array test against Xb4d. Annuli were centred on the HCO+ peak at that beam near each nucleus. N is at (0.00,+2.63)", S at (-0.15,-2.64)" from 156.96346,-43.90461. Noise comes from null windows (200 km/s, baseline refit excluding the window, scaled by sqrt(width)) with annulus sums at 300 random centres. I added an injection with a flat intrinsic ratio, a compact-excess fit, native-beam peak spectra, and an own-code Xb55 disk annulus calculation.

Scripts (re-runnable with `python3 -I`, in this order):
- `scripts/q7_cache.py`, then `q7_lib.py` (library)
- `q7_profiles.py` -> `out/q7/q7_profiles.csv`, `q7_profiles.log`
- `q7_inject.py` -> `out/q7/q7_inject.log`
- `q7_core.py` -> `out/q7/q7_core.log`
- `q7_disk.py` -> `out/q7/q7_disk.log`

## Results
Convention: ratio of integrated Jy km/s (same beam for both lines), ±350 km/s about 2780 km/s radio LSRK, 40 Mpc (1" = 194 pc).

1. **The N inward rise is reproduced.** Xb4b at 0.30", annuli 0-0.15 / 0.15-0.25 / 0.25-0.5 / 0.5-1 / 1-1.5 / 1.5-2 / 2-3":
   - HCN/HCO+ = 0.98 / 0.95 / 0.87 / 0.64 / 0.45 / 0.42 / 0.30 (±0.02 inside 1", ±0.05 floor outside).
   - Cumulative: 0.98 (r<0.15"), 0.90 (<0.5"), 0.74 (<1"), 0.57 (<3").
   - These agree with team-b p0005 (0.97, 0.87, 0.65, 0.41, 0.39).
   - The Xb4b 2-3" value (0.30) is lower than Xb4d's (0.40) because Xb4b loses extended flux. Use Xb4d for r>1.5".

2. **S is flat.** Xb4b at 0.30": 0.42, 0.41, 0.40, 0.36, 0.33, 0.37 out to 2". The 2-3" annulus (0.46) is blended with Bridge and Soutflow.
   - Xb4d at 1.05" (0.5" bins): 0.39, 0.36, 0.36, 0.39, 0.385 to 3".
   - The N inner/outer contrast is a factor of 2-2.5. The S contrast is none.

3. **The rise is not an array artefact.** Xb4b smoothed to 1.05" against Xb4d native, for r<0.5 / 0.5-1 / 1-1.5 / 1.5-2 / 2-3":
   - Xb4b: 0.786, 0.698, 0.555, 0.434, 0.329.
   - Xb4d: 0.805, 0.726, 0.599, 0.474, 0.399.
   - The two agree to 2-4% inside 1". The offset grows to 8% at 1-2" and 18% at 2-3", where Xb4b loses extended flux.

4. **Injection.** A flat ratio R0 = 0.4, 0.6 or 1.0 is returned flat. Offsets are ≤0.02 for r<0.5", and ±0.03 and ±0.01 at S and N for 0.5-1". The real null noise at N 1.5-2" adds +0.09 (about 2.5σ), so the outer-annulus error floor is ±0.05. The pipeline does not create the gradient.

5. **Velocity windows.** At r<0.25" the ratio is 0.95-1.00 for ±100, ±165 and ±350 km/s. The wing windows are noise-limited (N r<0.5" at 1.05": blue 0.78±0.15, core 0.77, red 1.13±0.21).
   - Q3's N blue 1.13 vs core 0.74 was a 2" aperture and a different window. It is not reproduced in the nuclear aperture, but it is not excluded either.

6. **Size of the HCN excess.** The excess is not point-like. Fitting HCN = k·HCO+ + A·Gaussian at 0.30" over r<2":
   - The best FWHM is about 0.6" (~115 pc).
   - Point-source rss is 1.56x the minimum. s=0.3" gives 1.30x; s=1.0" gives 1.12x.
   - The HCN excess over 0.38·HCO+ is 3.4 / 6.2 / 6.7 Jy km/s within r<0.5 / 1 / 1.5". Half of it lies within about 0.65".
   - So the high-ratio region is a gradient on ~0.5-1" (100-200 pc) scale. The ratio is ≈1 only at the peak pixel, so "compact (<0.25")" overstates it.

7. **Native-beam peak spectra (Xb4b, 0.25").**
   - N: HCN 8.9 mJy/beam (2.69 K, v=2823 km/s); HCO+ 9.4 mJy/beam (2.87 K, v=2841 km/s). Peak T ratio 0.94, matching team-b.
   - S: HCN 4.4 mJy/beam (1.33 K, v=2691 km/s); HCO+ 8.5 mJy/beam (2.58 K, v=2723 km/s). Peak T ratio 0.52, while the integrated ratio is 0.42 (the S HCN line is narrower or shifted).
   - rms per 4.4 km/s channel: 0.31-0.44 mJy/beam.

8. **Disk 2-1 (Xb55, 10" beam, far-field zero level).** Annuli about the N-S midpoint: 0.580 (0-4") / 0.555 (4-8") / 0.492 (8-12"). This agrees with p0005 (0.575/0.554/0.489).
   - The zero-level scatter is ±0.004. The 12% drop is real in these data.
   - It was measured in 2-1 only (7-m only). It is consistent with the radial decline seen in N, since the interior ratio is weighted towards the nuclear gas.

## Verdict on Q7
Confirmed (partly corrected):
- HCN/HCO+ 3-2 rises inward in N from ~0.35-0.4 at 2-3" to ~1.0 at the peak.
- S is flat at 0.33-0.43.
- The disk 2-1 is ~12% lower at 8-12".
- Correction to p0005: the excess has a ~0.6" FWHM (~115 pc), not <0.25", and it is a ratio of fluxes, not an abundance ratio. This feeds Q4: opacity vs abundance is still undecided.

## Still uncertain
- Annulus errors are null-window based. Outer annuli carry a ±0.05 floor and Xb4b loses extended flux there.
- The Gaussian-excess fit has no formal errors and degenerate k.
- Per-channel beams were taken from the header.
- I did not do velocity-resolved profiles; that is Q9.
- The disk ratio has an unquantified resolved-out-flux systematic.
