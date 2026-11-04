# Q22 work summary: negative ramp at the Xb55 spw22 edge

## Answer
The ramp is a **continuum-proportional, edge-localised residual of the continuum subtraction** (edge amplitude response not followed by the low-order continuum model), not CH3OH emission and not a sky feature.
- Xb55 spw22 low-frequency edge (163.244 GHz): first 10 channels (1.953 MHz) -3.1/-3.5/-5.0 mJy/beam at N/S/Bridge (rms 2.1/channel), relaxing roughly linearly to ~0 at ch ~140-170 (CH3OH blend centre ch ~153, line spans ch ~120-190).
- Ramp / continuum (line-free cont image, pbcor): ch0-30 -0.21+-0.02 (through origin, chi2 23/6; E1 is an outlier, ~0), ch0-100 -0.157+-0.011 (chi2 16/6); zero at the continuum-free E2/Sblue.
- Controls (same pixels, same continuum): spw22 high edge 0.0+-0.03 x cont; spw16 and spw20, both edges, 0+-0.04 x cont. spw18 low edge shows the same ramp (-0.10 to -0.15 x cont, ch 0-80). So it is spw-specific (spw18 and spw22 low edges), not position-specific.
- The spw20/spw22 sky-frequency overlap (164.54-165.11 GHz) agrees (diff +0.14+-0.44, +0.13+-0.40, +0.11+-0.51 mJy at N/S/Bridge). The spw16/spw18 overlap is only 53 channels (one bin; spw18 lower by 0.6-1.6 mJy, not significant).
- Cannot test with uv data (image products only), so "edge bandpass roll-off x continuum, not followed by the continuum fit" is an inference from scaling + localisation.

## Consequence for the CH3OH blend (Q9 baseline systematic)
Old s120b grid used +-400 km/s windows, which contain the ramp (ch 41-110), hence the factor 2-4. With an empirical ramp model and +-200 km/s window (mean of R ramp-extrapolation, L line, Q quadratic; J=4,5 masked in the anchors; sys = half-range over methods and +-150..300 km/s windows), beam-aperture integrals in Jy/beam km/s:
| region | adopted | stat | sys (half-range) | range | zero-baseline bound |
|---|---|---|---|---|---|
| N | 1.13 | 0.08 | 0.42 (37%) | 0.71-1.55 | 0.77 |
| S | 1.61 | 0.07 | 0.38 (24%) | 1.23-1.99 | 1.15 |
| Bridge | 1.69 | 0.07 | 0.33 (20%) | 1.37-2.04 | 1.19 |
| W | 0.72 | 0.07 | 0.10 (14%) | 0.69-0.89 | 0.43 |
Range factors 2.2/1.6/1.5/1.3 (were 3.5/2.1/2.6/2.0). Adopted values within 1-8% of the s120b medians, so the Q9 central values stand but the uncertainty is now +-14-37%, not factor 2-4.
Trot (thin LTE, 254 GHz line at the Xb55 beam, MC over the baseline range): S R_K 0.23 (0.19-0.29), Trot 13 (10-20) K; Bridge R_K 0.16, Trot 30 K (16%: 17 K, upper end unbounded); N unconstrained. Slightly warmer than Q9's 10 (+7/-5) K at S. Blend/HCO+ ratios scale with the blend (S x1.08, Bridge x1.01, W x0.98, N x0.95): contrast conclusions unchanged.

## Scripts (python3 -I, in order) and outputs
- scripts/s130_q22_edge_ramp.py -> out/q22/s130.json, s130_edges.png (10-ch binned edge profiles, all regions/spws)
- scripts/s131_q22_edges_allous.py -> out/q22/s131.txt/json (edge levels vs continuum for all OUS/spws; Xb4f spw18 and Xb47 spw16 low edges also show negative edges of -0.06 to -0.27 x cont; reference windows contaminated by lines in several spws, qualitative only)
- scripts/s132_q22_ramp_fit.py -> out/q22/s132.txt/json (overlaps, exponential fits (poorly constrained, ramp is ~linear), continuum values)
- scripts/s133_q22_template.py -> out/q22/s133.txt/json (stacked ramp/continuum, slope vs continuum)
- scripts/s134_q22_blend_empirical.py -> out/q22/s134.txt/json (blend integrals, methods x windows)
- scripts/s135_q22_adopted_trot.py -> out/q22/s135.txt/json (adopted values, Trot MC)
- scripts/s136_update_report_q22.py -> /output/report.json (2 results, 1 recipe, 1 caveat added; backup workspace/report_backup_q22.json)

## Uncertain
- Shape of the ramp under the line (ch 120-190) is extrapolated; this is now the dominant systematic (R vs Q baselines differ by up to 2x at N).
- The ramp is only 2 spws x 1 edge; "edge roll-off" mechanism is inferred, not shown.
- Peak-pixel spectra are used for the ramp fits; aperture spectra for the integrals; continuum at region pixels from the 9" cont image (N, S, Bridge are blended within the beam, so their continuum values are not independent; the fit leans on E1, W, E2, Sblue for the zero point).
- Xb55 spw16/18 edge profiles are contaminated by HCN/HCO+ 2-1 wings beyond ~100-150 channels.
