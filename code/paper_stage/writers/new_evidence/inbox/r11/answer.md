# r11: estimator spread at 1.5" (s192 uv-tapered HCN, HCO+ 3-2, N/S)

Scripts: `r11_apertures.py` (→ `r11_apertures.txt/.json`), `r11_xb4b_imgtaper.py` (→ `r11_xb4b_imgtaper.txt/.json`), `r11_fig.py`; figure `estimator_spread.png`.
Inputs: the s192 cubes `/workspace/out/q34/{HCN32,HCOp32}_{1p5,2p2}_{comb,c12}.pbcor.fits` (made by `scripts/s191_q34_casa_image.py`). Measurement conventions are those of `scripts/s192_q34_measure.py`: centre = HCO+ moment-0 peak within 0.4" of the nominal position, ±350 km/s window, baseline from line-free channels, RJ T_b. The same pixel and the same window are used for HCN and HCO+, for peak and for every aperture. The cubes have 0.3" cells, so apertures are evaluated on a 0.05" grid by cubic-spline interpolation.

## (i) N/S contrast (HCN/HCO+ N over S), 1-sigma from line-free channels

| beam | data | peak pixel | r=0.25" | r=0.5" | r=0.75" | r=beam |
|---|---|---|---|---|---|---|
| 1.5" | 12m+7m | 1.79±0.04 | 1.76±0.04 | 1.73±0.03 | 1.70±0.03 | 1.54±0.03 |
| 1.5" | 12m only | 1.90±0.04 | 1.87±0.04 | 1.83±0.04 | 1.78±0.04 | 1.59±0.03 |
| 2.2" | 12m+7m | 1.66±0.03 | 1.67±0.03 | 1.67±0.03 | 1.66±0.03 | 1.45±0.02 |
| 2.2" | 12m only | 1.69±0.03 | 1.70±0.03 | 1.69±0.03 | 1.67±0.03 | 1.46±0.03 |

Per-nucleus HCN/HCO+ ratios (N, S) for every cell are in `r11_apertures.txt`. For 1.5" 12m+7m: N 0.742, 0.740, 0.733, 0.719, 0.661 and S 0.414, 0.420, 0.423, 0.424, 0.429 for the same estimators in order.

My r=beam values reproduce s192 (1.532 and 1.578) to 0.006–0.008 (sub-pixel interpolation versus pixel mask). The s192 values in `out/q34/s192.txt` stand.

## Reading the spread
- The 1.91 (team-a peak pixel, 12m+ACA tapered images) vs 1.53 (s192, r=1.5") spread is mostly aperture overlap. Within r≤0.75" the contrast falls slowly (1.79 → 1.70 at 1.5", 12m+7m), which is dilution by beam smoothing and finite aperture. Going from 0.75" to 1.5" it drops a further 0.16 (12m+7m) or 0.19 (12m-only). The aperture centres are 1.45" apart, so a 1.5" radius aperture on S contains N: N's ratio falls (0.72 → 0.66) while S's rises (0.42 → 0.43).
- At 2.2" the contrast is flat at 1.66–1.70 for r≤0.75" and drops to 1.45 at r=2.2". The aperture overlap alone produces that drop, which is a clean demonstration because there is no dilution trend.
- The remaining 1.79 vs 1.91 at 1.5" is the array term: 12m-only is higher than 12m+7m by 0.08–0.11 at 1.5" (0.03 at 2.2"). Team-a's 1.91 matches my 12m-only peak (1.90) better than my 12m+7m one. I did not reproduce team-a's images, so that is a numerical match, not an identification. Any residual difference could come from weighting, uv range, or the restoring beam.
- Recommendation: quote the peak-pixel or r≤0.5" value as the "contrast at the beam" and the r=beam value only as a beam-averaged one. The estimator choice alone moves the 1.5" number by 0.25–0.30.

## 12-m-only 1.5" value 1.578
Confirmed: 1.578 in s192 is the 3-2 12m-only uv-tapered N/S contrast at 1.5" (circular beam; taper 1.1", r=1.5" aperture, ±350 km/s; `out/q34/s192.txt`, line "c12 w350 N/S"). From `s192.json`: HCN/HCO+ N = 0.660 (HCN 38.37, HCO+ 58.11 K km/s), S = 0.418 (HCN 16.75, HCO+ 40.04 K km/s), N/S = 1.578±0.027. The 12m+7m counterpart is N 0.658, S 0.430, N/S 1.532±0.026. Whether the 1-0 value at 2.5" is also 1.577 I have not checked (no 2.5" 1-0 measurement in my s183 output); the agreement would be a coincidence between different lines and beams, since the 3-2 number is fixed by the cubes above.

## Convolution vs taper term (Xb4b)
- Not done as a true uv-test. A Gaussian weighting applied in the image plane to a restored image is mathematically the same as convolving it with a Gaussian, so it cannot separate "convolution" from "taper". The difference arises in the visibilities (re-weighting, re-cleaning, and the residual/dirty-beam handling), and I have no Xb4b visibilities in the workspace (the Xb4b measurement sets were not produced; the q64 cubes at `out/q64` are Xb4d 12m/12m+7m only). A tclean uvtaper of Xb4b would need a CASA restore of Xb4b raw data, which I did not start.
- What I did instead: the Xb4b restored 0.24" cube, image-plane convolved to circular 1.05" and 1.5" (T_b conserved, same estimators). Results in `r11_xb4b_imgtaper.txt`.

| beam | peak | r=0.25" | r=0.5" | r=0.75" | r=beam |
|---|---|---|---|---|---|
| 1.05" | 2.13±0.05 | 2.11±0.05 | 2.08±0.05 | 2.03±0.05 | 1.98±0.05 |
| 1.5" | 2.03±0.05 | 2.00±0.05 | 1.98±0.04 | 1.95±0.04 | 1.83±0.04 |

- The convolution-vs-taper term at 1.5" (same data, so same array): Xb4d image convolved to 1.5" gives 1.809, versus 1.578 (12m only) and 1.532 (12m+7m) from the uv taper (`s192.txt`, "Xb4d image 1.809"). The term is 0.23 to 0.28 in the r=beam estimator, consistent with team-a's 0.2–0.3. The Xb4b convolved value (1.83) matches the Xb4d convolved value (1.81), so the array difference between Xb4b and Xb4d is small there. The same comparison at peak pixel: Xb4b convolved 2.03 vs uv-taper 1.79 (12m+7m) or 1.90 (12m only), a term of 0.13–0.24.
- At 1.05" I have only the convolved side. A matched 1.05" uv-taper exists for Xb4d 12m (`out/q64/*_c12_1p05.image.fits`) but I did not compare it, as the request needs Xb4b and the contrast definitions in q64 differ (Jy km/s annuli, r<1.0"). The Xb4d image-plane vs uv-taper at 1.05" is therefore not tested here.

## Caveats
- Uncertainties are statistical (line-free channels); the peak-pixel and aperture values within one beam use the same noisy pixels and are strongly correlated, so the trend is far more significant than the individual error bars suggest.
- The peak pixel is a position choice (HCO+ peak), not the position of each line's own peak.
- Some of the convolution-vs-taper term can be differences in cleaning (auto-multithresh) between the cubes; I did not separate that.
