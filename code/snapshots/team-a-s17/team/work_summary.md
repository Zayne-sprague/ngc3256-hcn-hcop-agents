# Work summary: Q23 (CO(2-1) T_b at the S and N cores), round 2 after reviewer fixes

## What was done
All four reviewer fixes were addressed with the archived CO(2-1) cubes. Core positions are the Q5 HCN/HCO+ centres, N (-0.06,+2.65") and S (-0.15,-2.59"). The cubes are 2015.1.00714.S 12m 0.18" (v<2844 km/s), 2015.1.00902.S 12m 0.82x0.68", the 2015.1.00714.S 7m mosaic and the SD cube. Beams were matched with radio_beam elliptical deconvolution kernels, not the circular geometric-mean shortcut.

Scripts and outputs (all re-runnable with `python3 -I`):
- `scripts/q23b_vmatch_pos.py` -> `out/q23b/vmatch_pos.json`: velocity-matched T_b, offset grid, Monte Carlo of the core centre.
- `scripts/q23b_missing_flux.py` -> `out/q23b/missing_flux.json`.
- `scripts/q23b_corrected.py` -> `out/q23b/corrected.json`. I fixed a crash: a circular 0.82" target cannot be deconvolved from the 0.82x0.68" beam, so it now uses only 1.0" and 1.7".
- `scripts/q23c_tk_ff.py` -> `out/q23b/tk_ff.json`: Tk_min versus ff and tau, plus the beam-dilution fit.
- `scripts/build_report_q23b.py`: writes Q23d-h to /output/report.json. The report validates (193 results).
- Round 1 (unchanged): `scripts/q23_co_hires.py`, `out/q23/co_hires.json`, `team/checks/q23_check.py`.

## Results (RJ T_b unless stated)
**1. Velocity-matched T_b.**
- S: the HCN line has v0 = 2711 km/s and FWHM = 117 km/s, and the CO peak is at 2712 km/s, so CO-bright and dense gas share the velocity.
- S mean T_b over the HCN FWHM is 20.1/12.7/7.8/4.8 K at 0.26/0.57/1.0/1.7". The peak channel is 31.2/19.2/11.8/7.0 K, so the mean is 0.64-0.69 of the peak. HCO+ gives the same values.
- N: the HCN line has v0 = 2797 km/s and FWHM = 166 km/s, so the 0.18" cube covers only 75% of the FWHM. The FWHM mean is 11.6/8.1/7.1/5.7 K against a peak of 13.3/8.5/7.5/6.0 K. With the 0.82" cube plus 7m+SD, the full-FWHM mean at 1.0/1.7" is 7.9/6.5 K.
- N CO peak channel (2763-2790 km/s) is offset from the HCN centroid by 7-35 km/s. The hottest CO within 0.3" of N (21.2 K RJ, 26.3 K Planck) is at (+0.23,+0.17") and v = 2727 km/s, so it is not at the HCN core.

**2. Position budget.** Offset grid at 0.26", min-max over 8 azimuths:

| Offset | S (K) | N (K) |
|---|---|---|
| 0" | 31.2 | 13.3 |
| 0.03" | 29.7-32.3 | 12.8-13.8 |
| 0.06" | 27.9-32.9 | 12.7-14.5 |
| 0.10" | 25.0-33.0 | 12.7-16.1 |
| 0.15" | 20.9-32.2 | 12.0-17.7 |

- Monte Carlo of the centre, median [16,84], at 0.26":
  - sigma_pos = 0.05": S 30.3 [27.6,32.2], N 13.8 [13.2,14.2].
  - sigma_pos = 0.10": S 28.1 [21.5,31.9], N 13.8 [12.6,15.4].
- At 0.57" with sigma_pos = 0.05": S 19.1 [17.7,20.0], N 8.9 [8.5,9.7]. The noise-only error of 0.3-0.7 K was too small.
- The reviewer's N excursion to 16-20 K is real, but only along one azimuth. CO brightens toward a peak 0.3" from the HCN centre.
- CO peak within 0.3" of the core at 0.26": S 33.4 K RJ (38.7 K Planck), N 21.2 K RJ (26.3 K Planck).
- Brunetti+2021 abstract (2021MNRAS.500.4730B): "as large as 37 K". This is reproduced within ~10% in either RJ (33.4 K) or Planck (38.7 K). The paper text I read does not state which convention it uses, so the check cannot say which.
- Table 2 comparison:
  - T_peak 16/50/84th percentiles at 55 pc are 6.5/13/18 K for S and 5.1/11/20 K for N.
  - The S core exceeds the 84th percentile by 1.7x. The N core sits in the 50th-84th range.
  - Table 2 is a distribution over the whole nuclear region, so the S core is its upper tail.

**3. Missing flux.**
- The 0.18" cube convolved to 0.82x0.68" recovers 0.86 (N, peak) and 0.91 (S, line mean ±100 km/s) of the 0.82" cube. The 12m-only 0.18" data are therefore 9-14% low at ~0.8", consistent with the earlier circular estimate of 10-12%.
- The 0.82" 12m cube convolved to the 7m beam recovers 0.83-0.89 of the 7m mosaic. The 7m mosaic convolved to the SD beam recovers 0.71-0.72 of SD.
- The additive feather correction (7m minus 12m, plus SD minus 7m) is +0.4-0.6 K at the line, i.e. +2% to +9% depending on beam. It does not include the 9-14% 0.18"-vs-0.82" deficit.
- T_b at 0.57" and above are lower limits by ~10-15%, and at 0.26" by about 5-10%.

**4. Conditional Tk (replaces Q23c).**
- Model: T_b = ff [J(Tex)-J(2.73)](1-e^-tau), with Tk >= Tex.
- At 0.26" with ff = 1 and tau >= 3:

| Core | Peak channel | HCN-FWHM mean |
|---|---|---|
| S | Tk >= 38 K | Tk >= 26 K |
| N | Tk >= 19 K | Tk >= 17 K |

- At tau = 3, Tk_min for S (peak/HCN-mean) is 71/48 K for ff = 0.5, 115/76 K for ff = 0.3, and 170/111 K for ff = 0.2. The N values are 34/30, 52/46 and 76/66 K for the same ff. At tau = 1 the limits are 40-50% higher.
- Gaussian-source fit of T_b versus beam, 0.18-1.7":
  - S: theta = 0.71 ± 0.08" (~150 pc), Ts = 33.9 ± 2.3 K, ff = 0.88/0.61/0.15 at 0.26/0.57/1.7".
  - N: theta = 1.1-1.3 ± 0.3", Ts = 11.6-13.4 K, ff = 0.95/0.79/0.3.
- So at 0.26" the CO core is nearly beam-filling.
- The statement "S cannot be a 17-25 K space-filling core" holds only conditionally:
  - It needs CO(2-1) thick and co-spatial with the dense gas at the same velocity, and no unmodelled 12m flux loss.
  - The velocity-matched limit (26 K) is only marginally above the 17-25 K S RADEX solution.
  - The peak-channel limit (38 K) is in tension if HCN fills the CO core.
- Q22 prior: Tk(S) >= 26-38 K at ff ~ 1 (>= ~50-70 K at ff ~ 0.5); Tk(N) >= 17-19 K.
- 13CO and C18O are not available, so tau(CO) is unconstrained.

## Literature
- Brunetti+2021 (2021MNRAS.500.4730B): abstract, section 4.3, Table 2.
- Table-2 HTML vs PDF: the arXiv text layer and ar5iv are the same source, so equality (12 rows x 15 values) is a null test. The typeset MNRAS PDF was not read separately.

## Uncertain
- Whether Brunetti's T_b is RJ or Planck.
- Co-spatiality of CO and HCN at the 0.05" level. The astrometry of the two cubes is assumed common.
- The N red wing above 2844 km/s is missing from the 0.18" cube.
- The missing-flux correction is additive and uniform over the beam, which is an assumption.
- The Tk limits depend on tau(CO), which is unconstrained.
- Dense-gas v0 and FWHM come from out/q5/q5e_width.csv.
