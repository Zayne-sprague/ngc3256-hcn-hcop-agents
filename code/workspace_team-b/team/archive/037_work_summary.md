# Q9 work summary (worker, team-b)

## What I did
Re-examined the CH3OH evidence behind p0001 (team-a) and our p0016: the Xb55 163.5 GHz J_-1-J_0 E blend (J=1-3), the 254.015 GHz 2(0)-1(-1)E line, Trot and the "Q-branch excess". Most earlier extractions (s17-s27, s97d) were reused; new work addresses the weakest link, the blend integral's baseline systematic.

## Scripts / outputs (re-runnable, `python3 -I`, in order)
- scripts/s120a_q9_spectra.py -> out/q9/s120a_spectra.png (Xb55 spw22 spectra at N,S,Bridge,E1,E2,Sblue,W; rms 2.1-2.5 mJy/beam per 1.95 MHz channel)
- scripts/s120b_q9_blend_integral.py -> out/q9/s120b.{txt,json} (model-independent blend integral, 5 baselines x 4 windows)
- scripts/s120c_q9_trot.py -> out/q9/s120c.json (Trot Monte Carlo from W254/Wblend, baseline range propagated)
- scripts/s120d_q9_blend_vs_hco.py -> out/q9/s120d.json (blend/HCO+ 2-1, Xb55 aperture, K km/s)
- scripts/s120e_q9_column.py -> out/q9/s120e.json (beam-averaged E-CH3OH column)
- scripts/s121_update_report_q9.py -> /output/report.json (now 180 results; 3 new "Q9:" entries, recipe, caveat)

## Results
1. **Blend integral (Jy/beam km/s, +-400 km/s, aperture)** depends on the baseline by a factor 2-4 (red wing + negative ramp toward the spw edge): N 0.46 (blue-side const) .. 1.6 (linear); S 0.86..1.8; Bridge 0.78..2.0; W 0.5..1.0; statistical error only 0.1-0.15. Adopted medians N 1.3, S 1.3, Bridge 1.6, W 0.86. team-a's 2.23 (N) and reviewer's 1.84 sit at the linear-baseline end; our earlier 0.97/1.38 is mid-range. The earlier "+-30% systematic" was too small. The line itself is real (4-6 sigma per channel at S/Bridge/N; absent at E1, E2, Sblue).
2. **Trot** (254/blend, thin LTE, Monte Carlo over baseline): S 10 (+7/-4) K; Bridge median 22 K but only >~7 K; N unconstrained (254 not detected). "Trot~9 K" is reproduced as a central value but is not a constraint. Beam-averaged N(E-CH3OH) ~1e13 cm^-2 at S (0.8e13 Bridge), x2 for A+E, no beam-filling correction.
3. **Q-branch excess / chi2=53/6**: not confirmed. Channel-space LTE fit (s18) at 8-10 K fits N and S (chi2 285/398, 303/396 channels); p0001's own J=4-8 single-line amplitudes are negative (J=5 -0.37+-0.08, -4.4 sigma), a baseline artefact, and their blend sits at the high-baseline end; a factor-2 swing in the blend equals the claimed 1.6x. Individual J=1,2,3 amplitudes are degenerate (only the sum is measured).
4. **Contrast**: at Xb55 (9") blend/HCO+ 2-1 = 0.05 (N), 0.047 (S), 0.051 (Bridge), 0.063 (W) (K km/s, errors +-40-80% from baseline) -- no contrast; the beam holds N, S and Bridge together. The Bridge/W enhancement (0.04-0.07 vs N<0.022) is a 1-4" result from the 254 line only, 3-4 sigma, single transition; it is not resolved-out (Xb4d/Xb4f recovery 0.77+-0.3) and partly reflects CH3OH being extended while HCO+ is compact.
5. **Kinematics**: CH3OH 254 centroid vs HCO+ 3-2: S +16+-5, Bridge +10+-7 / +7+-2 km/s, Sblue +43+-23 (2 sigma); no outflow-velocity offset; NS excluded.
6. Band 7 CH3OH coverage: no Eu<100 K, logA>-4.5 line available beyond those already tested (s23: 266.838 and the 251.7-251.9 A ladder; the latter is blended with SO 251.826 and gives only limits).

## Uncertain
- Which baseline is right: the ramp near the spw22 edge is unexplained (continuum-subtraction residual or spw edge); the true blend lies 0.5-2 Jy/beam km/s.
- LTE only; RADEX (s18) allows 8-40 K for n>=1e5 and excludes n<=1e4 for the blend alone. No second independent Eu constraint beyond 254 and the limits J=4,5.
- N/S contrast needs a second transition at <=4" resolution (none exists in the data); Q_E partition function approximate (+-30%).
- I did not re-derive team-a's centroid test for 251.66 independently of s20 (reused: CH3OH +7..+16 vs NS -38..-47 km/s).
