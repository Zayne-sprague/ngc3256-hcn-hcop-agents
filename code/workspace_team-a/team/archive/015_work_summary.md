# Q32 work summary — is the compact CO(2-1) core resolved at 0.26"? (ff range for the Tk bound)

## What I did
No new data. Used the Q23 matched-beam T_b of the 0.18" 12m CO(2-1) cube (2015.1.00714.S, convolved to 0.26-1.7"; one cube, so uv coverage is identical across beams) and tested the 0.82" cube (2015.1.00902.S) for missing flux.
- `scripts/q32_resolved.py` -> `out/q32/curve_fit.json/.log`: fit T_b(theta)=Te+T0 s^2/(s^2+theta^2) (Gaussian core + pedestal), MC for ff.
- `scripts/q32_radial.py` -> `out/q32/radial.json/.png`: radial profile of the native map about the peak vs the CLEAN beam (model-independent resolution test).
- `scripts/q32_tk.py` -> `out/q32/tk_ff.json`: Tk limits with fitted ff, tau from Q30 (R=40-90); 11% missing-flux variant.
- `scripts/q32_map_fit.py` (2D Gaussian+const over +-1.5") is NOT used: it locks on the larger clump (0.8-1.8"), not the core. Kept for transparency.
- `scripts/build_report_q32.py` -> /output/report.json (Q32a-c; 200 results).

## Results
- Resolved: native-map radial HWHM 0.27" (S) / 0.21" (N) vs beam 0.09".
- Gaussian+pedestal fit: S s=0.47(0.43-0.52)" (100 pc), T0=33(30-36) K RJ, Te=4.4 K; N s=0.27(0.22-0.32)" (57 pc), T0=13.5(11.6-16.1) K, Te=6.1 K.
- ff(0.26")=0.77 (0.73-0.80) S, 0.52 (0.42-0.60) N; ff(1.0")=0.21 S, 0.08 N.
- Tk lower limits with measured ff and tau(R=40-90): S>=44 K at 0.26" (38.4-38.9 K for the deconvolved compact T0 alone), N>=29-30 K at 0.26" (18-19 K for deconvolved T0). Replaces ff=1 S>=37, N>=19 K. ff=0.5/0.3 (S>=68/110 K) not supported at S.
- 13CO: beam 0.85" > core, ff of 13CO untestable (Q32c).

## Literature
Brunetti+2021 Table 2 (cloud-scale T_peak, used in Q23); Pereira-Santaella+2024 (JWST, S nucleus) only noted. No published CO core size for NGC 3256 S/N at 0.2" found in the lit.py search.

## Uncertain
- Pedestal vs no-pedestal model changes s (0.47 vs 0.75" S; N 0.27 vs 1.4"); the pedestal model fits far better (chi2 0.4 vs 7.5) but the pedestal is assumed constant and additive.
- Peak-channel (~2.5 km/s) extent; HCN-window mean is 0.64x lower in T_b at S (Tk>=26 K at ff=1).
- Sub-0.18" clumping not excluded; tau bound conditional on R>=40 (Q30); 12m-only missing flux (10-12% at >=0.75").
