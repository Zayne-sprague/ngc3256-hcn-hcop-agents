# Q6 work summary (worker, team-b): sizes, line widths, surface brightness, Tb and beam filling at 0.24" (Xb4b)

## What was done
All on the Xb4b pbcor products (HCO+ 3-2 spw25, HCN 3-2 spw27 cubes, 4.4 km/s channels, beam 0.246x0.232" / 0.246x0.237"; 4-spw continuum 257.76 GHz, beam 0.245x0.235"), moment-0 over |v-2780|<350 km/s radio-LSRK (per-pixel line-free baseline). Distance 40 Mpc, 1" = 194 pc (team convention; Sakamoto+2014 use 35 Mpc, 170 pc/": pc values scale by 0.88).

Scripts (python3 -I, re-runnable, outputs in /workspace/out/q6/):
- `scripts/s150_q6_sizes.py`: beam-deconvolved 2D elliptical Gaussian + constant fits, boxes 1.2/2/3", Condon-style covariance inflation (s150.json).
- `scripts/s151_q6_size_width_tb.py`: non-parametric half-max-area size (noise-injected errors), peak-pixel and r<0.25" spectra, half-max/W20/1-Gaussian widths, velocity-gradient beam smearing, T_b (RJ and Planck-equivalent) (s151.json; `s151_print.py` prints).
- `scripts/s152_q6_aperture_flux.py`: aperture fluxes r<0.5/1/2/4" at N, S in Xb4b/Xb4d/Xb4f (s152.json).
- `scripts/s153_q6_radex_filling.py`: RADEX LVG T_R(3-2) table and thermalised-limit filling (s153_radex_filling.csv).
- `scripts/s154_q6_table.py`: assembled table (s154_table.txt), dilution-corrected T_b, virial-style masses.
- `scripts/s155_q6_cont_tb.py`: continuum T_b and dust-only surface density (s155.json).
- `scripts/s156_q6_injection.py`: injection test of the half-max size estimator (unbiased to +-0.03" for D>=0.15"; +0.05" for an unresolved source).
- `scripts/s157_update_report_q6.py`: writes 6 Q6 results + recipe + caveat to /output/report.json (validated against the schema).

## Results (details and uncertainties in report.json, entries starting "Q6:")
1. **Sizes (deconvolved half-max diameters)**: HCO+ 3-2: N 0.52" (100 pc), S 0.50" (97 pc), Bridge 0.82" (160), E1 0.69" (135), E2 0.60" (116); HCN 3-2: N 0.41" (80 pc), S 0.43" (83), Bridge 0.70", E1 0.54". Continuum: S 0.23" (45 pc, marginally resolved), SE 0.14", E2 0.28", E1 0.48", N 0.81" (extended, faint). Nuclei are resolved (about 2 beams); HCN is more compact than HCO+ at N (ratio 0.79+-0.03); the S continuum is half the size of the S line core. Both nuclei also carry a 1.5-2" pedestal; single Gaussians are box-dependent and the core+pedestal 2-Gaussian fit is unconstrained (not used). Core (r<0.5") holds 16-31% of the line flux inside r<4", r<1" holds 28-68%.
2. **Line widths (peak pixel, native beam)**: HCO+ FWHM (half-max/1G): N 111/132, S 93/99, Bridge 86/92, E1 59/60, E2 46/54 km/s; HCN N 158/162, S 98/114, Bridge 85/82. Off-nuclear clumps have sigma 23-28 km/s = the CO cloud median 23 km/s of Brunetti+2022 at similar size (D 116-135 pc vs their R 88 pc), nuclei are 1.7-3x broader. Beam smearing from the measured velocity gradients is only 4-16 km/s, so nuclear widths are intrinsic at 50 pc. Virial-style masses (1040 sigma^2 R): clumps 5-9e7 Msun (Brunetti median 5.4e7), nuclei 1.5-3e8 (upper bounds; rotation/outflow inflate sigma).
3. **Tb and beam filling**: peak T_b(RJ) HCO+ N 2.75, S 2.81, Bridge 1.87, E1 2.05, E2 2.04 K; HCN N 2.46, S 1.18 (Planck-equivalent 5-7.4 K, = Tex lower limit). Beam dilution for the resolved structures is only 1.2-1.35. Compared with the thermalised limit (24/54/114 K for Tk 30/60/120), filling >=0.12/0.05/0.02 (HCO+ S); with n=1e5, N=1e15 (RADEX) filling 0.56. HCO+ 3-2 is 0.075 of the CO 2-1 peak (37 K, Brunetti+2021). So the emission is sub-thermal and/or fills <~10% of the 50 pc beam; filling and excitation are degenerate (no second line at 0.25").
4. **Missing flux**: Xb4b recovers Xb4d at r<2" (HCO+ N 1.00, S 0.83; HCN N 0.97, S 0.81), 0.8-0.87 at r<4"; vs Xb4f r<4" 0.75-0.86 (lines), 0.65-0.76 (continuum); consistent with p0003. Compact cores are unaffected at the 10% level; the pedestal is a lower limit.
5. **Continuum**: S peak 0.53 K, N 0.13 K; dust-only upper bound Sigma ~2-4e3 (S), 0.5-1e3 (N) Msun/pc2 (kappa and Td assumed, low confidence), same order as CO-based nuclear Sigma (Sakamoto+2014, Brunetti+2021).

## Literature used
Sakamoto+2014 (2014ApJ...797...90S): N disk ~2" diameter, S compact ~0.5" peak (angular sizes agree with ours), M_dyn(S, r<200 pc) ~2e9 Msun; Brunetti+2021 (2021MNRAS.500.4730B): CO 2-1 0.25", peak Tb 37 K; Brunetti+2022 (2022MNRAS.515.2928B): cloud sigma 23 km/s, R 88 pc, M_vir 5.4e7 (taken from abstracts/passages, not re-derived); Rosolowsky & Leroy 2006 virial convention.

## Uncertain / not done
- Sizes are statistical-error-dominated by non-Gaussian structure and blending (systematic 0.05-0.1"); Sblue, W, SE and HCN E1 have S/N<10 (tentative).
- Moment-0 beam taken from the cube header (one beam; no per-channel table in Xb4b).
- Filling factor cannot be separated from excitation with one transition at 0.25". Only a matched-beam Band 7 / Band 5 (4-3, 2-1) at this resolution would do it; not in the data.
- Xb4b manual (26.6 MHz) cube was not used as an independent imaging cross-check of the sizes (cheap follow-up).
- Distance scale (40 vs 35 Mpc) scales pc values by 0.88; no ramp (Q23) affects spw25/27 in Xb4b.
