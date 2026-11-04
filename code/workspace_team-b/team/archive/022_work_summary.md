# Work summary: Q1 round 2 (worker, team-b)

Responds to the review in team/review.json. Disagreement noted in advance: round-1 deblend (s31) already used line-specific templates (T[mol]); the problem was misfit, not a shared template.

## What I did (scripts in /workspace/scripts, outputs in /workspace/out; all re-runnable with `python3 -I scripts/<name>.py`)
1. **Disk / radial trend** `s36_disk_annuli.py` -> `out/s36_disk.txt/json`. 10" beam, annuli about the N-S midpoint, Xb55/Xb4f/Xb4d/Xb47, Gaussian 12-m PB model (line-free rms of these pbcor cubes does not rise as 1/PB, so it is not usable as PB), far-field zero level, continuum-vs-line annulus fractions, sectors. Go/no-go: PB>=0.5, flux >3 sigma both lines, bowl <30%.
2. **Deblend v2** `s37_deblend2.py` -> `out/s37_deblend2.txt/json`. Line-specific templates, signed N/S/Ext + flat bowl + 20" Gaussian, 150-draw template bootstrap, per-zone ladder ratios vs Xb4f at the same beam, injection test with real correlated noise.
3. **Sblue wing** `s38_sblue_wing2.py` -> `out/s38_sblue2.txt/json`, `out/s38_spectra.png/txt`. Null (400 random + bright positions), posterior, windows/apertures/centre shifts, 50 km/s bins, Xb4b-smoothed and Xb4f replication, neighbouring lines, peak-pixel bias injection, RADEX abundance consequences.
4. **Ratio maps and profiles** `s39_ratio_maps.py` -> `out/s39_ratio_maps.png/txt/json`.
5. **Report** `s40_update_report_q1r2.py` -> /output/report.json (104 results; Q1 entries amended, v1 deblend and v1 Sblue entries replaced, 5 new 'Q1r2:' entries, interpretation and caveats rewritten).

## Results
- **N/S radial profile (3-2):** N ratio rises inward: Xb4b annuli 0.39 (2-3") ,0.41 (1-2"), 0.65, 0.87, 0.97 (<0.25"); cumulative 0.57 -> 0.97. Xb4b and Xb4d fluxes agree to 3-10% inside 2", so this is beam dilution of a compact (<~0.5", ~100 pc) HCN-bright core, not filtering. S is flat at 0.37-0.43 from 0.25" to 3" in all arrays (Xb4f 0.38-0.40). HCN and HCO+ peaks coincide with the continuum peaks to <=0.07". Outflow-side region (Sblue r<1.5", |v|<350): 0.51+-0.03 vs S r<1.5" 0.37+-0.01.
- **Disk:** only Xb55 (2-1) is a GO: 0.575 (0-4"), 0.554 (4-8"), **0.489+-0.004 stat (+-~0.1 sys) at 8-12"**, i.e. ~12% below the interior. Annuli beyond 12" have negative net flux (bowl) in all arrays. The round-1 "0.30" (8<r<20") was an artefact of the negative 17-24" ring and is withdrawn. Xb4f/Xb4d/Xb47 are GO only for r<8" (ratio drops 10-20% from 0-4" to 4-8" at every J). 3-2 and 4-3 disk ratios are unconstrained (missing short spacings vs excitation degenerate). Fraction of r<8" flux in the 8-12" annulus: continuum 0.41/0.27/0.27/0.14 (B5/B6f/B6d/B7), HCO+ 0.33/0.23/0.20/0.09: same ordering = increasing resolve-out. Clumps E1/E2/W: R(2-1) 0.51-0.54, R(3-2) 0.25-0.35.
- **Deblend v2** (chi2_red still 9-99, so errors are inflated and the model is approximate): R(2-1) N 0.64+-0.12, S 0.46+-0.08 (10"); R(4-3) N 0.46+-0.14/S 0.35+-0.14 (10"), N 0.46+-0.07/S 0.35+-0.07 (6.5"); S/N ~0.7+-0.2 at all J. Zone ladders (Jy): 2-1/3-2 N 0.79/0.72 (HCN/HCO+), S 0.96/0.90 (+-0.15-0.2) -> T_b ratio 1.6-2.2 (>1); 4-3/3-2 N 0.81/1.02, S 0.74/0.93 (+-0.2-0.3) -> T_b 0.41-0.57. Injection: S coefficient biased high vs N by 10-25% -> adopt 15-25% systematic. Ext is unconstrained.
- **Sblue wing:** HCN 0.449+-0.062 (7.3 sigma), HCO+ 0.320+-0.060 (5.3 sigma); R=1.40, posterior 68% 1.12-1.79, 95% 0.90-2.35; P(R>1)=0.93; R_wing-R_core>=+0.36 (95%). 0/400 random positions have both wings >3 sigma; bright positions show no comparable wing. Window dependence: -350..-100 gives 1.0; R rises with blueshift (0.55 core -> ~1 at -200..-150 -> 1.3-2.6 beyond, bins only 1-3 sigma). **Replication weak:** Xb4b smoothed to Xb4d beam 0.86+-0.25 (HCN recovered 53%, HCO+ 91%); Xb4f blended 0.83+-0.25. So R_wing>R_core is robust, R_wing>1 is not confirmed. The peak-pixel bias statement about p0001 is dropped (injection shows no bias). RADEX: same excitation needs A(HCN/HCO+) x2.6; a denser/hotter wing (n 1e6->1e7) reproduces x2.6 without abundance change: degenerate.
- **Systematics** are separated in a report entry: within-band R(J) is calibration-free (5%); R(J) trends add ~10%; flux ladders add 10-35% (no overlapping line in Xb55/Xb47 to test directly).

## Still uncertain
- Disk ratio zero level (+-0.1); no 3-2/4-3 disk ratio.
- Deblend model misfit and S-vs-N bias; Ext zone free.
- Sblue R>1 not replicated by Xb4b; HDO/c-C3H possible blends in the HCN window; spectra at out/s38_spectra.png not visually inspected by me beyond numbers (bins above).
- HCN-vib and opacity for the N compact core (Q4) not tested; T_b(2-1)/T_b(3-2)>1 may reflect missing 3-2 flux.
- Forum: not posted (nothing outward-facing done).
