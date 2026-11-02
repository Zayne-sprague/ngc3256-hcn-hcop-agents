# Q34 work summary: uv-matched test of N/S HCN/HCO+ and r31 (1-0 vs 3-2)

## What was done
- Restores (rawdata.py): Band 6 2017.1.00379.S Xb4d 12-m EB Xe07 (0.81 GB split MS) and 3 Xb4f 7-m EBs (0.28 GB). Band 3 visibilities could NOT be restored (X961/X145: "no pipeline manifest"; X94a 57.7 GB > 30 GB limit), so the 1-0 side stays the short-spacing-complete archive cubes (X145 1.5", X961 2.2"; agree <=5%).
- s191_q34_casa_image.py (CASA on host): HCN/HCO+ 3-2 cubes, mosaic gridder, natural weights, uvtaper + restoringbeam = circular 1.5" and 2.2", 20 km/s, pbcor; variants `comb` (12m+7m) and `c12` (12m only). Outputs /workspace/out/q34/*.pbcor.fits (+.image, .pb).
- s192_q34_measure.py: p0034 conventions (r = beam, +-350 and +-150 km/s, RJ T_b, line-free noise) -> out/q34/s192.txt/json.
- s193_q34_radex.py: s186 RADEX evidence framework with the new 3-2 -> out/q34/s193.txt/json. s194: report.json (validated, 0 errors).

## Results
- N/S contrast of HCN/HCO+ at 2.2": 1-0 1.467+-0.022; 3-2 (12m+7m) 1.441+-0.023 (Xb4d image alone 1.650). D = 0.983+-0.022 (+-350), 1.038+-0.016 (+-150); at 1.5": 0.948+-0.038 / 1.011+-0.029. Adopted D = 1.00+-0.04. The earlier D = 1.10+-0.06 (p0034) was caused by 12-m flux loss that is larger at S.
- Flux recovery (12m+7m / Xb4d image, 2.2"): HCN N 1.06, S 1.34, Bridge 1.15; HCO+ 1.07, 1.17, 1.10: a measured two-sided correction, in agreement with the Xb4f/Xb4d bound (S 1.30-1.36).
- r31 at 2.2": HCN N 0.263, S 0.261, Bridge 0.258; HCO+ 0.351, 0.343, 0.333 (abs. +-14%). N/S of r31: HCN 1.007+-0.020, HCO+ 1.025+-0.010 (previously 1.27 / 1.09): no excitation difference N vs S. r31(HCN)/r31(HCO+) ~ 0.75 everywhere.
- RADEX: posteriors (Tk, n, A) essentially unchanged; N vs S ln BF abundance-only (same Tk,n) +1.0 vs excitation-only (same A) +0.2 at 2.2" (+0.3 vs -0.4 at 1.5"): weak (<1) preference for abundance over excitation, not decisive. Pure "filling only" is disfavoured less strongly than before (-4.6 vs -14 at 2.2").
- If one stays within 2 sigma: nothing requires an excitation difference; the N/S HCN/HCO+ difference is present equally at J=1-0 and 3-2.

## Caveats
- Only 1 of 2 Xb4d 12-m EBs and 3 of 14 Xb4f 7-m EBs restored; Band 3 not re-imaged from visibilities and no common kλ cut with Band 3 enforced; no taper-to-Band-3 whole-disk r31 control. CASA printed "cannot convolve to new beam" warnings in 6 of 8 runs; final cubes have the exact target beam.
- Window dependence (350 vs 150 km/s) of D is 0.05 at 2.2", included in the error.
- 1.0" matching with X94a not attempted (resolves out 30-50%); an X165 (2015.1.00412) independent 3-2 set and the X94a visibilities would be the next steps.

## Literature
Used Harada+2018 (Table 13 1-0 ratios, via Q29), Michiyama+2018 (via Q29); no new literature search was needed for this uv-matching test.
