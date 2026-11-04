# r10 — answer (team-a liaison)

Summary: (a) the 1.48-1.49 is Harada's own tabulated number, not something our team measured; the script below recomputes it from the Table 13 entries, and our own matched-beam value is 1.41-1.45. (b) Neither requested quantity exists in our work, and I did not invent them; the closest outputs are listed. (c) Bibcodes are given; two are inferred from the standard format and marked as such.

## (a) Harada et al. 2018 Table 13, HCN/HCO+ 1-0, N/S at 1.7"

Script and output: `r10a_harada_ratio.py`, `r10a_harada_ratio.txt` (same directory).

Important: the 1.48/1.49 in p0058 is a **published** number (Harada+2018, ApJ 855, 49, 2018ApJ...855...49H, Table 13). p0058 (team-b) quotes it; neither team re-measured it at 1.7". The script only divides the tabulated entries. Cite Harada Table 13, not our script, as the source of the 1.7" value.

Inputs (Table 13, cubes convolved to 1.7", ∫ over ±150 km/s about 2775 km/s, K km/s; Δ(I∆V) in brackets):

| | HCN 1-0 | HCO+ 1-0 | HCN/HCO+ |
|---|---|---|---|
| N | 197.23 (0.96) | 205.25 (0.95) | 0.961 ± 0.006 |
| S | 56.68 (0.96) | 87.87 (0.95) | 0.645 ± 0.013 |

N/S = 0.961/0.645 = **1.49 ± 0.03** (quadrature, statistical errors only as tabulated). The 1.48 follows from dividing the two-decimal ratios 0.96/0.65 = 1.477, so the 1.48 vs 1.49 difference is rounding. The peak-intensity version (Ipeak, mK) gives 0.947/0.660, so N/S = 1.43. The ratio is a T_B ratio, since both lines are in K.

Our own measurement (not Harada's): 2016.1.00965 X94a spw27, HCN and HCO+ from the same cube, fixed window ±350 km/s, T_B ratio = Jy ratio × 1.0126 (`scripts/q55f_r10_core.py`, `out/q55/r10_core.csv`, quoted in p0058-r01). We have no 1.7" beam; the neighbouring beams are:

| beam | N | S | N/S |
|---|---|---|---|
| 1.4" | 0.974 | 0.693 | 1.41 ± 0.13 |
| 1.9" | 0.933 | 0.654 | 1.43 ± 0.10 |
| 2.5" | 0.899 | 0.623 | 1.44 ± 0.09 |
| 1.0" | 1.059 | 0.731 | 1.45 ± 0.22 |

Interpolating to 1.7" gives about 1.42, consistent with Harada's 1.49 ± 0.03 within our errors. team-b's 2.5" value is 1.56 ± 0.03 (p0058, different cubes, 2015.1.00993). Differences in window (±350 vs ±150 km/s), positions and baseline are not controlled, so the agreement is at the 5-10% level, and the paper should say "consistent with", not "reproduces". The 1.4-2.5" ratios are the 1-0 line, a different quantity from the 3-2 core contrast of ~2.1.

## (b) Xb4b 0.30" CLEAN-only vs restored core flux; per-channel errors on velocity-resolved ratio

**CLEAN-only vs restored core flux for Xb4b at 0.30": not available.** The Xb4b 0.30" maps were made from the archive pipeline cubes (`scripts/q64r2c_images.py`, `out/q64r2/`: Xb4b convolved to 0.3-1.5" circular beams). Those products contain no CLEAN model, and we never re-imaged Xb4b from visibilities. The raw-data imaging with model files (`scripts/q64c_image.py`, `out/q64/img/*.model.fits`) covers Xb4d (12-m, 0.9") and 12-m+ACA, not Xb4b. Producing the split would need a new CASA run on Xb4b visibilities, which I did not start; it is outside what this request needs to cite. Do not use any CLEAN/restored split for Xb4b in the paper.

What exists instead, for the Xb4b 0.30" HCN/HCO+ 3-2 ratio (`out/q64r2/summary.txt`, `profile.txt`; pipeline cubes, median over baseline windows, [min,max]): N 0.986 [0.932, 1.006], S 0.462 [0.448, 0.485], N/S 2.13 [1.92, 2.25]. Uv-plane model fits give N 0.92-1.12, S 0.43-0.46 at 0.4-0.6" (models A/B, same file), which is a spread over models, not a CLEAN-only vs restored comparison. Per-line fluxes are in `out/q64r2/img_ratios.csv`.

**Per-channel errors on the velocity-resolved ratio, 0.95-1.09 range:** I could not trace the source of this range in our outputs or on the forum, so I cannot supply errors for those numbers. Please tell me which post and figure it comes from. The only velocity-resolved HCN/HCO+ 3-2 profile in our work is `out/q1/q1_ratio_vs_velocity.csv` (`scripts/q1_velratio.py`). It is Xb4d at a common 1.05" beam, not the 0.30" Xb4b maps, in 50 km/s bins about v_sys, with errors propagated from channel noise (statistical only). The ratio there is not 0.95-1.09:

| v − v_sys bin (km/s) | N R ± e | S R ± e |
|---|---|---|
| −150 to −100 | 0.818 ± 0.055 | 0.365 ± 0.015 |
| −100 to −50 | 0.786 ± 0.008 | 0.386 ± 0.006 |
| −50 to 0 | 0.802 ± 0.007 | 0.413 ± 0.009 |
| 0 to 50 | 0.808 ± 0.006 | 0.360 ± 0.012 |
| 50 to 100 | 0.836 ± 0.008 | 0.361 ± 0.014 |
| 100 to 150 | 1.270 ± 0.077 | 0.447 ± 0.056 |

Bins at |v| > 150 km/s are lower limits (HCO+ below 3σ in the N bins at 150-250 km/s) and the S bin at −200 to −150 is 0.603 ± 0.120. The 100-150 km/s N bin is the red wing, where the HCN v2=1f line (+401 km/s relative to HCO+) was checked in the same script and where HCO+ is weak, so treat the 1.27 with care. Errors are statistical only and do not include baseline or bin-edge systematics (the integrated ratio varies by ≤0.01 with window in `q1_velratio.log`).

## (c) References

- **Tang et al. 2019**: Tang, X. D., Henkel, C., Menten, K. M., et al. 2019, "ALMA view of the 12C/13C isotopic ratio in starburst galaxies", A&A 629, A6, DOI 10.1051/0004-6361/201935603. Bibcode by the standard ADS pattern: **2019A&A...629A...6T**. I could not open ADS to confirm the string; the title, volume, article number and DOI were confirmed from the A&A page via web search. Team-b cites it for the 12C/13C prior (p0044, p0050, p0061). Note it measures 12CN/13CN in NGC 253, NGC 1068 and NGC 4945, not NGC 3256, so say "starburst-nucleus range" in the paper.
- **Falstad et al. 2021**: Falstad, N., Aalto, S., König, S., et al. 2021, "CON-quest: Searching for the most obscured galaxy nuclei", A&A 649, A105, DOI 10.1051/0004-6361/202039291, arXiv:2102.13563 (journal ref and DOI confirmed from the arXiv record). Bibcode by pattern: **2021A&A...649A.105F** (not opened on ADS). Our team used it only for the CON criterion, via p0057/p0059; it was not read in full by us.
- **Jimenez-Donaire et al. 2017**: Jiménez-Donaire, M. J., Cormier, D., Bigiel, F., et al. 2017, "Optical depth estimates and effective critical densities of dense gas tracers in the inner parts of nearby galaxy discs", MNRAS 466, 49, bibcode **2017MNRAS.466...49J** (bibcode, title and first author confirmed in the paper library; volume/page follow from the bibcode; the full author list is from memory, so check it).
- **Harada et al. 2018**: ApJ 855, 49, **2018ApJ...855...49H** (confirmed in the paper library).
