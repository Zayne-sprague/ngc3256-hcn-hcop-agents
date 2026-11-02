# Q35 round 2 — HCN-rich or HCO+-poor at N? (worker summary)

## What I did (reviewer fixes, in order)
1. **CS excitation contrast N/S** (`scripts/s210a_q35_cs_radex_grid.py` -> `out/q35/grid_cs.npz`; `s210c_q35_cs_radex.py` -> `out/q35/s210c.txt/json`; `s210d` thin tracers; `s210e` C18O check).
   - CS 5-4/2-1 (Harada 1.7", K km/s): N 0.315, S 0.202, C 0.180; N/S = 1.57 [1.28,1.94]. HCN 3-2/1-0 N/S 1.16, HCO+ 0.90. Also C18O 2-1/1-0 N/S 1.48, 13CO 1.34: the N excitation excess is shared by CS and the CO isotopologues, absent in HCN/HCO+ r31 (Q34). So it is a property of the gas (low-density/thin tracers more excited at N), not of CS alone.
   - Joint RADEX fit (HCN, HCO+, CS; common Tk, n, eta; Harada 6 lines + 4.3" CS7-6/HCN4-3 and HCN4-3/HCO+4-3): N Tk 199 [178,223] K, log n 4.50; S Tk 100 [63,158], log n 4.62; Bridge(C) 178 K, 4.25. CS Tex(2-1/5-4/7-6): N 11/9/17 K, S 13/8/14, C 7/9/17 (sub-thermal): the contrast is Tk/tau, not Tex.
   - Tension test: HCN+HCO+ (Tk,n) predicts CS 5-4/2-1 = 0.10 [0.04,0.23] at N vs 0.315 observed; S 0.10 vs 0.20 (OK). One-gas best-fit chi2 10.2 (N, ~0 dof; 11.6 with dof 2) vs 2.1/3.5 at S -> CS is a questionable abundance reference at N.
2. **Column densities / abundances** (same scripts). Adopted fit, N/S: N(HCN)/N(HCO+) 1.64 [1.31,2.09]; N(CS)/N(HCN) 0.96 [0.71,1.28]; N(CS)/N(HCO+) 1.58 [1.19,2.12]; N(HCO+)/N(CS) 0.63 [0.47,0.84] (P>1 = 0.06). Bridge/S ~1 for all. Thin limit (`s210d`): H13CN/H13CO+ 1-0 N 0.75, S 0.23 [0.09,0.42] (truncated), N/S = 3.3 [1.7,9.0], P>1 0.97, 5th pct 1.19 — weak (S H13CN is 1 sigma); with 12C/13C 30-100 (cancels) + 0.15 dex per-species fractionation N/S of the column ratio 3.6 [1.3,11.8]. The previous "N/S >1.0 at 3 sigma" is withdrawn.
3. **Systematics** (`s210b_q35_injection.py` -> `out/q35/s210b_box.txt`, `s210b_field.txt`; `s210d`).
   - Monte Carlo injection into the beam-shaped deblend: bias <=5% in 4 scenarios; true scatter of the double ratios is 1.5-2x the old sqrt(8) estimate: HCO+/cont N/S 0.67 [0.56,0.78] (2.8 sigma, was "4 sigma"), HCO+/CS 0.68 [0.45,0.96] (1.2 sigma), HCN/CS 1.10 [0.72,1.67], HCN/HCO+ 1.62 [1.23,2.15].
   - Free-free from H40a (Harada): 2-5 mJy/beam at both N and S (5-12% / 9-20% of continuum); dust-only N/S 1.9-2.0 (2.3 generous) -> HCO+/cont 0.64-0.52, HCN/cont 1.05-0.84. HCN/HCO+ unchanged. Dust Td/kappa difference N vs S not constrained.
   - Velocity-resolved deblend (100 km/s): HCN and HCO+ 4-3 single peaked, same widths at N and S, HCN/HCO+ constant with v (N 0.55, S 0.32-0.36): no self-absorption dip seen at this resolution.
   - Alternative normalisers (Harada 1.7"): HCN/HCO+ N/S vs C18O 1-0: 1.48/0.99; CS 2-1: 1.63/1.09; C18O 2-1 with 3-2: 1.16/0.60; CS 5-4: 1.20/0.62; 13CO is thick (N/S 0.96) and unusable.
4. **Consistent conclusion + plot**: `s211_q35_plot.py` -> `out/q35/s211_ns_vs_J.png`; report entries rewritten by `s212_update_report_q35.py` (backup `/workspace/report_backup_q35r2.json`; validator passes).

## Result
HCN/HCO+ is 1.5-1.9x higher at N at every J. The HCO+/CS sign flip (1.09 -> 0.63-0.68) comes from CS(2-1) vs CS(high J) because N is more excited in low-density tracers (C18O, CS) while HCN/HCO+ r31 are equal at N and S. Excitation-corrected columns with CS as reference: HCN tracks CS (N/S 1.04), HCO+ is low at N (0.63, ~2.5 sigma) -> "HCO+-poor", not "HCN-rich"; low-J references give the opposite. Confidence **low**: CS is excitation-discrepant at N (chi2 10), the continuum deficit is 2.8 sigma, Td differences unconstrained. Disk (W, Sblue) has only CS limits; C/OS/SW at 1-0/2-1 do not depart from S except OS x1.6 in both ratios.

## Literature
Harada+2018 (Table 13 line intensities, parsed in `lit/harada_table13_parsed.csv`); it gives no abundances for NGC 3256, so I compared line ratios only. Rohlfs & Wilson RRL/continuum relation for H40a free-free; LAMDA (CS, HCN, HCO+, CO) for RADEX.

## Remaining uncertainties / not done
- Harada cubes are not in /data, so the 1.7" and 4.3" sets are not on one beam; the 4.3" constraints enter as line ratios. A convolved CS 2-1/5-4 at 4.3" would need the Harada or ALMA archive cubes (CS 2-1 Band 3, CS 5-4 Band 6 12m+7m restores).
- RADEX grid cells are 0.2 dex (posteriors jittered +-0.1 dex); MC bias scenarios are idealised; H13CN blend check (SO/SiO/CH3OH) only qualitative; H13CN 4-3 not used.
- Dust temperature difference N vs S (3 mm continuum or 345 vs 230 GHz spectral index) would test the continuum normaliser; not done.
- C18O 2-1/1-0 observed (1.1-1.6) is far below the dense-gas RADEX prediction (3.7-4.0): CO isotopologues come from lower-density gas; used only for the sign of the N-S contrast.
