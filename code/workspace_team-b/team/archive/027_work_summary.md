# Q42 work summary (worker, team-b)

## What I did
1. Added the H13CN(4-3)/HCN(4-3) ratio to the Q40 one-zone RADEX fit as a same-line 13C/12C ratio term, with N(H13CN)=N(HCN)/X and the same X nuisance as Q40 (X marginalised, band-scale priors 10/15, 20/25, 40/40 %). Three data versions: censored (alpha = 0 +- 0.04, i.e. the 3-sigma limit 0.12 at N from Q4/p0007), fitted (0.045 +- 0.04), and a dilution variant (f = 0.35, sigma/f). H13CN was modelled with its own LAMDA file (h13cn@xpol.dat; new grid).
2. Re-extracted/tabulated H13CO+/HCO+ and H13CN/HCN at N, S, C at ONE common 2.2" taper (Q29 s185 products; native beam 2.1x1.7", so 1.7" itself is unreachable) in three extractions (template aperture, +-150 km/s aperture, +-150 km/s peak pixel) and put them into the RADEX fit as ratio terms. This gives X_N, X_S, X_C at one taper and the S/N contrast, compared with the Q40 numbers (Harada 1.7" at N vs 2.2" at S).
3. Checked what the common-taper 13C data do to the Q40 conclusion (cold/dense thick branch at N).

## Scripts and outputs (re-runnable, `python3 -I`)
- scripts/s250a_q42_h13cn_grid.py -> out/q42/grid_h13cn.npz (RADEX H13CN 1-0/3-2/4-3, same grid as s220a)
- scripts/s250b_q42_ratio_fit.py (fits; env XSTEP=1 0.1-dex X grid, XEXT=1 X=112-178) -> out/q42/fit_<pos>_<var>_s<k>_lX<..>.npz; run scripts out/q42/run_all*.sh, logs run_all*.log. Variants: H4c/H4g/H4d (Harada + H13CN 4-3), T/Tp/T150/T32/T432 (common-taper ratios)
- scripts/s250c_q42_analysis.py -> out/q42/analysis.txt/.json (X posteriors, branch probabilities, BF, N(HCN)/N(HCO+), X_S/X_N)
- scripts/s250d_q42_branch_predictions.py -> out/q42/branch_predictions.txt (predicted alpha, tau, ratios per branch), branch_ratio_terms.txt (per-datum chi2)
- scripts/s250e_q42_common_taper_ratios.py -> out/q42/common_taper_ratios.txt (model-free ratios and contrasts)
- scripts/s251_update_report_q42.py -> /output/report.json (221 results: new Q42a, Q42b, Q42c; entries Q40 (216) and the H13CN opacity-bound entry (104) annotated with a Q42 note)

## Results
**A. H13CN 4-3 does not constrain anything in this model (answer to the first request).**
- With the censored limit, the fitted alpha, or the dilution variant, the N result is unchanged: P(Tk<=30 K) = 0.0026 / 0.0018 / (0.020 vs 0.016 for 20/25 % priors) vs 0.0021 without it (10/15 % priors, X prior LN(45, 0.2 dex)); X_N = 32-33 [27,42]; N(HCN)/N(HCO+) = 3.49 [3.15,3.85]; BF hot:cold 380-540 vs 457.
- Why (branch_predictions.txt): the RADEX alpha(4-3) = 0.016-0.035 in every Tk branch and every X = 32-63, also at tau(HCN 4-3) = 7-19 in the cold branch. Radiative trapping raises the 4-3 excitation of the thick HCN line while the thin H13CN is subthermal, so the Q4 equal-Tex relation ratio ~ tau/X fails; the limit 0.12 is 4-7x above any prediction, and a measurement would need ~+-0.005 to separate branches that differ by <40 %. The Q4 tau(HCN 4-3) bounds (entry 104) are now annotated as equal-Tex statements only.
- H13CN 1-0/HCN 1-0 is the discriminating H13CN datum: model 0.12-0.13 cold, 0.085-0.09 mid, 0.04-0.06 hot vs 0.0207 +- 0.0034 observed (2.2" aperture), 0.030 +- 0.006 (peak).
- Stated plainly: H13CN 3-2 (259.012 GHz rest, 256.6 GHz sky) and H13CO+ 3-2 (260.255 GHz rest, 257.8 GHz sky) are not covered by the Band 6 tuning (249.0-252.6 and 262.8-266.4 GHz sky). HCN 4-3 (spw16) and H13CN 4-3 (spw20) are different spws of one OUS (beams 4.14x2.70" and 4.26x2.77", 3 % difference, not corrected), not one cube. H13CN 4-3 sits on the CO 3-2 blue wing (Q4: wing-model-dependent, alpha centres 0-0.095).

**B. Common-taper (2.2") ratios and the X contrast (second request).** Table (R = I(13C)/I(12C), K km/s):

| | N | S | C |
|---|---|---|---|
| H13CO+/HCO+ template aperture | 0.0231 +- 0.0031 | 0.0184 +- 0.0040 | 0.0185 +- 0.0033 |
| +-150 aperture | 0.0271 +- 0.0049 | 0.0150 +- 0.0047 | 0.0182 +- 0.0045 |
| peak pixel | 0.0372 +- 0.0066 | 0.0141 +- 0.0088 | 0.0222 +- 0.0087 |
| H13CN/HCN template aperture | 0.0207 +- 0.0034 | 0.0057 +- 0.0059 | 0.0237 +- 0.0044 |
| H13CN/HCN peak pixel | 0.0302 +- 0.0058 | -0.0006 +- 0.0130 | 0.0265 +- 0.0117 |

Harada 1.7" peak: N 0.0458 +- 0.0049 (HCO+/H13CO+ 22 +- 3), C 0.0247 +- 0.0069. N falls 2x (12C/13C(thin) 22 -> 43) from the 1.7" peak to the 2.2" aperture.
- Model-free R_N - R_S: H13CO+ +0.9 / +1.8 / +2.1 sigma (template / +-150 / peak); H13CN +2.2 sigma in each. N vs C: 0-1.4 sigma (H13CO+), consistent with 0 in H13CN.
- RADEX (LN(45, 0.2 dex) prior, 10/15 % band priors): X_N = 58 [48,68] (aperture), 41 [34,53] (peak), 50 [43,59] (+Harada 3-2) vs 33 [27,41] with Harada intensities; X_S = 84 [70,95] (X<=100; 99 [76,119] to X<=178), 76 [60,93] peak; X_C = 52 [45,64]. **X_S/X_N = 1.3 [1.0,1.6] (aperture), 1.8 [1.2,2.2] (peak), 1.5 [1.2,1.8] (+3-2); P(X_S>X_N) = 0.89-0.95; with 20/25 % band priors 1.0-1.2 and P = 0.6-0.75.** The Q40 contrast (X_S 78-85 vs X_N 31-35, ~2.5) was mostly a beam/aperture artefact; at one taper it is <=1.8 and <=2 sigma. X_S is only a lower bound (>~65) because H13CN/HCN at S is not detected (<0.0177). H13CN/HCN (N vs S 2.2 sigma) is the stronger discriminator than H13CO+; opacity (HCN thicker at N) vs a real isotope gradient cannot be separated.

**C. Consequence for Q40 (cold/dense thick branch at N).** P(Tk<=30 K) (X<=100, LN(45,0.2 dex)): Harada intensities 0.002 / 0.016 / 0.025 (10/15, 20/25, 40/40 % band priors); common-taper aperture ratios 0.005-0.014 / 0.20-0.28; peak-pixel 0.04-0.06 / 0.47-0.51; aperture + Harada 3-2 0.02-0.07 / 0.46-0.58. So the thick branch is disfavoured (P<=0.07) only with tight band-scale priors, and NOT excluded with 20/25 % scales or aperture-averaged 13C/12C; Q40/p0044 must be reduced accordingly. The hot branch (Tk>=100 K) stays the best fit in every variant; N(HCN)/N(HCO+) = 3.4 [3.1,3.7] there (1.5-3 when the cold branch carries weight). Bridge: P(cold) <= 0.0007 (X_C 52). S: 0.01-0.11. In the hot best-fit cell the model overpredicts H13CN/HCN 1-0 by ~2x at X = 45 (tau(HCN 1-0) ~ 2), which pulls X_N to ~58.

## Literature used
- Harada+2018 Table 13 (via lit/harada_table13_parsed.csv): 1.7" peak-pixel H13CO+ 1-0/3-2 at N, C.
- Tang+2019 (12C/13C 30-67 in NGC 253 nuclei, 24-62 in NGC 1068, as used in Q40) and Martin+2019 (A&A 624, A125; NGC 253 12C/13C ~ 21 +- 6 from 13CO/C18O): X = 21-67 is the starburst-nucleus range; our N/C values 41-58 and S >65 are at/above its upper end. Jimenez-Donaire+2017 (MNRAS 466, 49) for the tau(13C)=tau(12C)/X relation (equal-Tex, which Q42a shows is not valid for the 4-3 line).

## Uncertain
- The 13C/12C ratios are 2.2" while the 12C lines in the fit are Harada 1.7" peak pixels; aperture vs peak moves X_N by 1.4x and P(cold) from 0.005 to 0.06 (tight priors). A fully matched set would need the Band 3 12C and 13C lines at one aperture AND the 3-2 12C at the same aperture (not done; Q34 has 12m+7m 3-2 at 2.2" for HCN/HCO+ but no H13CO+ 3-2 exists in the data).
- X_S is a lower bound; X > 100 is outside the prior used for most numbers (extended to 178: P(cold) at N rises to 0.15-0.56).
- H13CO+ uses the HCO+ LAMDA file (no H13CO+ file installed); H13CN uses h13cn@xpol (HCN-file convention would lower H13CN 4-3 by 1.15-1.75x, making it even less informative).
- One zone, 8-12 data points; the alpha wing model envelope (+-0.04) comes from Q4 (s41).
