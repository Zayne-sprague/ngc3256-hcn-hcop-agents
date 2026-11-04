# Q40 (round 2) work summary

## Bottom line
With the matched-beam Harada 1.7" H13CO+ 1-0 and 3-2 at N added to the Q39 one-zone RADEX fit, the cold/dense optically thick branch (Tk<=30 K, log n>=5.4) is **disfavoured, not excluded**. The statement is conditional on 12C/13C <~60 and on the one-zone model. The reviewer was right that the round-1 "excluded at >=4 sigma" and "hot branch 170-300 K" wordings were not supported; they are withdrawn.

## What was done (scripts in /workspace/scripts, outputs in /workspace/out/q40r2)
- `s235_q40r2_Xmarg.py`: re-runs the s230 machinery (imports s220b/s230) at 14 values of X=12C/13C (log X 1.35-2.0, step 0.05; X=22-100) for three band-scale prior widths (10%/15% = Q39 default, 20%/25%, 40%/40%). Variant H = Harada H13CO+ only (all lines 1.7" peak pixel, same beam). Variant A2 = H + H13CN 1-0 and S/C H13CO+ from the Q29 2.2" maps, given an extra 30% beam-mismatch systematic. Run for N (H, A2) and S (A2). About 18-50 s per fit.
- `s236_q40r2_analysis.py`: (A) 13C ratio alone; (B) per-datum decomposition of chi2 at the best cell of each Tk branch; (C) marginalisation over X, branch probabilities, Tk-cap sensitivity. Outputs `analysis_N_H.txt`, `analysis_N_A2.txt`, `analysis_S_A2.txt` (the `analysis_*.txt` files in out/q40r2 may hold the last run of each).
- `s237_q40r2_ratio_marg.py`: N(HCN)/N(HCO+) marginalised over X, N/S -> `ratio_marg.txt`. `s238_update_report_q40r2.py` writes /output/report.json (backup of the previous version: /workspace/report_backup_q40r1.json).
- Bug found and fixed: the round-1 `lnZ` (out/q40/summary.txt) used exp(-0.5 chibest) but Wc is built from three species each shifted by chibest, so Z = Wc exp(-1.5 chibest). It is harmless at fixed X (round-1 branch probabilities stand) but wrong when comparing X. The round-2 numbers use the right factor.

## Results
1. **Decomposition (reviewer fix 1)**, N, variant H, band priors 10/15%, best cell per branch (cold Tk<=30 vs hot Tk>=100):

   | X | total dchi2 | band/eta prior | 12C HCN | 12C HCO+ | 12C CS | 13C pair |
   |---|---|---|---|---|---|---|
   | 40 | 17.2 | 3.1 | 1.0 | 6.6 | 4.3 | 2.3 |
   | 63 | 4.2 | 0 (9.3 in both) | -0.1 | 1.4 | 1.2 | 1.7 |

   So at X=40 only ~2 of the 17 comes from the 13C lines themselves. The rest is the cold branch's inability to fit the 12C lines and the 13C columns together. The per-datum pulls are in `analysis_N_H.txt` (the 13C 3-2 is over-predicted by the cold branch: 3.9 vs 2.6+-0.49 at X=40, 3.6 at X=63; the hot branch gives 2.5 and 2.3). The band-scale penalty does matter: at X>=63 both branches pay 9-13 in prior, so these fits are strained by the 12C line ratios regardless of branch. Widening the band scales to 20/25% or 40/40% removes the penalty and the best chi2 drops to 4.6-5.1 or 2.3-3.3.
   **13C ratio alone** (R(3-2/1-0)=0.277+-0.06 stat, 30% total with band scales): dchi2 = 0 between branches. Some (n, N13) reproduces the ratio at every Tk, so the ratio by itself discriminates nothing. The discrimination comes from tying N(H13CO+)=N(HCO+)/X.
2. **X as a nuisance (fix 2)**: log-normal priors centred at X=45 or 60 with 0.2 dex (this brackets the starburst-nucleus values of Tang+19 A&A 629 A6: NGC 253 30-67, mean 42; NGC 1068 24-62; NGC 4945 6-44) or flat in log X. The data themselves pick X: lnZ(X) peaks at 30-50 and falls by 5 by X=100. Posterior median X_N = 31-35 (10/15% scales), 40-51 (20-40%). 13CO/C18O=2.4 at N (Harada Table 13; 5.9 at S, 3.3 at C) does not constrain X. A low value at N points to enhanced 18O (cf. Sliwa+17, Brown+19 for ULIRGs) rather than low 12C/13C, because low 12C/13C would raise 13CO/C18O.
   Marginal posterior of Tk branches at N (prior mass of each is ~0.32/0.35/0.32 for <=30 / 30-100 / >=100 K; the grid is log-uniform in Tk):
   | band-scale prior | P(Tk<=30) | P(30-100) | P(>=100) | BF(>=100 : <=30) |
   |---|---|---|---|---|
   | 10/15% | 0.002-0.006 | 0.02-0.04 | 0.95-0.98 | 170-1200 |
   | 20/25% | 0.007-0.032 | 0.11-0.19 | 0.77-0.88 | 24-135 |
   | 40/40% | 0.011-0.047 | 0.16-0.25 | 0.70-0.83 | 16-73 |
   The ranges span the X priors above. These are Bayes factors on a 3-4 parameter model with 8-9 data and should be read as "moderate to strong, model-conditional".
3. **Beam (fix 3)**: the primary result uses only the Harada 1.7" pair for the 13C data. The A2 variant (Q29 H13CN 1-0 added at 2.2", extra 30% systematic) gives P(Tk<=30)=0.000-0.018, i.e. the same result. S with A2 is acceptable (chi2 2-8) but wants a higher X (median 65-85) than N. That is either a real gradient (Tang+19 find the lowest 12C/13C at nuclei and the highest in outskirts) or the 2.2"/1.7" aperture mismatch, and I cannot separate them without a common-taper re-extraction (not done; raw cubes at a common 1.7" beam would settle it). The round-1 variant-A "S fits badly" is a fixed-X artefact: S is fine once X_S is free. One shared eta per gas for HCO+ 1-0 and 3-2 is retained.
4. **Tk cap and intermediate branch (fix 4)**: the hot solution sits at the 300 K cap (log n 4.6-4.8), so only a **lower limit** on Tk is established (Tk >~ 50 K). Best-over-X chi2 cost of capping Tk: <=150 K: 5.6 / 1.0 / 0.4 (10/15, 20/25, 40/40 priors); <=100 K: 8.2 / 1.5 / 0.5; <=50 K: 10.5 / 2.7 / 2.0; <=30 K: 11.8 / 4.2 / 3.3. With loose band scales intermediate (50-150 K) solutions are allowed. The grid was not extended above 300 K (the HCN/HCO+ LAMDA files are fine to that range; no time for a new RADEX grid).
5. **N(HCN)/N(HCO+) and N/S** (marginal over X, Tk, n): N = 3.5 [3.15,3.8] (10/15%) to 3.7-3.85 [3.2,4.7] (loose); S = 2.2 [1.8,2.6]; N/S = 1.6-1.8 [1.3,2.4], P(N/S>1)=0.98-0.99. The profile significance stays ~2.1 sigma, limited by 12C line S/N. N(HCN)/N(HCO+) varies little with the branch assumption, which is why the abundance contrast is more robust than the temperature statement. See `ratio_marg.txt`.

## Literature
Harada+2018 Table 13 (12C and H13CO+ intensities at 1.7"); Tang+19, A&A 629, A6 (12C/13C 6-67 in starburst nuclei from 12CN/13CN; prior); Sliwa+17 (ApJ 840 L11) and Brown & Wilson 2019 (extreme 13CO/C18O in ULIRGs: 18O enhancement) for the 13CO/C18O interpretation. Henkel+14 (NGC 253) was not read in this round.

## What remains uncertain
- One-zone model; the thick N core and envelope are not separated (a two-gas fit with the 13C data is the natural next step). The 3-2/1-0 relative calibration enters through the band-scale prior and drives the strain at high X.
- Tk above ~100 K is not constrained (grid cap 300 K).
- X_N and X_S are inferred, not measured: an independent measure (13CN/CN, 13CO/C18O with a 16O/18O prior, H13CN/HC15N) is needed. A common-taper re-extraction of H13CN/H13CO+ at S and C would remove the beam caveat.
