# Q40 work summary (round 2): clump vs nucleus HCN/HCO+ (reply to team-b p0045, p0032-r01, p0033-r01)

**Reply to team-b: CONFIRMED (a clump deficit exists; the 2-1->3-2 decline is common to nucleus and clumps). REVISED ('mostly in 1-0' is not supported: the E1 deficit stays or deepens with J, W is flat; vs S alone the 1-0 deficit is small). 4-3 and E2: limits only.**

## What changed after the review
1. **Harada vs Michiyama resolved (uv coverage, not astrophysics).** `scripts/q40c_uvcmp.py` (out/q40b/q40c_uv.log): beam-deconvolved FT of both 1-0 maps on a common grid. Harada (beam 1.38x1.14") has H/M visibility amplitude 0.3-0.4 below 4 kl, 0.6-0.75 at 4-12 kl, ~0.8 above (both lines; the small common offset cancels in R). Plain 10" HCN: Harada lower by 5% N, 9% S, 33% E1, 20% W. `scripts/q40d_uvfilter.py` forward-models a common uv high-pass (u>umin, applied to the maps and to the component templates) at one common 8.64x6.0" beam; umin=0 reproduces the earlier native decomposition (control). At u>8-10 kl Harada and Michiyama agree: nuclei within 1-3%, E1 R 0.39+-0.10 vs 0.50+-0.02 (1.1 sigma), W 0.54+-0.05 vs 0.49+-0.02. The old 'E1 R10 = 0.44 (H)' was spatial filtering; Michiyama (shorter baselines) is the 1-0 map to use. Restoring raw visibilities was not needed (Harada EB 58 GB, ~4 h); the image-plane forward model reproduces the uv dependence. The matched range u>8-10 kl is set by the 3-2 7-m shortest baseline (~8 kl).
2. **Common beam, position, size, injection** (`scripts/q40b_fit.py`, `out/q40b/q40b_fit.log`): E1/(N+S) at 1-0 (M) / 1-0 (H) / 2-1 / 3-2: native 0.79/0.62/0.64/0.49; common 8.64x6.0 0.77/0.51/0.64/0.32; circular 10" 0.76/0.46/0.65/0.24; matched uv (u>8 kl) 0.70/0.51/0.61/0.43, (u>10) 0.74/0.60/0.59/0.51. +-1" shifts change E1 R(3-2) 0.13-0.17 (central 0.155); 1-0, 2-1 <=0.04. Clump size 0-8": 1-0, 2-1 contrast changes <=0.02; 3-2 0.32->0.28. Nuclear envelope size 3" lowers the residual rms 10-30% and E1/(N+S)(3-2) to 0.24: **E1 3-2 is the least robust number** (0.24-0.5), but it is below the 1-0 value in every variant. Injection-recovery with real residuals: bias <=0.04 (E1), <=0.06 (W) for sizes 3-8".
3. **Single contrast table** (`scripts/q40e_sys.py`, `q40f_table.py` -> out/q40b/q40f_table.md/json). Central = mean of matched-uv u>8 and 10 kl; noise = filtered null windows added and refitted; sys = RMS over 14 variants x umin 6/8/10; same-J contrasts carry no calibration term; 10% only on absolute cross-band double ratios.

| R=HCN/HCO+ (T_B) | 1-0 Michiyama | 1-0 Harada | 2-1 | 3-2 |
|---|---|---|---|---|
| N | 0.82 | 0.84 | 0.71 | 0.58 |
| S | 0.58 | 0.59 | 0.53 | 0.42 |
| N+S | 0.70 | 0.71 | 0.62 | 0.50 |
| E1 | 0.50+-0.02 | 0.39+-0.11 | 0.37+-0.03 | 0.24+-0.05 |
| W | 0.49+-0.02 | 0.54+-0.06 | 0.47+-0.02 | 0.31+-0.03 |

| contrast (noise\|sys) | 1-0 M | 1-0 H | 2-1 | 3-2 |
|---|---|---|---|---|
| E1/N | 0.62 (.01\|.02)* | 0.47 (.12\|.06)* | 0.52 (.03\|.03)* | 0.41 (.05\|.07)* |
| E1/S | 0.86 (.02\|.03)* | 0.67 (.17\|.08) | 0.70 (.04\|.03)* | 0.56 (.07\|.10)* |
| E1/(N+S) | 0.72 (.01\|.02)* | 0.56 (.14\|.07) | 0.60 (.03\|.03)* | 0.47 (.06\|.08)* |
| W/N | 0.60 (.03\|.02)* | 0.64 (.06\|.04)* | 0.66 (.03\|.01)* | 0.54 (.02\|.05)* |
| W/S | 0.84 (.04\|.03)* | 0.92 (.08\|.05) | 0.88 (.04\|.01)* | 0.74 (.02\|.07)* |
| W/(N+S) | 0.70 (.03\|.02)* | 0.76 (.07\|.05) | 0.75 (.03\|.01)* | 0.63 (.02\|.06)* |

(* = >3 sigma from 1, noise and sys combined.) Plain 10" position values (E1/N 0.89 M, 0.74 H, 0.84, 0.65) and the native-beam decomposition (E1/(N+S) 0.79, 0.62, 0.64, 0.49; W 0.69, 0.76, 0.74, 0.75) are in the same table file; they show the same trend but the plain maps dilute the contrast.

Double ratios: R21/R10 N+S 0.89+-0.09, E1 0.74+-0.10, W 0.96+-0.11; R32/R21 N+S 0.80+-0.08 (N 0.81, S 0.79), E1 0.63+-0.16, W 0.67+-0.09; R32/R10 N+S 0.72+-0.07, E1 0.47+-0.10, W 0.64+-0.09. Clump/(N+S), calibration-free: E1 R21/R10 0.83+-0.07 (-2.3 sigma), R32/R21 0.79+-0.19 (-1.1), R32/R10 0.65+-0.13 (-2.7); W 1.07+-0.07, 0.83+-0.08 (-2.0), 0.89+-0.09 (-1.4). Nothing exceeds 3 sigma.

**Conclusion.** (i) R(clump) < R(nucleus) at every J, >3 sigma for E1 and W against N and N+S. (ii) The J dependence is clump specific: E1 contrast 0.72 -> 0.60 -> 0.47 (deepens, -2.7 sigma 1-0 to 3-2, not significant at 3 sigma and robust only qualitatively); W flat (0.70/0.75/0.63). So 'mostly in 1-0' is rejected. (iii) The 2-1->3-2 decline (0.8 nucleus, 0.63-0.67 clumps) is common to within errors: team-b's point is confirmed. (iv) Against S alone E1 is close to normal at 1-0 (0.86+-0.04) and W 16% low: the 'deficit' is largely relative to the N nucleus (R10 0.82 vs S 0.58). Literature: Harada+2018 (ApJ 855, 49, Sec. 3.2.1) also find the S nucleus has the lower 1-0 HCN/HCO+; no published clump values at 2-1/3-2 to compare.

## 4-3 and E2 (`scripts/q40g_43.py`, out/q40b/q40g_43.log)
HCN 4-3 at E1 0.00+-0.17 Jy km/s (empirical sigma at r=7.8"), 3-sigma limit 0.50; HCO+ 4-3 4.65+-0.19 -> R43(E1) < 0.13-0.17 (T_B, model variants), versus R32(E1) 0.24+-0.05: R43/R32 < 0.55-0.7 against 0.74 (N), 0.62 (S): a limit, no evidence of a faster decline. W: HCO+ 4-3 0.44+-0.19, not measurable. E2: R(1-0) 0.25+-0.03 (+-0.07 sys; Michiyama) vs ~0 (Harada), 3-2 HCN 0; unconstrained, not used.

## Physical meaning (RADEX, `scripts/q40h_radex.py`, `q40i_radexfit.py`, out/q40b/q40i_radexfit.log)
Grid T 20-150 K, log n 3-7.5, log N(HCN) 13-17.5 (dv 150), log X = N(HCN)/N(HCO+) -1..+1; chi2 on R10, R21, R32 (E1 also R43<0.15). Individually every region is unconstrained (3 data, 4 parameters; thin and thick models both fit). Sharing all parameters between nucleus and clump is excluded (chi2 115 N-E1, 76 (N+S)-E1, 49 (N+S)-W for 6 data), but changing ONE parameter is enough: X chi2 1.0-3.8, n 0.5-2.9, N/dv 0.7-3.0, T 1.0-11 (T worst vs N). Required clump minus (N+S) offsets: dlogX -1.5..-0.25 (E1; W -0.5..-0.25), dlog n -2..-0.5 (median -0.6, both), dlog N -4.5..-0.5. So the 0.5-0.8 factor can be abundance (HCN/HCO+ lower by 2-4), 0.5-1 dex lower density, or lower column/opacity; the data do not separate them (delta chi2 < 3). E1's deepening with J favours density/excitation or opacity; W (flat) favours a J-independent effect. H13CN: only N, S, Bridge have limits (HCN/H13CN(4-3) > 8.3 N, > 4 S; tau(HCN 4-3) <= 3.2/5.5/10 for 12C/13C = 30/50/90 at N); there is no limit at the clumps (HCN 4-3 itself undetected at E1/W), so clump opacity is not bounded by present data. What would: H13CN or HC15N (4-3/3-2) stacked at E1/W with a Band 7 or 6 pointing centred on them, or 3-2 HCN/HCO+ at 5.7" with a uv-matched 7-m+12-m 1-0 (Michiyama-quality, with ACA) map.

## Uncertain
- E1 3-2 R and its contrast (0.24-0.5 across nuclear-size/beam/uv variants). Positions fixed from team-b (best-fit shifts of E1 (0,-0.5") and W (+1,+1") improve the residual by 2% and 11%, no effect beyond the quoted variants).
- Matched-uv from images (not visibilities): clean maps contain model power below the shortest baseline; the u>8-10 kl cut was chosen to avoid it. Harada E1 is noisy (+-0.1).
- Each ladder uses different dishes/arrays; calibration 5% (R) and 10% (double ratios) assumed.
- RADEX one-zone, equal beam filling; ratio information is not enough to give n, T, N, X (see above).

## Files
scripts: q40b_fit.py, q40c_uvcmp.py, q40d_uvfilter.py, q40e_sys.py, q40f_table.py, q40g_43.py, q40h_radex.py, q40i_radexfit.py, q40j_report.py (earlier q40_clump_nucleus.py, q40_decontam.py); outputs in out/q40b/ and out/q40/. Report: /output/report.json (218 results; backup out/report.backup.q40r2pre.json).
