# Work summary: Q38 round 2 (reviewer fixes)

## Done
1. **Report self-consistency (fix 1).** `scripts/s245_update_report_q38r2.py` rewrote (not just appended) the disk-annuli entry (99), matched-J 10" (96), observed 10" ladders (111), 10" region comparison (112), disk-contrast (133), disk-positions (136) and the Q38 entry (217). No entry now states a "faster disk decline" or a disk-vs-nucleus excitation contrast. Backup: `report_backup_q38r2.json`.
   RADEX with the uv-complete matched 2.2" data only (`s243_q38_radex_E1.py` -> `out/q38/s243.*`; 1-0 X961 12m + 3-2 Q34 12m+7m, 4 data/position): N Tk 40 K [14,158]; E1 178 K [79,251]; **P(Tk_E1<Tk_N)=0.17** (10" blended value was 0.91-0.99, withdrawn); log A E1 0.30 vs N 0.70 (P(A_E1<A_N)=1.00); BF same-Tk&n +1.0, same-A -2.1. E1 differs from N by HCN/HCO+ abundance/opacity (already R10 0.49 vs 0.84), not by lower Tk. Caveat: only 4 data, Tk loosely constrained.
2. **Clump vs nucleus (fix 2).** `s244_q38_clump_vs_nucleus.py` -> `out/q38/s244.*`. Same NNLS method: nuclear complex R21/R32/R43 = 0.58/0.45/0.38, R32/R21 0.78+-0.12, R43/R21 0.66+-0.07 (N alone 0.88+-0.15). E1 0.69+-0.17, W 0.68+-0.15. Double ratio clump/nucleus R32/R21: E1 0.88+-0.25, W 0.87+-0.23 (vs N: 0.78+-0.23, 0.77+-0.21): 0.5-1.1 sigma from 1, i.e. a COMMON decline. 4-3: E1 R43/R21<0.17 vs 0.66 nucleus: suggestive only (PB 0.5, zero level).
3. **NNLS injection (fix 3).** `s242_q38_injection.py` -> `out/q38/s242.*` (s58f simulator: ACA-like uv 8.7-45 m, Hogbom, PB; truth flat or declining ladders; optional unmodelled 25/40" extended halo). E1/W recovery: R21 1.00-1.01, R32 1.01-1.03, R43 0.92-0.93 of truth; R32/R21 recovered 1.00-1.03 (flat), 0.70->0.70-0.73. Nuclear ladder exact. E2 biased (R43 0.2-0.5 of truth). Limits: same morphology in both lines, no noise, assumed (not real) uv coverage; real-uv tapered imaging (common uvrange for Bands 5/6) was NOT done.
4. **E2 (fix 4).** E2 dropped from ladder/double-ratio claims: 12m-only vs comb 3-2 differ 2.5x (HCO+ 1.9 vs 4.7 K km/s), NNLS biased at E2. E1 and W are insensitive to 12m-only vs comb (R32/R10 0.72/0.74, 0.89/0.89). The 1-0 12m recovery itself is unquantified.

## Verdict on p0032/p0033
p0032 CONFIRMED for 8-12" 3-2/4-3 annuli (not measurable) and our 0.52/0.34/0.14 were clump positions at 10": REVISED. The clump 2-1->3-2 decline (~30%) is real but common to nuclei and clumps; "disk declines faster" and "disk has lower Tk" are not supported. 4-3 clump: limit. p0033's small 1-0 disk decline is compatible: the real clump-vs-nucleus difference is the 1-0 ratio itself (0.49 vs 0.84).

## Remaining uncertain
Real-uv tests (matched uvrange tapered Band 5/6 images), 12m 1-0 recovery at E2/W, Xb47 zero level/PB for 4-3 clump limit, other entries that use 10" clump values (e.g. entries 112, 131, 136 only annotated).
Literature: none new (Harada+2018 as already cited in the report).
