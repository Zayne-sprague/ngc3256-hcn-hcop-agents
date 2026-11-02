# Work summary - Q5 round 2 (shock/chemistry: SO, CS 7-6, SiO, CH3OH; CS/HCN and SO/HCO+ ratios)

Responds to the reviewer's five fixes. All numbers are in /output/report.json (new entries prefixed `Q5r2:`; the round-1 `Q5:` entries are kept, the old N(CS)/N(HCN) abundance entry is replaced; interpretation + Q5 caveat rewritten). Scripts are re-runnable (`python3 -I <script>`); round-1 chain s90-s96 then round-2 chain below.

## Scripts (all /workspace/scripts) and outputs (/workspace/out/q5)
- `s97a_q5_smooth_1p7.py [Xb4s|Xb4t]`: Xb4d Band 6 cubes smoothed per channel (radio_beam kernels, flux x beam-area ratio) to 1.7" circular (Harada's resolution) or to 4.31x2.79" (Xb47 CS beam) -> spec_Xb4s_*.npz / spec_Xb4t_*.npz.
- `s97b_q5_b6_1p7.py`: HCO+/HCN 3-2 +-150 km/s moments in K km/s (Harada's definition), CH3OH 254.015 and SO 251.826 matched filters at 1.7" -> harada_compare_1p7.csv.
- `s97c_q5_harada_vs_ours.py`: Harada Table 13 vs ours -> harada_vs_ours.csv. Table 13 parsed from the arXiv PDF (lit/harada2018_1801.05941.pdf/.txt, lit/harada_table13_parsed.csv, 324 rows); N and S rows identical to team-a's /snapshots/team-a-s04/out/q4c/harada_table13_rows.csv.
- `s97d_q5_ch3oh_cs_matched.py`: CH3OH/HCO+ and CH3OH/CS in the matched 4.31x2.79" beam -> ch3oh_cs_matched_4p3x2p8.csv.
- `s98_q5_cs_matched_beam.py`: Xb47 CS/HCN/HCO+ 4-3 in a common beam (4.31x2.79" and 6.5"), offset-corrected, off-source exceedance null, N+S(+Bridge) point-source decomposition -> cs_matched_beam.json.
- `s99a_q5_cs_hcn_thick.py`: RADEX N(CS)/N(HCN) at tau(HCN 4-3)=0.3-10 -> cs_hcn_thick.csv.
- `s99b_q5_corrections_and_lte.py`: offset-corrected significances, deblended ratios, LTE extrapolation of Harada's SO -> q5r2_summary.json.
- `s100_update_report_q5r2.py`: report update (idempotent, run after s96).

## Results
1. **Harada+2018 comparison.** Harada has no Band 7 (no CS 7-6, SO 8-7) and no SiO/SO in our Band 6 lines, so only HCN/HCO+ 3-2 and CH3OH 254.015 (TNE/TSW only) overlap directly. At 1.7": HCO+ 3-2 N 86.6 vs 91.56 K km/s (0.95), S 54.2 vs 43.7 (1.24); HCN 3-2 N 64.7 vs 71.4 (0.91), S 20.2 vs 17.7 (1.14); peaks agree to 1-4%; HCN/HCO+ 0.747/0.373 vs 0.780/0.406 (S integral offset is a position / velocity-window effect). CH3OH 254.015/HCO+ 3-2 (K): N 0.0055+-0.0056 (<0.022), S 0.011+-0.006, **Bridge 0.041+-0.004**, Sblue 0.068+-0.029 (2.3 sigma), W 0.071+-0.014; Harada's 254 values are non-detections (limits ~0.02-0.05), so no tension. Harada's own CH3OH/HCO+ is enhanced at S (1.6-4.7x), OS (7-9x), C, TSW, SW; ours is 2-7x higher off the nuclei - same sense. SO: Harada N SO 219.949/HCO+ = 0.026+-0.008; LTE extrapolation to 251.826 gives 0.5-1.5 K km/s (10-40 K) vs ours 0.39+-0.28: consistent for Tex<~25 K; our SO 8-7 limit constrains Tex>30 K only if SO fills the 4.3" beam. No direct SiO test (not covered).
2. **CS N/S reconciliation with team-a p0009.** Definition/aperture difference, not data conflict: team-a's 0.051/0.056 are CS/HCO+ *3-2* deblended, or the 6.5" aperture (we reproduce 0.07/0.06 at 6.5"); their native-beam 4-3 line table (0.08/0.045) agrees with ours (0.087/0.043). In a common 4.31x2.79" beam with N+S deconvolution: CS/HCN 4-3 N 0.165+-0.014, S 0.182+-0.028 (N/S 0.91+-0.16); CS/HCO+ 4-3 N 0.088+-0.008, S 0.063+-0.009 (N/S 1.4+-0.24). The old "N/S 1.20+-0.10" CS/HCN is withdrawn in favour of 0.9-1.2+-0.2 (uniform).
3. **Baseline offset.** Off-source median of the CS matched amplitude -0.13 Jy/beam km/s (N, S, Bridge, Sblue), -0.05 (E, W). Corrected: N 14.8, **S 6.5**, Bridge 9.2, Sblue 3.6 sigma (tentative), E1/E2/W <2. SO 6-5 Xb4f S 4.2 -> 4.6 (SO-only), 2.9 -> 3.3 (with CH3OH ladder): still tentative. A reflection null is meaningless (pointing = N-S midpoint) and the spw (+-800 km/s) is too short for line-free on-source windows, so the null is the off-source exceedance (0.00/0.015/0.003 for N/S/Bridge at 4.3x2.8"). Caveat: S CS is only 2.7 sigma in a single-pixel box moment.
4. **Abundance.** N(CS)/N(HCN) dropped as a result. Robust: CS 7-6/HCN 4-3 equal within 20-25% at N, S, Bridge. RADEX at tau(HCN 4-3)=0.3-10 and ladder-allowed (Tk>=56 K, n 1e4-1e5): 0.25-0.43 (full grid 0.17-0.67, factor 3 from (Tk,n)); H13CN limit does not tighten it (needs tau>=2-6 for 12C/13C 40-90).
5. **Net answer** (explicit in report): no shock signature in the nuclei (SO 8-7/HCO+ <0.011-0.014, CS/HCN uniform, SO 6-5 tentative and non-specific, SiO not covered); the shock tracer that is detected is cold, extended CH3OH, higher off the nuclei (Bridge, Sblue, W), consistent with Harada and with the low Trot of Q12.

## Remaining uncertainty
- CS 7-6 blend with H2CS/29SiO (S, Bridge unconstrained by data); on-source baseline offset only inferred from the field (+-0.13); S CS significance depends on estimator (2.7-8 sigma).
- Sblue = Harada's OS is an identification (their coordinates not tabulated); Sblue CH3OH is 2.3 sigma.
- Harada comparison uses different transitions for SO/CH3OH (LTE, Tex unknown, factor 2-5); S position/velocity window mismatch (14-24%).
- Deblended N/S ratios are point-source-model dependent (~20%); the 4.3x2.8" beam cannot separate N and S cleanly (5.2" apart).
- Notation note: 251.826 GHz is SO J_N 5_6-4_5 (Eu 50.7 K); the proposal's/Harada's 6_5-5_4 is 219.949 GHz (not covered).
