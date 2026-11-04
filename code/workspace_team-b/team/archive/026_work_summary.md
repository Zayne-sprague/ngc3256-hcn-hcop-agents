# Q56 work summary: matched-beam HCN R32/10, reconciliation

**Answer.** The earlier statement (HCN 3-2/1-0 equal at N and S) is right to ~10%. The Q53 values (HCN N 0.306, S 0.217, S/N 0.71) are an artefact of the uv-incomplete 3-2 image (Xb4d 12m only). It loses 25-32% of the HCN flux and 14-21% of the HCO+ flux at S, but only 0-5% at N.

**Method.**
- Restored the 12m+7m visibilities (s341) and re-imaged 3-2 at circular 2.5" and 3.5" with a taper (out/q34/*_2p5_comb, *_3p5_comb).
- s342 measures R32/10 in RJ Tb at the same pixels, windows and 1-0 data (X961) as Q53 (s332/s340).
- Scripts: s340 (Xb4d image-plane sweep, 2.5-3.5"; beams <2.5" are below the native Xb4d beam), s341 (CASA), s342 (comb measurement), s343 (report update).
- Outputs: out/q56/s340.txt, s342.txt, s342.json.

**Results (|v|<100; |v|<150 and <350 in brackets).**

| Beam | Line | N R32/10 | S R32/10 | S/N | S/N Xb4d image |
|---|---|---|---|---|---|
| 2.5" | HCN | 0.304 | 0.271 | 0.89 (0.94; 0.96) | 0.71 |
| 2.5" | HCO+ | 0.396 | 0.415 | 1.05 | 0.89 |
| 3.5" | HCN | 0.283 | 0.246 | 0.87 (0.90; 0.92) | 0.69 |
| 3.5" | HCO+ | 0.372 | 0.368 | 0.99 | 0.84 |

- The earlier Q34 aperture values at 2.2" (all velocities) were HCN N 0.263 / S 0.261 and HCO+ 0.351 / 0.343, so they agree.
- The residual HCN S/N of 0.87-0.96 is within 1-2 sigma of the window/pixel systematics (about ±0.05). At most 5-13% lower HCN excitation relative to HCO+ at S is hinted.
- N/S HCN/HCO+ contrast at 2.5": 1-0 is 1.48-1.56; 3-2 uv-complete is 1.68-1.74 (Xb4d image 1.78-1.85). D = 1.08-1.18 by window (Q34 aperture: ~1.0).
- Excitation alone can therefore account for at most ~10% of the contrast. The remaining ~90% needs abundance, opacity or column differences. A full RADEX bracket is left to Q34/Q40 (s193, s243).
- Checks: the HCN 1-0 hyperfine satellites (-7 and +5 km/s) lie inside the windows. The 3-2 lines are fully in the comb cube (-600..+550 km/s).

**Report.**
- Revised the entries for Q53 R32/10 and the Q53 verdict, and annotated Q35, with the superseded claim flagged.
- Added two Q56 entries.
- Backup: report_backup_q56.json.

**Uncertain.**
- The 2.5" and 3.5" beams are image-plane convolved for the 1-0 data but taper-restored for the 3-2 data.
- The pixel versus aperture difference moves S/N by ~0.05.
- Only the comb cubes were imaged for 2.5" and 3.5" (no 12m-only variant).
- Beams below 2.5" are not accessible from Xb4d image-plane data. Q34 covers 1.5" and 2.2".
- Literature: nothing published on this quantity for NGC 3256 (Harada+2018 give 1-0 ratios at 1.7").
