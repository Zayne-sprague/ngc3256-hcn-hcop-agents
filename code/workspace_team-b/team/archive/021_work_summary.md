# Q53 work summary (worker, team-b)
Question: is the S deficit of HCN/HNC/CN relative to HCO+ a depletion of N-species or an HCO+ boost from outflow/shock gas?

## Done
- scripts/s330_q53_velbins.py (reuses s320 loader/baseline; velocity bins core |v|<100, blue/red 100-350, wings, core50; 1-0 at 2.5/3.5", 3-2 at 1.4/2.5") -> out/q53/s330.txt/.json
- scripts/s331_q53_13co_ref.py: 13CO(2-1) (2015.1.00412 X155 spw0, 0.85") as column reference -> out/q53/s331.txt. (C18O 2-1 and CO 1-0 isotopologues are not in the cached archive files; C18O 2-1 falls just outside X155's band.)
- scripts/s332_update_report_q53.py: 4 entries added to /output/report.json (backup report_backup_q53.json).

## Results (N/S of X/HCO+, 2.5" J=1-0 unless noted)
- Core |v|<100: HCN 1.35+-0.01, HNC 1.78+-0.05, CN(G2) 1.88+-0.02 (all-velocity 1.56, 1.90, 2.07). 3-2 1.4": HCN 2.10 (all 2.14), HNC 1.91 (1.84).
- HCO+ wings (100-350) = 7.6% of 1-0 flux at N, 21% at S (3-2: 3%, 15%): real S outflow/shock contribution, but removing it lowers the contrast only 13% (HCN, 1-0), 6% (HNC), 9% (CN), <=2% (3-2). HCO+ boost explains <=15%.
- Wing ratios: HCN/HCO+ N 1.72+-0.14 vs S 0.43+-0.03 (N/S 4.1, 1-0); S wings are also N-poor, N wings N-rich. N wings are weak (HCO+ 3-2 wing 0.7-2 K km/s): errors large at 3-2.
- 13CO(2-1): N 60.3, S 57.8 K km/s at 2.5" (N/S 1.04). X/13CO N/S: HCO+ 1.36, HCN 2.12, HNC 2.58, CN 2.80 (1-0, 2.5"); 3-2 at 1.4": 1.66, 3.55, 3.05. Same pattern as vs 1.1 mm continuum, so not tied to continuum contamination at S.
## Conclusion
Predominantly a deficit of HCN/HNC/CN at S (or enhancement at N), present in the line core; outflow HCO+ boost contributes <=15% at 1-0, ~0 at 3-2.
## Uncertain
Abundance vs excitation/opacity within N-bearing group not separated (needs H13CN/H13CO+ matched beam, CN hyperfine RADEX, joint multi-ladder fit); 13CO(2-1) opacity and excitation; no uv matching; wing baselines at N are noise-limited. No new literature used beyond Q52 (Harada+2018 Table 13, Michiyama+2018 outflow, quoted from earlier notes, not re-read).
