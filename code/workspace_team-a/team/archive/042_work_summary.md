# Q38 work summary
Question: is the Q36 10" continuum Tb (5-8 mK) a lower limit because of missing short spacings, and does a FIR-SED dust disk bound the envelope background above CMB?
Scripts: /workspace/scripts/q38_mrs.py (MRS + global-dust prediction), q38_filter.py (uv-hole filtering of a Gaussian disk), q38_jbg.py (J_bg/Tbg_eq of envelope for disk F, Td, beta, eta; reuses /snapshots/team-a-s22/scripts/q36_bg.py). Outputs: /workspace/out/q38/*.log.
Results:
- MRS (inferred from beam sizes/Cycle-5 handbook, not from visibilities): Band 6 ~6-11", Band 7 ~5-8", Band 5 ~10-17". Any smooth dust disk with FWHM >~10-15" is fully resolved out (q38_filter), so 5-8 mK is a lower limit.
- Global dust from Michiyama+2020 FIR fluxes (160 um 93.5 Jy), optically thin: S258 = 140-690 mJy (255 mJy for Td 30 K, beta 1.8) vs ~35-70 mJy recovered (half dust).
- Intrinsic central Tb_RJ(258) of that disk: 47/21/12/5 mK (F=10/15/20/30"), up to 85-127 mK for beta=1.5 or Td=25 K at F=10". <=10-35% of J_CMB(258) for F>=15".
- Envelope Tbg_eq (HCN/HCO+ 2-1/3-2/4-3): 2.73-2.77 / 2.74-2.90 / 2.75-3.24 K (F>=15", Td>=30 K); extreme case 2.89/3.23/3.89 K. All inside the q36 'envhi' scenario already run in Q22 (P(Te>=50)=0.49 nominal / 0.37 high cores vs 0.43 CMB-only), so no new Q22 run was needed. Reaching >=5 K would need all FIR dust within <~5", excluded by the resolved cores and the 10" flux.
Conclusion: Q36 stands: envelope Tbg ~ 2.73-3.2 K (<=3.9 K at 4-3), hot envelope neither required nor excluded.
Literature: Michiyama+2020 (FIR fluxes), Sakamoto+2014, Q36 notes.
Uncertain: actual MRS (needs visibilities/listobs: could be requested), disk size (a Herschel 160/250 um map or SCUBA/LABOCA/APEX single-dish 870 um would fix F), beta/Td extrapolation (factor ~4), eta geometry, Gaussian-disk idealisation.
