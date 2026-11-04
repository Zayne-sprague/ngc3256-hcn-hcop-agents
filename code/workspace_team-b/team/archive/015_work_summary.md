# Q51 work summary: is N CON-like?

**Answer: no.** N (and S) show no HCN-vib and have none of the CON continuum properties.

## Scripts / outputs
- scripts/s310_q51_con_like.py -> out/q51/s310.txt, out/q51/falstad_table3.csv, out/q51/nishimura_split.csv
- scripts/s311_update_report_q51.py (4 entries added to /output/report.json; backup report_backup_q51.json)
- reused: scripts/s42_hcn_vib.py -> out/s42_hcn_vib.csv (3 sigma HCN v2=1f limits, Q4)

## Results
1. HCN-vib: Falstad+2021 (Paper I) CON criterion Sigma_vib(3-2) > 1 Lsun/pc2. Our 3 sigma limits at N: L_vib < 29 Lsun in the 48 pc beam, Sigma < 0.017 (>=60x below threshold); <52 Lsun at 180 pc; <95 Lsun at 1.1 kpc; 4-3 <207 Lsun at 800 pc. S: <25, <43, <102, <310 Lsun. Falstad's own Table 3 gives NGC 3256 N, S L_vib < 65 Lsun (Sigma <0.001, <0.004), while the detected nuclei have 1.2-62 x10^3 Lsun.
2. vib/HCN L' ratio in detections (Falstad L_vib / Nishimura HCN L'): CONs 0.41-0.69 (beam), 0.06-0.18 (10 beams); ESO 173 0.66/0.058, F14378 0.10/0.022. Ours N: <0.025-0.05 (55 pc), <0.02 (250 pc), <0.01 (1.1 kpc); S <0.085-0.10 / <0.03-0.05 / <0.015-0.02. N is >=8x below all CONs at 55 pc, >=3x at kpc scale.
3. Continuum (Falstad Table 3, 1.1 mm): N 17.6 mJy, 347x250 pc, T_b 0.13 K; S 9.3 mJy, 242x106 pc, 0.24 K. CONs: 9.8-79 K (median 20 K), 8-60 pc. Non-CON median 0.25 K. Our 0.24" peaks 0.18 K (N), 0.59 K (S). Dust gas surface density bound (Q6) N <=1.5e3 Msun/pc2 vs 2.3e4 for N_H2=1e24: >=15x below.
4. Nishimura HCN/HCO+ split by Paper-I class (flags from Falstad Table 3 criterion, not the garbled col. 8): CON (n=3) resolved median 1.74 (1.34-2.34); vib-detected non-CON (n=2) 1.13-1.32; non-detected (n=11+1 limit) 1.03; 6 galaxies not in Paper I unclassified. CON vs non-detected Mann-Whitney p=0.23 (resolved), 0.09 (global): suggestive, not significant. N 0.99 (55 pc) / 0.69 (250 pc) matches the non-detected median; S 0.44 below all CONs.

## Literature
Falstad+2021 2021A&A...649A.105F (Table 3, criterion); Nishimura+2024 2024A&A...686A..48N (Tables 3, 4); Aalto+2015 background.

## Uncertain
- Our limit does not exclude moderate pumping (Tvib 150-300 K, alpha < 0.04, Q4); it excludes CON-like luminous HCN-vib.
- Falstad L_vib is whole-nucleus; Nishimura HCN is one beam or 10 beams: the literature vib/HCN is a 2-3x aperture-dependent range.
- Literature tables transcribed from garbled PDF text; checked by recomputing Falstad's Sigma for ESO 173 and NGC 3256 N. n=3 CONs in the Nishimura overlap; Falstad's Paper I lacks 6 Nishimura galaxies.
- Not done: 860 um (Band 7) continuum brightness separately (Xb47 346 GHz peaks N 36, S 22 mJy/beam at 4.2x2.8" are too coarse; T_b at 0.24" 1.1 mm already shows no compact core; ratio of Band 7/6 spectral indices in Q6).
