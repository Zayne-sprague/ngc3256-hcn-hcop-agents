# Q46 work summary: N and S HCN/HCO+ among published LIRG/AGN values

**Script:** /workspace/scripts/s300_q46_lit_context.py -> /workspace/out/q46/s300.txt (transcribed tables embedded; re-runnable). Lit text dumps: /workspace/lit46/.
**Report:** four entries Q46a-d plus a conclusion entry added to /output/report.json (backup report_backup_q46.json).

## Literature used
Nishimura+2024 (CON-quest II; Tables 1, 2, 4; 3-2, 23 U/LIRGs, 0.3" = 50-420 pc), Izumi+2016 (Table 1; 4-3 diagram), Harada+2018 (Table 13, 1.7"), Michiyama+2018 (Table 4, 2.1"), Ueda+2021 (quotes Privon+2015 GOALS 1-0 means), Imanishi+2016/2018/2023 and Saito+2018 (qualitative: AGN > SB at 3-2/4-3; VV 114 E AGN >1 flat in J), Ledger+2021 (NGC 3256 HCN 3-2/1-0 0.39 at N vs 0.27 elsewhere; no HCN/HCO+). Privon+2015 and Krips+2008 are not in the library; their values are second-hand. NGC 3256 is not in Nishimura's sample.

## Results (all T_b ratios)
- 3-2: Nishimura resolved median 1.13 (0.45-2.67, n=21), global median 1.05. N 0.99 (55 pc) = 38th percentile, mid-range; 0.69 at 250 pc = 19th. S 0.44 (55 pc), 0.37 (250 pc): at/below the minimum of all 21. For beams <110 pc (n=8) median 1.08 [0.56-2.34].
- No correlation with alpha_AGN in the sample (rho=-0.36, p=0.12; medians 1.32 for alpha<0.1, 0.79 for 0.1-0.25).
- 4-3 (Izumi): ours N 0.52, S 0.29 (800 pc, Xb47); N+S 10" 0.38 reproduces Izumi's APEX NGC 3256 value 0.39+-0.11. SB range 0.24-1.03 (median 0.48); AGN 1.1-3.2. N and S both SB-like.
- 1-0: Harada N 0.96/S 0.645; Michiyama 0.94/0.60; ours 0.84-0.94/0.57-0.58; GOALS mean 0.9-1.0, 0.77 without AGN.
- Conclusion: N ~1.0 is ordinary for SB-dominated LIRGs, not AGN/CON-like; S (0.4-0.6) is the unusual (low) nucleus, starburst-like despite hosting the AGN. HCN/HCO+ at one J and 50-800 pc fails as AGN diagnostic here; the literature placement does not by itself decide the cause (Nishimura: outflow ratio enhancement is abundance/shock-heating driven).

## Uncertain
Beams not matched (factor 2-4); N/S at 4.1x2.7" partly blended; alpha_AGN unreliable for obscured nuclei; Privon+2015 and Krips+2008 via secondary quotes; classification of N=SB, S=AGN taken from literature; Nishimura CON flags were garbled in the extracted text so a CON vs non-CON split was not done (answer to "CON-like?": ratio 1-2 occurs without luminous AGN; no CON-specific value could be tested).
