# Q30 work summary
Answer: yes, tau(12CO 2-1) is constrained, given an assumed 12C/13C. 13CO(2-1) (2015.1.00412.S spw0) is detected at both cores, S/N 18-70.
Matched 1.0" 12CO/13CO T_b ratio: S 10.8 (peak ch) / 13.3 (HCN-FWHM mean); N 13.6 / 15.2 (+-20% systematic). At 2.5": S 15.9, N 14.6.
tau(12CO): S 3.8/5.8/8.8 (peak) for R=40/60/90; N 2.9/4.5/6.9; R=25 gives 2.1 / 1.5. So tau>=3 holds unless R<~35.
Tk_min (ff=1, Q23 0.26" T_b, R 40-90): S 37 K (peak) / 26 K (HCN-mean); N 19 / 17 K. ff fit 0.88-0.96: S 41, N 20. ff=0.5: S 68, N 33; ff=0.3: S 110, N 51. The Q23 bound is confirmed with tau measured. ff is the open parameter.
CO(1-0), CO(3-2) at 2.5" (1-0 native beam 1.92x1.68" so 1" matching is impossible): r32/21 ~0.78-0.81, r10/21 ~1.16 -> thermalised thick CO, Tex ~12-16 K for the 2.5" gas. Consistent, but it does not probe the compact hot component.
Scripts: scripts/q30_13co.py (convolve+sample), q30_ratio.py (1" ratios), q30_other_co.py, q30_final.py (tau, Tk, CO(1-0)/(3-2)), build_report_q30.py. Outputs: out/q30/ (ratios.json, tau_tk_table.csv, other_co_2p5.json, spectra).
Literature: lit.py search found no 12C/13C or 13CO for NGC 3256 (Harada+2018 has none; Sliwa+2017 only other LIRGs). R=25-90 is an assumption (NGC 253 ~25-40 from Martin+2019, Arp 220-like >~100 from memory, unverified).
Uncertain: R; 13CO/12CO same Tex and ff (13CO more extended -> tau lower limit logic); missing flux differs between projects; CH3OH 220.078 GHz blend in the 13CO red wing not excluded; my 1" S peak T_b (12.6 K) is 6% above the Q23 1.7"-fit-table value (11.8 K) from sampling position/pixel choice; 2.5" S/N values are formal (correlated noise). Forum replies to p0026/p0028/p0011 not posted.
