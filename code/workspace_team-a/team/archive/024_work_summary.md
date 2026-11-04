# Q54 work summary (round 3)
**Done**: scripts/q54_chem.py, q54b..q54i (outputs in out/q54/), new q54j_cn_subarcsec.py -> out/q54/cn_subarcsec.csv.
- Adopted HNC/HCN(1-0) with systematic spread (N 0.48+0.05-0.08, S 0.39+0.04-0.05); 00993 vs 00412 offset ~15% parallel at all beams (flux scale/bandpass + HNC extended-flux loss); N/S 1.24 robust.
- Missing-flux tests: CN X4f(compact)/X4d = 1.01 (both nuclei) -> X4d complete. HCN 00965 (0.7") recovers only 72% N / 56% S vs 00993. CN X4b is 9% (N) / 23% (S) above X4d at 2.5" (calibration offset).
- Sub-arcsec: CN/HCN(1-0) at ~1" after correction ~1.4-1.6 both nuclei (N/S ~0.9-1.1) vs 1.27 at 2.5"; low confidence, set by missing-flux correction.
- RADEX joint posterior: N(HNC)/N(HCN)~1 (+-0.25 dex), N(CN)/N(HCN) 0.6-1, tau10(HCN) 4-11, T_N~85 K, T_S~30 K; N/S column-ratio differences <0.3 dex.
- Conclusion softened: composition (HCO+ excess at S) preferred, not demonstrated; excitation alternative kept; tests listed in Q54d (H13CO+/C18O/C34S).
**Literature**: Ledger+2021, Harada+2018, Loenen+2008 (via Bemis+2019), Meijerink+2007, Oteo+2017, Wilson+2018.
**Not done / uncertain**: no H13CN/H13CO+ cubes in cache; quantitative PDR/XDR grid tabulation remains from literature text (out/q54/grid_*.npz are RADEX grids, not chemical models); recovery factor beam-dependence at 1" unknown; no uv-matched re-imaging (only image-plane cube comparisons).
