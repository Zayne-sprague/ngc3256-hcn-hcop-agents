# Q54: HNC/HCN, CN/HCN at N and S (chemistry)

## Did
- Literature: Ledger+2021 (CN(1-0)/HCN(1-0) ~1.7 at N, ~1.1 elsewhere), Harada+2018 Table 13 (HNC/HCN 0.56 N, 0.50 S at 1.7"), Meijerink+2007/Loenen+2008/Boger&Sternberg 2005 expectations (via lit.py snippets; model values not recomputed).
- Archive (requested and cached): CN 3-2 (2018.1.00493 X1df1 spw25), CN 1-0 (2018.1.00223 X4d spw27), HNC 3-2 (2015.1.00412 X169), HNC 1-0 (2015.1.00412 X149); already cached: HCN/HCO+/HNC 1-0 from 2015.1.00993 and 2015.1.00412 X145; HCN 3-2/HCO+ 3-2 from /data Xb4d; HCN 4-3 X1df1.
- Fixed-window Tb integrals at N/S in common circular beams (same recipe as Q53); CN integrated over the whole hyperfine complex (and main group alone); window/baseline variants for N/S double ratios.

## Scripts (/workspace/scripts) and outputs (/workspace/out/q54)
q54_chem.py (line table, HCN/HNC integrals -> chem.csv), q54b_cn_ratios.py (CN and pair ratios -> ratios.csv; importing runs it), q54c_double_ratios.py (variants -> double_ratios.csv), q54d_vs_hco.py (X/HCO+ -> vs_hco.csv), cn_spec.png (CN spectra at 1.9"/2.5").

## Results (full numbers in report.json Q54a-d)
- CN(1-0)/HCN(1-0), 2.5": N 1.76, S 1.39 (N/S 1.27+-0.03); reproduces Ledger (1.7 at N). HNC/HCN(1-0): N 0.52, S 0.42 (0.44/0.36 with 00412). 3-2: HNC/HCN ~0.6 at both (N/S ~1.0), CN/HCN(3-2) N 1.0, S 0.9 (1.1-1.4") falling to 0.6 at S at 2.5" (beam dependent).
- N/S of X/HCO+ (1-0): HCN 1.58, HNC 1.96, CN 2.00; (3-2): 2.0, 1.95, 2.5. CN/HCN and HNC/HCN do not rise at S.
- Conclusion: the S deficit is shared by HCN, HNC, CN (and CS, Q14e) relative to HCO+; no HCN-specific destruction signature (XDR/PDR would raise CN/HCN at S), HNC/HCN < 1 and N=S so no sign of different mechanical heating; excitation cannot explain the J=1-0 contrast. N is more CN-bright (PDR/UV-like) than S. Confidence low-medium.

## Caveats
- Cross-project flux scale (10%/project assumed) affects single ratios, not N/S double ratios. CN 3-2 at S and HCN 4-3 at S are window/beam sensitive; CN(3-2)/HCN(4-3) not used.
- 2016.1.00965 (0.55-0.65") cubes lose 25-45% flux and the S HNC disagrees with 2015.1.00412 (flagged, unused). No sub-arcsec chemistry ratios were made: CN 1-0 at 0.29" (2018.1.00223 X4b, 8.9 GB) plus HCN 1-0 at 0.65" would give the N core; CN 2-1 and HNC 3-2 beyond 266 GHz not in the products used.
- Opacity (HCN tau10 ~5-12; CN complex main:J=1/2 group ratio ~1.6 vs 2 thin, tau~0.8) and abundance are not separated; models not re-run; no H13CO+/H13CN-normalised columns.
