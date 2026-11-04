# r11 — N/S contrast of HNC/HCO+, CN/HCO+, HCN/HCO+ (T_b), one pipeline

Source: team-b Q52 pipeline (`scripts/s320_q52_chem_ratios.py` → `out/q52/s320.json`; ratios `scripts/s321_q52_ratios.py` → `out/q52/s321.{txt,json}`; table by `r11_table.py` → `r11_table.md/.csv`). Nothing re-run; this is the existing single-pipeline result, repackaged.

## Table (RJ T_b ratios; N/S = R_N/R_S)
| ratio (T_b) | 1-0 2.5" N | S | N/S | 3-2 1.4" N | S | N/S | 3-2 2.5" N | S | N/S |
|---|---|---|---|---|---|---|---|---|---|
| HNC/HCO+ | 0.474±0.013 | 0.249±0.010 | **1.90±0.11** | 0.477±0.011 | 0.259±0.014 | **1.84±0.14** | 0.395±0.009 | 0.175±0.013 | **2.26±0.21** |
| CN/HCO+ | 1.590±0.011 | 0.755±0.010 | **2.10±0.04** | 0.893±0.007 | 0.348±0.008 | **2.56±0.09** | 0.725±0.005 | 0.242±0.006 | **2.99±0.14** |
| HCN/HCO+ | 0.912±0.008 | 0.585±0.008 | **1.56±0.03** | 0.833±0.006 | 0.389±0.005 | **2.14±0.05** | 0.737±0.004 | 0.382±0.004 | **1.93±0.05** |
| HNC/HCN | 0.520±0.014 | 0.427±0.018 | **1.22±0.07** | 0.572±0.013 | 0.666±0.037 | **0.86±0.06** | 0.536±0.012 | 0.458±0.035 | **1.17±0.10** |
| CN/HCN | 1.743±0.014 | 1.292±0.020 | **1.35±0.03** | 1.072±0.008 | 0.894±0.022 | **1.20±0.04** | 0.983±0.007 | 0.633±0.018 | **1.55±0.05** |

Method (same code for every line/beam): pbcor cube → RJ T_b at native beam → unit-area kernel to a circular beam (2.5" for 1-0; 1.4" and 2.5" for 3-2) → beam-weighted value in the pixel at the fixed N/S position (not the per-line peak) → linear baseline from line-free channels → integral ±350 km/s about v_sys (N 2785, S 2750 km/s). Errors = line-free-channel noise only (no 10% flux-scale terms: they cancel in N/S and in ratios of lines within one band). ±250 km/s windows change the ratios by ≤0.03 and N/S by ≤0.06 (see `s321.txt`). 3.5" 1-0 values are in `s321.txt` (N/S: HNC/HCO+ 1.83, CN/HCO+ 1.98, HCN/HCO+ 1.48).

## Lines and CN components
- 1-0: HCN, HCO+ from 2015.1.00993.S X961 (spw0), HNC from X961 spw1 (HNC also from 2015.1.00412 X149: HNC/HCN N/S 1.23, HNC/HCO+ 1.91, consistent). CN 1-0 = **sum of both hyperfine groups**: G1 (N=1-0 J=1/2-1/2, 113.144 GHz) + G2 (J=3/2-1/2, 113.491 GHz), 2018.1.00223.S X4d spw27, each integrated over its own ±350 km/s window (G1 −650..+350 to include its hyperfine components), the other group masked in the baseline.
- 3-2: HCN, HCO+ from 2017.1.00379.S Xb4d (12m); HNC 3-2 from 2015.1.00412.S X169 (12m); CN 3-2 from 2018.1.00493.S X1df1 spw25 (whole N=3-2 complex: J=7/2-5/2 group at 340.248 GHz plus the J=5/2-3/2 group at 340.032 GHz, one window −350..+450 km/s about the 340.248 GHz reference, line at ≈+700 km/s masked). So CN is the full hyperfine complex in both transitions, not a single component.

## Caveats that matter for the paper
1. **The 3-2 entries are 12m-only, not flux-complete.** The 7m data were combined only for HCN/HCO+ (`s192`, 2.2" aperture, uv-matched with 1-0): there the HCN/HCO+ N/S falls from 1.65-1.68 (Xb4d image) to **1.44-1.50** (12m+7m). The table's 3-2 2.5" HCN/HCO+ contrast 1.93 (pixel estimator, 12m only) is therefore an upper estimate by ≈0.3-0.4 and not comparable with 1-0 1.56; the 1.4" value 2.14 is probably also inflated by a smaller amount (r07: −0.2 to −0.3 at 1.5"). No 7m data exist in the cache for HNC 3-2 (00412 X169) or CN 3-2 (00493 X1df1) and I did not restore them, so **HNC/HCO+ and CN/HCO+ 3-2 contrasts (1.84-2.26, 2.56-2.99) cannot be corrected here** and must be quoted as 12m-only upper values; S, the more extended source, loses more flux. Only the 1-0 column (2015.1.00993, X961 short-spacing-complete at 2.5") is a flux-complete one-pipeline set; for 3-2 the matched-uv statement can only be made for HCN/HCO+. A matched-uv (taper to common uv range) 3-2 set for HNC/CN would need restoring raw data.
2. The 1-0 and 3-2 columns use pixel values at the fixed nuclear positions, whereas r07/s192 use apertures centred on the HCO+ peak; the 1-0 HCN/HCO+ N/S is 1.56 (pixel) vs 1.45-1.47 (aperture, ±150/±350, 2.2" beam).
3. Window/centroid: 3-2 N/S ratios in the 1.4" row use native beams of 0.93" (HCN/HCO+), 0.75" (HNC), 0.92" (CN), all smaller than 1.4": OK. 2.5" (and 1-0 native 2.1") fine.

## Summary of what can be stated
- 1-0 at 2.5" (flux-complete): N/S = 1.56±0.03 (HCN/HCO+), **1.90±0.11 (HNC/HCO+), 2.10±0.04 (CN/HCO+)**; HNC/HCN 1.22±0.07, CN/HCN 1.35±0.03 → all three N/S contrasts are >1.5, and relative to HCO+ the S deficit is as large or larger in CN and HNC than in HCN (HCN is not special at 1-0 at this beam).
- 3-2 (12m-only, upper values): HCN/HCO+ 2.14±0.05 (1.4") and 1.93±0.05 (2.5"); HNC/HCO+ 1.84±0.14 / 2.26±0.21; CN/HCO+ 2.56±0.09 / 2.99±0.14. HNC/HCN N/S = 0.86±0.06 (1.4") and 1.17±0.10 (2.5") — consistent with 1 within the 12m-only uncertainty; CN/HCN 1.20±0.04 / 1.55±0.05.
- Relative to dust (1.1 mm continuum, `s322`): the S deficit applies to HCN, HNC and CN alike (p0058).
