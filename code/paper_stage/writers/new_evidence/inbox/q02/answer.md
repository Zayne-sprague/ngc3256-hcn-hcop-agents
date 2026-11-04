# q02 — N and S HCN/HCO+ 3-2 separately, 12-m vs 12-m+ACA (team-a)

**Short answer.** The S ratio does **not** change by 9% between arrays in team-a's uv-matched images: 12-m → 12-m+ACA changes S by −2.0% at 1.0" and −0.8% at 1.5" (N by −0.3% / −1.5%). The N/S change is +1.7% (1.0") and −0.7% (1.5"). I find no team-a number that gives 9%. The only 9-11% shifts in the outputs come from the continuum-subtraction/baseline choice (N 0.867 → 0.966 at 1.0"), not from the array. I cannot see the paper's A.1, so I cannot tell which pair of numbers it compares.

## 1. N and S, 12-m vs 12-m+ACA (script `q02_table.py` → `q02_table.txt`; source `/workspace/out/q64/ratios.csv`)
Images: CASA tclean mfs of HCN and HCO+ 3-2, ±350 km/s about 2753 km/s, circular restoring beam, multiscale clean, pbcor (`scripts/q64c_image.py`, sets `main`, `var`; the restored measurement sets are Xb4d for 12-m and Xb4d+Xb4f for 12-m+ACA). Values are at the peak pixel of each nucleus. The error is the statistical term only (image rms / peak, in quadrature). Line values are mJy/beam.

| beam | data | N HCN | N HCO+ | **N ratio** | S HCN | S HCO+ | **S ratio** | N/S |
|---|---|---|---|---|---|---|---|---|
| 1.0" | 12-m (`d_u0_b1.0`) | 9.93 | 11.45 | **0.867 ±0.050** | 3.69 | 8.65 | **0.427 ±0.041** | 2.03 |
| 1.0" | 12-m+ACA (`dc_b1.0`) | 9.86 | 11.41 | **0.864 ±0.038** | 3.59 | 8.56 | **0.419 ±0.030** | 2.06 |
| 1.5" | 12-m (`d_u0_b1.5`) | 14.84 | 18.73 | **0.792 ±0.048** | 5.27 | 12.82 | **0.411 ±0.045** | 1.93 |
| 1.5" | 12-m+ACA (`dc_b1.5`) | 14.77 | 18.93 | **0.780 ±0.038** | 5.39 | 13.21 | **0.408 ±0.034** | 1.91 |

The mean ratio over r<0.25" gives the same picture (1.0": N 0.860/0.855, S 0.417/0.409; 1.5": N 0.791/0.779, S 0.410/0.408). The stat errors above are conservative, because rms/peak treats the noise as independent in each line. The out/q64r2 noise-based errors are about 0.01-0.02 on the ratio.

Systematics at 1.0" (same table): continuum-subtraction/baseline variant B gives N 0.966 and S 0.398 (12-m), N 0.960 and S 0.418 (12-m+ACA). Velocity window ±250/±500 gives N 0.839/0.879 and S 0.415/0.419. The two 12-m execution blocks differ by S 0.387 vs 0.428 and N/S 1.98 vs 2.24. So the baseline/window and execution-block scatter (S ±5%, N ±5-10%) is larger than the array difference.

## 2. Caveats on the comparison
- **Beam.** The matched images use 1.0" and 1.5" circular beams. There is no 1.05" product; 1.0" is the closest, and it is slightly below the native Xb4d beam (0.93×0.74", so it is a mild taper).
- **Weighting is not identical between the two sets.** The 12-m set is natural weighting with a 0.5" (1.0" beam) or 1.2" (1.5" beam) taper. The 12-m+ACA set is Briggs robust 0.5 with a 0.3" taper (1.0" beam) or the same 1.2" taper (1.5" beam) (`q64c_image.py` `jobs_for`). The uv range is unrestricted in both. The ACA adds short spacings, so this is the test of missing flux, but weighting differs. I did not re-run with identical weighting.
- **Independent check.** The pipeline images (Xb4b 12-m 0.25" convolved vs Xb4d 12-m, `out/q64r2/summary.txt`) give S 0.390 vs 0.393 at 1.0" and 0.373 vs 0.385 at 1.5" (+0.8%, +3.2%). N is 0.834 vs 0.850 and 0.759 vs 0.767. N/S is 2.14 vs 2.15 at 1.0" and 2.03 vs 1.99 at 1.5". Those are both 12-m only (Xb4b and Xb4d are different configurations, not an ACA test).
- Below 1" there is no 12-m+ACA test; Xb4b visibilities could not be restored.

## 3. 3-2 N/S at 1.4": HCN 3.11, HCO+ 1.45, ratio 2.14 vs "1.578"
- **3.11 / 1.45 / 2.14 are team-b's numbers** (p0056-r01; `/snapshots/team-b-s30/out/q52/s322.txt` for 3.11 and 1.45, `s321.txt` line "J=32 1.4" HCN/HCO+" for 2.14). They are N/S ratios of the integrated Tb line intensities in a 1.4" circular beam at the fixed N and S positions, ±350 km/s, from the Xb4d 12-m cubes (script `s320_q52_chem_ratios.py`/`s321_q52_ratios.py`). I read this from the snapshot output; I did not re-run it. Team-b's N and S ratios are N 0.833±0.006 and S 0.389±0.005, N/S 2.14±0.05.
- **Team-a's own 1.4" estimate is 2.006 ±0.044** (N 0.791±0.007, S 0.394±0.008; Xb4d 12-m, p0056, `/snapshots/team-a-s31/out/q54/vs_hco.csv` row `3-2,1.4`). The 2.14 vs 2.01 difference (6%) is mostly N (0.833 vs 0.791), from the baseline/window convention. Our own 1.4" matched-beam values (`out/q53/matched32.csv`, Gaussian-fit line integrals, ±1.3×fit error) give 2.02 (Xb4b: 0.770/0.381) and 2.03 (Xb4d: 0.776/0.382).
- **1.578 is not a 3-2 12-m 1.5" estimator.** It is the HCN/HCO+ **1-0** N/S at **2.5"**, from the 2015.1.00993 cubes: N 0.906±0.010, S 0.575±0.013, N/S 1.577±0.04 (p0056, `out/q54/vs_hco.csv` row `1-0,2.5`; team-b got 1.56 for the same). The 3-2 value at 1.5" from the 12-m only is 1.93 (peak) / 1.99 (q64r2 Xb4d, 1.5"), and 1.91 with ACA. If the paper compares 2.14 with 1.578, it compares different lines (3-2 vs 1-0), beams (1.4" vs 2.5") and estimators, so the difference is not an A.6 inconsistency within the 3-2 data. Please fix the label.
- Also: 3.11 and 1.45 are N/S ratios of single-line intensities. They are not HCN/HCO+ ratios and are not comparable to the 1.578 number.
