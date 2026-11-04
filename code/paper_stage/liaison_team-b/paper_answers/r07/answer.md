# r07 — Adopted N/S HCN/HCO+ 3-2 contrast, error budget, and beam dependence

All numbers are **T_b (RJ) ratios, HCN 3-2 / HCO+ 3-2 at equal circular beams** (the nu^2 factor is applied), aperture r = beam FWHM centred on the HCO+ 3-2 moment peak (within 0.4" of nominal), noise from line-free channels. Contrast = R_N / R_S.

## One-line statement
**At 50-120 pc (0.26"-0.62") the N/S contrast is 2.1 +0.1/-0.2 (adopted 2.1 ± 0.2; aperture estimator 2.12-2.32, other estimators 2.04-2.24). It falls smoothly with beam size to ≈2.0 at 1.0", 1.9 at 1.29" (250 pc), ≈1.6-1.8 at 2.2" with the compact-array data, and 1.44-1.50 at 2.2" with uv-matched 12m+7m data.** The 1.45 (1-0) and 1.50 (3-2) at 2.2" and the 2.1 at 0.26" are therefore different measurements of a beam-dependent contrast, not a conflict (details below).

## Table (script: `r07_budget.py`; inputs `out/q26/s170_*.json`, `s400_*.json`, `out/q34/s192.json`; `budget_table.md`)
| beam | aperture (Xb4b) | stat | pixel | peak-chan | Xb4d aperture | 12m+7m aperture (±350 / ±150) | uv term |
|---|---|---|---|---|---|---|---|
| 0.26" (50 pc) | 2.26 | 0.06 | 2.09 | 2.13 | – | – | ≤0.06 (assumed, see below) |
| 0.28" (55 pc) | 2.32 | 0.06 | 2.14 | 2.24 | – | – | ≤0.06 |
| 0.62" (120 pc) | 2.12 | 0.04 | 2.24 | 2.04 | – | – | ≤0.06 |
| **1.0" (194 pc)** | **1.99** | 0.04 | 2.16 | 1.84 | 1.94 | not made | −0.06 (Xb4d−Xb4b) |
| 1.29" (250 pc) | 1.90 | 0.04 | 2.10 | 1.83 | 1.84 | not made | −0.06 |
| 1.5" (291 pc) | 1.83 | 0.04 | 2.04 | 1.83 | 1.81 (s189) | 1.53 / 1.63 | −0.3 / −0.2 |
| **2.2" (427 pc)** | **1.64** | 0.03 | 1.88 | 2.00 | 1.65 | **1.44 / 1.50** | −0.20 / −0.14 |
| 6.5" blends | – | – | – | – | 1.09 (pixel 1.33) | – | Xb4f 1.06 |

The 1.0" and 2.2" rows are new, run with the same code as s171 (`scripts/s400_r07_extra_beams.py` = copy of s170 with extra beams; logs `s400_Xb4b.log`, `s400_Xb4d.log`; the Xb4b/Xb4d 2.2" values agree, 1.64 vs 1.65).

## Error budget
1. **Statistical**: 0.03-0.06 (MC of the line-free-channel errors, `s171_summary.txt` Table 2: 2.25 [2.20, 2.32] at 0.26").
2. **Estimator spread** (aperture vs same pixel vs peak channel): ±0.1 (half range) at 0.26-0.62"; ±0.15 at 1.0-2.2" (range 1.84-2.16 at 1.0", 1.64-2.00 at 2.2"). The aperture value sits at the top of the range at 0.26-0.28", the pixel value at 0.62". This is the "1.9-2.3" quoted in Fig. 2.
3. **Flux scale (14% = 2×10%)**: **cancels** in the N/S contrast. N and S are in the same two spectral windows, observations, and calibrations, so the scale multiplies both. It applies only to individual ratios (R_N 0.99 ± 0.14, R_S 0.44 ± 0.06 at 0.26"), not to the contrast (`s171` docstring, `s189` docstring).
4. **uv recovery** (the only systematic that does not cancel, because the missing flux differs between N and S; S is more extended/faint and loses more): 
   - 1.0" and 1.29": Xb4d (more compact array) vs Xb4b gives −0.06 (1.99→1.94; 1.90→1.84): small.
   - 1.5" and 2.2": 12m+7m (uv > 7.7 kλ, `s192`) gives a contrast lower than the 12m-only images by 0.14-0.30 (1.83→1.53/1.63 at 1.5"; 1.64→1.44/1.50 at 2.2"). The reason is that Xb4d recovers S 1.36× less than the compact array at 6" (R3 = Xb4f/Xb4d: N 1.16, S 1.36, `s189`), which inflates the Xb4b/Xb4d contrast at ≥1.5". The Xb4d recovery-bound-corrected 3-2 contrast at 2.2" is 1.56 (±150: 1.59).
   - **At ≤0.62" no matched-uv test exists** (the 12m+7m images were made only at 1.5" and 2.2"). I assumed ≤0.06 downward from the 1.0-1.29" Xb4d-Xb4b comparison, but a bias as large as that seen at 1.5" (−0.2 to −0.3) cannot be excluded by the data. Flag this as the main untested systematic of the 2.1.
5. **Window**: ±350 vs ±150 km/s about 2780 changes the 12m+7m 2.2" contrast by 0.06 (1.44 vs 1.50); at 1.5" by 0.10.
6. **Baseline**: constant baseline (per-pixel median of line-free channels), not linear (see correction below). Blue-only vs both-sided line-free channels changes R by ≤0.01 (`s171` Table 4).

**Total (adopted)**: 0.26-0.62": 2.1 ± 0.2 (stat ≈0.05, estimator ≈0.1, uv ≤0.06 downward, quadrature ≈0.12; widened to ±0.2 to include the untested uv term). 1.0": 2.0 ± 0.2 (+0, −0.06 uv). 1.29": 1.9 ± 0.15. 2.2": 1.5 ± 0.1 (12m+7m, ±150 km/s; 1.44 at ±350).

## Which numbers are T_b ratios in the same velocity window?
- Table above and Fig. 2: ±350 km/s about 2780 km/s for N, S, Bridge, same window in HCN and HCO+ (HCN regridded to the HCO+ velocity grid), T_b ratio.
- The p0062 values "1.45 (1-0) and 1.50 (3-2)" at 2.2" use **±150 km/s** about 2780, both lines, circular 2.2" beams: 1-0 1.448 ± 0.015 (2015.1.00993 X961); 3-2 1.503 ± 0.017 (12m+7m `s192`). At ±350 the same data give 1.467 ± 0.022 and 1.441 ± 0.023 (`out/q34/s192.txt`, `out/q29/s189.txt`). So the quoted 1.45 / 1.50 are window-dependent at the ±0.03-0.06 level.
- The 2.2" Xb4b/Xb4d 3-2 values (1.64-1.65 at ±350; 1.65-1.68 in `s189`) are *not* uv-matched to the 1-0 data and should not be compared with 1-0 1.45.

## Consistency: 2.1 at 0.26" vs 1.45-1.50 at 2.2" vs 2.08 at 1.1" (team-a)
- The beam dependence is monotonic in one pipeline: 2.26 (0.26") → 2.12 (0.62") → 1.99 (1.0") → 1.90 (1.29") → 1.83 (1.5") → 1.64 (2.2"), plus a uv-matched drop to 1.44-1.50 at 2.2".
- **The cause is the N core**: R_N falls 0.99 → 0.84 → 0.69 from 55 to 250 pc (power-law index −0.23 ± 0.03), whereas R_S stays flat (0.43 → 0.37, index −0.10 ± 0.03) (`s171` Table 5). This is the dilution of the HCN-rich ~100 pc core (ratio ≈1 within 30 pc, `s171` Table 9; profile 0.97 inside 0.25" to 0.39 at 2-3" from the radial analysis) into an envelope with ratio ≈0.4-0.6. I have not re-derived that radial profile here.
- **Not all of the 2.26 → 1.50 change is beam dilution**: with the same Xb4b pipeline the 2.2" contrast is 1.64, and the further −0.14 to −0.20 to 1.44-1.50 comes from uv matching (item 4). Pure beam dilution accounts for ≈ 2.26 → 1.64; the remainder is flux recovery of the compact-configuration data.
- Team-a's 2.08 at 1.1" lies between my 1.0" values (aperture 1.99 ± 0.04, pixel 2.16 ± 0.04, Xb4d 1.94 / 2.06) and the 1.29" aperture value (1.90) and within my estimator range (1.84-2.16). I have not seen their pipeline, so I only claim consistency with this range, not identical methods.
- 6.5" blends (Xb4d 1.09-1.10, Xb4f 1.06-1.07) contain N+Bridge+S in one beam with differential uv loss; they are consistent with the trend but are not a clean continuation of it.
- The 12m+7m contrasts at 1.5" (1.53-1.63) and 2.2" (1.44-1.50) are the ones to compare with the 1-0 contrasts (1.62 ± 0.05 at 1.5"; 1.45-1.47 at 2.2"): the 3-2/1-0 ratio of the contrasts D is 0.95-1.04 (`s192`), i.e. no J-dependence of the contrast at matched resolution.

## Correction to the request's caption facts (also relevant to r01/r02)
The request says "linear baseline". The pipeline behind Fig. 1/2 (s170/s171) uses a **constant baseline = per-pixel median of line-free channels** (v_rel < −450 or > 600 km/s), not a linear fit. Baseline choice (both sides vs blue side only) changes R by ≤0.01.
