# Work summary: Q1 (HCN/HCO+ at matched J), worker turn 17

## What I did
All scripts are in /workspace/scripts. Outputs are in /workspace/out. Re-run with `python3 -I scripts/<name>`.

| Script | Purpose | Output |
|---|---|---|
| s30_matchedJ_ratios.py | Crop each pbcor cube, subtract a per-pixel line-free baseline, convolve HCN and HCO+ to a common circular beam (10", 6.5") or keep native, integrate over \|v-VSYS\|<350 km/s (windows 250/500 as tests), read at region positions and in apertures | s30_ratios.txt/json, s30_maps_*.npz |
| s31_zone_decomp.py | Deblend N and S at J=2-1 and 4-3 by fitting zone templates (N, S, Ext) built from the Xb4d 3-2 map | s31_zones.txt/json |
| s32_sblue_wing.py | HCN/HCO+ of the Sblue blue-wing blob in fixed apertures, with noise from 400 random line-free positions | s32_sblue_wing.txt/json |
| s33_ratio_ladder_radex.py | RADEX (LVG, LAMDA) screen of the three-J ladder for density, temperature and abundance ratio | s33_radex_q1.txt/json |
| s34_xb4b_native_ratio.py | Native 0.25" and 0.94" 3-2 ratios | s34_native_ratios.txt |
| s35_update_report_q1.py | Writes the Q1 entries into /output/report.json (idempotent; report validates) | |

Conventions:
- Ratios are Jy km/s over Jy km/s. Multiply by 1.0126 for K km/s.
- Error is statistical, plus 5% relative HCN-vs-HCO+ calibration.
- p0003 showed that filtering cancels in HCN/HCO+ at N, S and Bridge, so no 10-35% recovery term is applied to these ratios. It does apply to the faint clumps.

## Results
1. **3-2 at native resolution, per region.**

   | Region | 0.25" (Xb4b) | 0.94" (Xb4d) |
   |---|---|---|
   | N | 1.03±0.03 | 0.88–0.90 |
   | S | 0.48±0.02 | 0.39–0.40 |
   | Bridge | 0.56±0.04 | 0.48–0.50 |
   | E1 | 0.25±0.06 | 0.30–0.37 |
   | E2 | 0.48±0.08 | 0.40–0.51 |
   | Sblue | 0.77±0.14 | 0.71–0.75 |

   - The Xb4d range is the region position vs the HCO+ peak pixel within 1".
   - The Xb4b W clump is not significant (HCN S/N 2).
   - N rises with resolution: 0.60 at 5.7×3.5", 0.90 at 0.94", 1.03 at 0.25". A compact HCN-bright core is diluted in larger beams.
   - S barely changes: 0.38, 0.39, 0.48.

2. **Matched beam, all J, 6.5" circular** (N and S are 5.2" apart, so each reading carries about 15–20% of the other nucleus).

   | Position | R(3-2) Xb4d / Xb4f | R(4-3) Xb47 |
   |---|---|---|
   | N | 0.544 / 0.540 | 0.425 |
   | S | 0.393 / 0.420 | 0.315 |
   | Bridge | 0.480 / 0.484 | 0.358 |

   R(4-3)/R(3-2) is about 0.75–0.80 everywhere.

3. **Matched beam, all J, 10" circular, N-S midpoint.**
   - R(2-1) = 0.581, R(3-2) = 0.451 (Xb4f) / 0.428 (Xb4d), R(4-3) = 0.331.
   - The ratio falls about 25% per J step.
   - The r<8" aperture gives 0.562, 0.415 / 0.377, 0.290.
   - At 10" the N and S readings are blends. Do not read them as nuclear values.

4. **Deblended N vs S** (zone templates from Xb4d 3-2).
   - Check at 3-2 (Xb4f, 6.5"): N 0.51±0.05, S 0.36±0.05. The template values are 0.558 and 0.375, so the method recovers them within 10%.
   - R(4-3): N 0.39±0.06, S 0.24±0.07 (6.5").
   - R(2-1): N 0.63±0.08, S 0.48±0.06 (10").
   - S/N ratio of HCN/HCO+ is about 0.6–0.75 at every J. **S is never enhanced relative to N.** This agrees with Harada+2018 (no enhancement at the AGN).
   - Confidence is medium. chi²_red of the fits is 35–240, because the zone templates are rough morphology. The quoted errors are inflated by that.

5. **Outflow blob (Sblue, wing -400 to -150 km/s), Xb4d native, fixed aperture.**
   - Wing R = 1.43±0.33 at r=0.5", and 1.38–1.49 at larger radii. HCN is detected at 6.9σ and HCO+ at 5.4σ.
   - Core velocities at the same position: R = 0.5–0.6. S-nucleus core: 0.38±0.02.
   - The wing is about 3σ above the S core. Confidence is tentative.
   - Team-a's 1.8±0.3 (p0001) used peak pixels picked independently per line. That is biased high. Our fixed-aperture value is 1.4±0.3, so we qualitatively agree on an enhanced wing but not on the size.
   - There is no red-wing flux at Sblue (<0.07 Jy/b km/s).
   - Wing flux is uncleaned and filtered, so it is a lower limit. The HCN/HCO+ ratio is probably less affected, but this is untested for the blob.

6. **RADEX on the 10" midpoint ladder** (R_K = 0.59, 0.45, 0.34).
   - 31 of 56 (T, n) cells fit (Tkin 15–200 K, n 1e4–1e8).
   - Abundance ratio A = N(HCN)/N(HCO+) is anticorrelated with n: A ≈ 3–4.5 at n ~1e4–1e5, A ≈ 1 at n ~1e6, A ≈ 0.6–0.9 at n ≥ 3e6.
   - The ladder cannot separate abundance from excitation. H13CN opacity (Q4) and the full Q3 analysis are needed.
   - Fits to the 6.5" N and S pairs have only two ratios and are not constraining.

7. **Clumps** (E1, E2, W).
   - At 10" their ratios are 0.5 (2-1), 0.3 (3-2) and 0.05–0.19 (4-3), i.e. a steeper decline with J than at the nuclei.
   - They are blended with the nuclear wings, and the Band 7 primary-beam correction at r~8–10" (PB ~0.5) is not in the error budget.
   - Treat the clump values as low confidence.

## Still uncertain
- **No clean "disk" ratio.** Annulus fluxes at 8–20" are negative in Xb4f and Xb47 (resolved-out bowl), so the extended-disk ratio is not measurable. Only the discrete clumps (E1, E2, W) are available, and they are noise-limited at 4-3.
- **Unresolved nuclei at J=2 and 4.** The deblending assumes the Xb4d 3-2 zone morphology holds at J=2 and 4. The N and S coefficients are anticorrelated (-0.5 at 6.5", -0.7 at 10"). The 10" R(4-3) for S is undetermined.
- **Sblue wing.** About 3σ and tentative. The pointing of the aperture centre is ±0.2". Whether the wing is filtered like the core is untested.
- **Abundance vs excitation.** Needs H13CN (Q4). The RADEX run is a screen, not the Q3 result.
- **Xb4f may itself miss flux.** Not yet compared with Xb47/Xb55 at a shared beam (Q3).

## Suggested next
Q4 (H13CN limits, opacity), then Q3 (full RADEX with the 10" ladder and a 10–35% cross-band systematic).
