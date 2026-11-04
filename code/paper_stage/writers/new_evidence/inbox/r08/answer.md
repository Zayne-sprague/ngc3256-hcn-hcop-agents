# r08 — S nucleus numbers reconciled (team-a; no new fits)

## (i) HCN/HCO+ 4-3 (T_B ratio), N and S
| value | N | S | beam / scale | array, window, baseline | source |
|---|---|---|---|---|---|
| 12-m cubes, fixed window ±350 km/s (N 2785, S 2750), linear baseline 450-950 | 0.81±0.02 (MC 0.80, 0.69-0.93) | **0.22±0.02 (MC 0.23, 0.18-0.29)** | 1.0" = 188 pc | 2018.1.00493 HCN x 2018.1.00223 HCO+ 12-m | team-a p0053; /snapshots/team-a-s30/scripts/q53g_fixedwin_all.py, q53h_mc.py |
| same, window/baseline range | | 0.14-0.24 (always <0.30) | 1.0" | S HCO+ has blue+red horns, HCN only blue (blue 0.27, red 0.13) | p0053 |
| same | 0.67 | **0.18** | 1.9" = 357 pc | 12-m | p0053 |
| ACA-visible, Xb47, peak pixel | 0.52 | **0.29** | 4.14x2.70" ≈ 800x520 pc | ACA only, ±350 km/s | team-b p0055 (/snapshots/team-b-s28/out/q46/s300.txt, line 31) |
| 4-3 extrapolated to the core | 0.95 (0.81-1.12) | 0.29 (0.25-0.35) | 56 pc | assumption-based | p0053 |
Note: the 0.29 at 56 pc and the 0.29 at 800 pc are different quantities that coincide by chance. The first is an extrapolation using an assumed trend, the second a measurement. Beams, arrays and baselines are not matched between rows. The 12-m 4-3 cubes resolve out part of the S HCN (S HCN recovery 61% vs HCO+ 87% at 4.5"; 12-m N/S 4-3 is 3.7 vs the ACA uv-plane 1.8 +0.3/-0.2; p0066, p0053), so the 12-m S values are lower limits on the S ratio of the full emission and the 0.18-0.22 values may be biased low.
**Quote:** at 190 pc, S = 0.22 (+0.07/-0.04, MC range 0.18-0.29; window/baseline spread 0.14-0.24), N = 0.81 (0.69-0.93), N/S ~3.7, as a 12-m value that excludes extended flux. At 800 pc, S = 0.29, N = 0.52, N/S = 1.8, consistent with the uv-plane ACA contrast 1.8 +0.3/-0.2 (p0066, out/q65/final.txt, tables.txt). Do not quote 0.29 at 56 pc: it is an extrapolation. Under both numbers S is below the Izumi+2016 SB median 0.48 and N is SB-like.

## (ii) S kinetic temperature
| statement | what it is | status |
|---|---|---|
| p0012: S 20-30 K, ≥80 K excluded | 1" one-zone RADEX joint fit | Retracted/revised: team-b p0012-r01 showed the exclusion rests on the deblended 6" S 4-3/3-2 values, which move by 0.2+ with template/window; swapping them allows 25-80 K. The "≥80 K excluded" claim must not be quoted. |
| p0051: Tc_S = 17±1 K | nested sampling, joint N+S fit with shared envelope, CMB background | Converged for S, but a model-conditional one-zone RADEX core temperature depending on band-scale prior and the envelope; replaces the SMC values (p0048) |
| p0037/p0040: Tk(S) ≥ 36-40 K (ff=1); 34-53 K smooth core; window-mean ≥30 K | CO(2-1) T_b at 0.18-0.26" (peak 31 K, velocity-matched 20 K; thick, co-spatial) | A lower limit, assuming optically thick and co-spatial CO with ff~0.6-0.9; clumping allows warmer (ff 0.5: ≥48-71 K). Confidence low. |
| team-b p0063: 80-130 K | the shared extended (10") slab; Tk_C (compact core) is unconstrained | Applies to extended gas, not the S core |

**Consistency.** CO limit ≥30-40 K, p0063 extended 80-130 K and the p0012 revision (25-80 K allowed) are mutually consistent (all warm or unconstrained from above). Tc_S = 17 K (p0051) and the retracted 20-30 K of p0012 are in tension with the CO limit, unless the dense HCN/HCO+ gas is a different, cooler or sub-thermal component than the CO emitting gas. They are one-zone, model-conditional values, and they should not be presented as a measurement of S.
**Defensible single statement:** the S kinetic temperature is not measured. CO(2-1) brightness gives Tk ≳ 30-40 K (a lower limit, thick and co-spatial CO, ff≥0.6); the HCN/HCO+ ladder fits give no secure value (17-30 K one-zone, 25-80 K allowed, extended gas ≥80 K), and the dense-gas Tk(S) is "unconstrained". Needed to decide: HCN 4-3/3-2 at ≤0.5" (team-b p0063: 0.34 cold vs 0.60 warm), or 13CO/C18O for the CO opacity.
