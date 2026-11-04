# Q53 round 3b — work summary (worker)
Reviewer fixes done; scripts re-runnable, outputs in out/q53/.

## What I did
- **Fixed common velocity window** instead of Gaussian fits (`scripts/q53g_fixedwin_all.py` -> `out/q53/fixedwin_all.csv`; 3-2 and 4-3, beams 0.3-2.5", all/blue/red halves). Window N 2785, S 2750 km/s ±350, linear baseline 450-950 km/s, same for both lines/J. The spectra are double-horned (`out/q53/spec43_overlay.png`, `q53f_plotspec.py`), so the Gaussian fits were wrong: HCO+ N was overestimated by 13%.
- **uv/missing-flux test** (`q53e_fixedwin43.py` -> `fixedwin43.csv`): new X1df1 (HCN) and X858 (HCO+) cubes vs same-project Xb47 4.1x2.7" cubes at 4.5". Recovery N HCN 88%, HCO+ 77%; S HCN 61%, HCO+ 87%. Ratio at 4.5": N 0.59 vs 0.52, S 0.22 vs 0.32. Adding the missing flux as a smooth additive T_B shifts the 1.0" ratios N 0.81->0.79, S 0.22->0.24. No visibilities re-imaged (not restored); assumption of a smooth missing component is unverified.
- **Pointing/peak** (`posoff43.csv`): 5x5 offsets ±0.4" at 1.0": ratio at line peaks equals the fixed position (N 0.81/0.79, S 0.22); rms over ±0.2" 0.02/0.016.
- **Error budget + MC** (`q53h_mc.py` -> `mc.json`): stat, 10% flux scale per project, missing flux, pointing, window, and extrapolation (flat..linear R43/R32).
- `report.json` updated (`q53i_report.py`; backup `report.pre_q53r3b.json`): old matched/core/classification/verdict entries replaced, 3 new entries added (uv test, S line shape, R43/R32 vs beam).

## Results (T_B ratio HCN/HCO+)
| scale | J | N | S |
|---|---|---|---|
| 56 pc (0.3") | 3-2 | 0.98±0.03 | 0.49±0.02 |
| 188 pc (1.0") | 3-2 | 0.84±0.01 | 0.40±0.01 |
| 188 pc | 4-3 | 0.81±0.02 stat; MC 0.80 (0.69-0.93) | 0.22±0.02; MC 0.23 (0.18-0.29) |
| 357 pc (1.9") | 3-2 / 4-3 | 0.72 / 0.67 | 0.37 / 0.18 |
| 56 pc core | 4-3 (extrapolated) | 0.95 (0.81-1.12) | 0.29 (0.25-0.35) |
- Izumi SB ≤1.03 / AGN ≥1.11. N 188 pc: P(SB)=0.95, P(AGN-like)=0.014. **N core: P(≤1.03)=0.68, P(≥1.11)=0.18 -> 4-3 diagnostic not discriminating at the N core**; clearly SB-like at ≥188 pc. S: P(below SB median 0.48)=1.0 at all beams: false negative. 3-2 N core 0.98 is at the low edge of the ULIRG range (Imanishi+2019 median 1.3).
- R43/R32 is not beam independent: N 0.97->0.89, S 0.56->0.45 from 1.0 to 2.5" (partly 4-3 missing flux). Earlier "N 0.88 core" superseded.
- S line shape: HCO+ blue+red horns, HCN only the blue (4-3 ratio S blue 0.27, red 0.13; N 0.82/0.81). Explains the 25 km/s centroid and 120 vs 192 km/s width mismatch; not vib-HCN (1480 km/s away), absorption or an H13CN blend. S ratio depends on window/baseline (0.14-0.24, always <0.30).
- N/S contrast persists at 357 pc (3-2 2.0, 4-3 3.7), not core only. Q47: tau10(HCN)~5-12, X_S/X_N unconstrained (median ~1, 95% 0.45-2.6) -> ratio not an abundance measure; XDR enhancement at S untested.
- 1-0 values (N 0.96/0.94, S 0.645/0.60) are Harada+2018 Table 13 and Michiyama+2018 Table 4 (`out/q53/lit_table.csv`), not measured here.

## Literature
Izumi+2016 Table 1 (4-3 SB 0.24-1.03, AGN 1.11-3.19); Imanishi+2019 (3-2 ULIRGs); Harada+2018; Michiyama+2018; Privon+2015/2017 (1.84 is a Privon 1-0 AGN mean, not an Izumi threshold).

## Uncertain
- Not re-imaged in uv; X1df1/X858 missing flux handled only as a smooth component. HCN 4-3 S flux is sensitive to baseline (spectral coverage ends +540 km/s from HCN).
- Core 4-3 is extrapolated; the beam trend is partly instrumental. Flux-scale 10% per project is assumed.
- The two projects' pointing centres/PB corrections not checked separately (inside the pointing term only via source-position offsets).
