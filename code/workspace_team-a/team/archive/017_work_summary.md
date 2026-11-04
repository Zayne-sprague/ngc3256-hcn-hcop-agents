# Q9 work summary (round 2, after reviewer "fix")

## What was done (all scripts re-runnable: `cd /workspace/scripts; python3 -I <script>`; caches in out/q9b/cache)
| Fix | Script | Outputs |
|---|---|---|
| 1 unmasked window-integrated ratios, null + baseline-fit errors, limits | `q9b_lib.py`, `q9b_windows.py` | `out/q9b/windows.csv`, `windows.log` |
| 2 beam/spatial-mixing test (matched beams, slopes, centroids, injection) | `q9b_beam.py` | `beam.log`, `beam_slopes.csv`, `centroids.csv`, `injection.csv` |
| 3 N/S separation of the 6.5" 4-3/3-2 ladder (Xb4d weights, forward model) | `q9b_ladder.py` | `ladder.log`, `ladder_bins.csv`, `ladder_fit.csv` |
| 4 RADEX + Tb per channel | `q9b_radex.py`, `q9b_tb.py` | `radex.log`, `radex_allowed.csv`, `radex_grid_ratios.csv`, `radex_{hcn,hco+}.npz`, `tb.log`, `tb_scatter.png` |
| 5 contamination (catalogue), aperture-centre shift, array recovery | `q9b_windows.py` (centre-shift block) | `centre_shift.csv`; catalogue from `out/inband_lines_curated.csv` |
| 6 wording | `build_report_q9b.py` | `/output/report.json` (Q9 entries replaced: 6 entries; backup `out/report.backup.q9b_pre.json`) |

## Results (HCN/HCO+ 3-2, Jy ratio; Tb = x1.0126; windows relative to VSYS 2780 km/s)
* **N, 0.30" r<0.25"**: core |v|<60 0.982+-0.013; -60..0 1.04+-0.02, 0..60 0.95+-0.01, 60..100 0.92+-0.02 (~12% blue->red decline); wings 100-150 1.05+-0.09 (red) / 1.04+-0.14 (blue); 100-260 stacked 1.05+-0.13. No significant wing rise at 0.30".
* **N, Xb4d 1.05" r<0.5"**: core 0.769+-0.004; 100-150 0.86+-0.05 (blue)/1.13+-0.08 (red); 100-260 stacked **1.19+-0.10** (4.2 sigma above core); 150-250 stacked HCN 7 sigma, HCO+ 1 sigma -> R>1.2 (2 sigma lower limit). Same pattern r 0.5-1" (1.15+-0.11) and 1-2" (1.0+-0.26).
* **Beam test**: slope per 100 km/s (40-km/s bins, +-140, unmasked): Xb4d@1.05 +0.082+-0.006; Xb4b->1.05" +0.072+-0.008; Xb4b@0.30" -0.003+-0.013. Same data, only the beam changes -> beam effect, not an array-recovery effect (array fluxes agree 94-105%, ratios within 0.5-2% in cores). Injection of the measured radial ratio profile with identical kinematics gives +0.049 (60-70% of the slope); residual 0.02-0.03 unexplained. HCN centroid is 0.04-0.13" nearer to N than HCO+ (dx +0.05 blue -> -0.125+-0.034" at +120 km/s, Xb4d).
* **S**: 0.375 (core) -> 0.44+-0.03 stacked 100-260 (1.05"); 0.414 -> 0.44+-0.05 at 0.30"; HCN detected out to ~260 km/s: deficit (2-3x below N) present at all velocities with HCN detected; beyond that only HCO+-limited. **Bridge**: core 0.439, wings 100-260 stacked 1.35+-0.32 (2.9 sigma rise).
* **Ladder at 6.5"** (unmasked): HCN 0.409+-0.003 core, 0.386+-0.036 red wing; HCO+ 0.500+-0.002 core, 0.396+-0.028 red wing. Xb4d weights: HCO+ red wing is 81% S flux. Forward model R_N/R_S: HCN 0.459/0.303, HCO+ 0.575/0.402 (chi2/dof 6-16: poor); explains ~2/3 of the HCO+ red-wing drop (model 0.434 vs 0.396+-0.028).
* **Physics (RADEX, equal filling)**: N core ratio 0.99 with per-channel Tb(HCN)/Tb(HCO+)=0.96-0.97 = two optically thick (tau>3), thermalised, equal-Tex lines (ratio->1 for X=N(HCN)/N(HCO+)<=3); peak Tb 1.5-2.8 K => filling 0.02-0.35. S (0.42): impossible for both thick with equal Tex (0/541 models); needs X<=1 with moderate HCN tau or Tex(HCN)<Tex(HCO+) (low n branch whose ladder 0.28/0.43 matches the fitted S R_S 0.30/0.40). 1.05" wing value 1.2 needs X>=3 or lower HCO+ opacity.
* **Contamination**: catalogue blends listed in report (CH3NH2 +144..160 on HCN red; CH3OH 9(0)-8(1) +171, CH3CHO -186 on HCO+). HCN excess is on both blue and red sides at N, which a red-only blend cannot explain. Aperture shifts +-0.1" change ratios <=0.02 (core) and <=0.07 (wings). HCN-HCO+ peak offset N 0.00", S 0.047".

## Remaining uncertain
* Whether the residual +0.02-0.03/100 km/s (not reproduced by injection) is a real HCN wing excess; HCN-wing blends not quantified (no individually detected CH3NH2/CH3CHO).
* Ladder N/S decomposition: poor chi2, depends on Xb4d recovery (78-80%) and an unmodelled extended component; cannot separate N/S directly at 6.5".
* RADEX is single-zone equal filling; X, n, T, filling degenerate; RADEX unreliable above tau~1e3 (cut at logN<=17.5/150 km/s).
* Wings beyond +-150 km/s are noise-limited in at least one line (limits only).
* Earlier q9_vel.py (masked) outputs are retained in out/q9 but their 'flat' wording is superseded.
