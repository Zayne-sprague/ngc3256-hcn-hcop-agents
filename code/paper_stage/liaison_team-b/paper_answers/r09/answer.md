# r09 — HCN/HCO+ 3-2 T_b ratio vs radius, N and S (Xb4b vs Xb4d)

Script `r09_profile.py` (re-runs the s170/s171 convolution/baseline/moment code, exec'd from `scripts/s170_q26_ratios.py`) → `r09_profile.json`; figure `r09_fig.py` → `fig_ratio_profile.png`. This re-derives the profile (not previously re-derived since r07).

## Definitions
- T_b (RJ) ratio = HCN/HCO+ with ν² applied; moment-0 over ±350 km/s about 2780 km/s, constant per-pixel baseline (line-free v<-450, >600), HCN regridded to HCO+ velocities; flux-conserving convolution to a circular beam.
- Annuli: 0-0.25, 0.25-0.5, 0.5-1.0, 1.0-2.0, 2.0-3.0" (plus 0.5-1.5" for the cross-check), centred on the Xb4b 0.26" HCO+ peak (N 156.96346, -43.90388; S 156.96339, -43.90533, deg), identical centres for Xb4d. Annulus value = mean T_b over pixels; error = rms per-pixel line-free noise (includes spatial correlation) / sqrt(N_ind), N_ind = pixels / beam area. Statistical only; 14% flux scale on absolute ratios.

## Xb4b, 0.26" circular beam (matched, essentially native 0.25")
| annulus (") | N R | S R | N HCN / HCO+ (K km/s) | S HCN / HCO+ (K km/s) |
|---|---|---|---|---|
| 0-0.25 | 0.998±0.019 | 0.443±0.014 | 235.8±2.6 / 236.3±3.7 | 94.1±2.6 / 212.3±3.5 |
| 0.25-0.5 | 0.877±0.017 | 0.386±0.017 | 118.4±1.5 / 134.9±2.1 | 37.2±1.5 / 96.3±2.0 |
| 0.5-1 | 0.648±0.011 | 0.331±0.021 | 57.2±0.7 / 88.3±1.0 | 12.7±0.7 / 38.2±1.0 |
| 1-2 | 0.432±0.017 | 0.340±0.033 | 11.0±0.4 / 25.4±0.5 | 4.3±0.4 / 12.5±0.5 |
| 2-3 | 0.377±0.054 | 0.446±0.027 | 2.4±0.3 / 6.2±0.4 | 5.7±0.3 / 12.9±0.4 |
| 0.5-1.5 | 0.573±0.010 | 0.318±0.020 | 32.4 / 56.5 | 7.8 / 24.4 |

- **N**: 1.00 (<0.25", r<49 pc) → 0.38 at 2-3": confirmed (earlier native-beam s39 numbers 0.97 → 0.39; difference only the beam match). **S**: flat, 0.44 → 0.33-0.34 over 0-2" (mild decline, index small), consistent with ~0.4. The S 2-3" annulus (0.45) rises: it is the outer, low-brightness region where other emission (bridge/extended gas; HCO+ T_b is flat at ~12 K km/s) dominates, and where Xb4b loses extended flux; it is not part of the compact S profile. The N 2-3" value is also affected by missing flux (see below) and is noisy.
- Contrast N/S in the same annuli: 2.25 (0-0.25"), 2.27, 1.96, 1.27, 0.85 — the contrast exists only inside ≈1" (≈190 pc).

## Convolved to 1.0" (both arrays, same beam)
| annulus | N Xb4b | N Xb4d | S Xb4b | S Xb4d |
|---|---|---|---|---|
| 0-0.25 | 0.842±0.008 | 0.844±0.007 | 0.394±0.007 | 0.406±0.006 |
| 0.25-0.5 | 0.807±0.009 | 0.812±0.007 | 0.386±0.009 | 0.399±0.007 |
| 0.5-1 | 0.702±0.008 | 0.721±0.007 | 0.362±0.010 | 0.381±0.008 |
| 1-2 | 0.503±0.009 | 0.535±0.008 | 0.366±0.014 | 0.379±0.010 |
| 2-3 | 0.368±0.026 | 0.395±0.017 | 0.435±0.016 | 0.377±0.012 |
| **0.5-1.5** | **0.631±0.006** | **0.656±0.006** | **0.353±0.009** | **0.376±0.007** |

## Xb4b vs Xb4d at 0.5-1.5" (1.0" beam)
- Ratio: Xb4d higher by +0.025 at N (+4%) and +0.023 at S (+6.5%); 0.5-1.0": +0.019 / +0.019; 1-2": +0.032 / +0.013. The N/S contrast in 0.5-1.5": 1.79 (Xb4b) vs 1.74 (Xb4d), a −0.05 (3%) difference, i.e. the Xb4d-Xb4b uv effect on the contrast is small, consistent with r07 (−0.06 at 1.0").
- Flux: HCO+ agrees at N (56.2 vs 56.4 K km/s) while HCN is 4% lower in Xb4b (35.5 vs 37.0). At S Xb4b is lower in both: HCO+ −16% (28.6 vs 34.1), HCN −21% (10.1 vs 12.8). Inside r<0.5" the arrays agree within 2-8% (N 0-0.25: HCN 110.6 vs 112.7; S 35.3 vs 38.9, −9%). So Xb4b loses ~15-20% of S flux at 0.5-1.5" and slightly more in HCN than HCO+ (the recovered extended component has ratio ≈0.38), which is why the S Xb4d ratio is slightly higher. The N central ratio is robust to this: 0.84 (Xb4b) vs 0.84 (Xb4d) at 1" beam.
- Caveat: Xb4b and Xb4d both lack zero spacings; beyond ~1.5-2" both lose flux (r07 item 4; the Xb4d−Xb4f comparison gives the rest), so the 1-3" annuli (and the N 0.38-0.43 outer values) are lower/upper limits in absolute terms, with ±0.03-0.05 level differences between arrays. The inner 0-1" profile is robust against array choice at the ≤0.03 level.
- Resolution caveat: at the 1.0" beam the 0-0.25 and 0.25-0.5" annuli are not independent (beam ≫ annulus); the 0.26" profile is the one that resolves the core.
