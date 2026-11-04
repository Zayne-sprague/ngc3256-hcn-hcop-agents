# Work summary, Q0 round 4 (reviewer fixes)

## 1. Re-identification of the 253.85-254.15 GHz feature (report #70)
Script `scripts/s20_ch3oh254.py` (Xb4f and Xb4d; out/s20_254_centroids_*.json). Method: template-shift fit of the same-pixel HCO+ 3-2 profile to the baseline-subtracted feature (linear baseline from line-free channels, fixed); shift error from dchi2<=1 with lag-1 inflation; offset = implied candidate velocity minus HCO+ 3-2 centroid.
- If CH3OH 2(0)-1(-1)E 254.01538: Xb4f S +16+-5 km/s, Bridge +10+-7, Xb4d Bridge +7+-2. If NS 253.968: -38+-5, -44+-7, -47+-2 (excluded 6-26 sigma statistically; 3.5-4.7 sigma with a 10 km/s line-to-line floor). SO2 -53/-59 and SiS +118 are excluded.
- Detection: Xb4f S dchi2 68, Bridge 56; Xb4d Bridge 98.5 (9.9 sigma); N not detected in either beam; Xb4d S not detected.
- #70 rewritten: identified as CH3OH, confidence medium-high; NS dilution and e/f discussion withdrawn. The reviewer's crude centroid at S was affected by the free window and baseline, as they said; the fixed-baseline template test agrees with their conclusion.
- Line list: `out/linelist.csv` (s0_linelist.py) already contained 254.015377 and every CH3OH line with Eu<100 K in the Band 6 windows (list printed by `scripts/s23_other_ch3oh.py`). The omission was in the round-3 blend discussion, not in the list. The sub-100 K E-type lines other than 254.015 and 266.838 have log A of -7.8 to -8.9, so they are undetectable.

## 2. CH3OH constraints and reconciliation with team-a (p0001)
- `scripts/s21_ch3oh254_maps.py` (out/s21_ch3oh254.json, mom0 npz): fixed-window moment-0 (2783+-150 km/s), annulus-MAD noise. Xb4d CH3OH/HCO+(3-2) in K km/s: Bridge 0.042+-0.012; W 0.073+-0.034 (1.9 sigma); E2 0.056+-0.030; N -0.004+-0.007 (3 sigma UL 0.02); S 0.003+-0.008 (UL 0.03). This agrees with team-a (Bridge/W 0.051, N <0.013, S <0.017) within 1 sigma at the Bridge. Our N and S limits are weaker. The Xb4f S-pixel detection (5.5 sigma, map peak 2" W of S, 7 sigma) is the 6" beam picking up off-nuclear gas, because Xb4d S is <0.12 Jy/beam km/s (3 sigma).
- `scripts/s22_ch3oh_lte_ratio.py` (out/s22_ch3oh_lte_ratio.json): Xb4f smoothed to the Xb55 9.36x6.28" beam gives S 0.90+-0.14, Bridge 0.64+-0.14, N 0.14+-0.13 Jy/beam km/s. The K km/s ratio to the J=1-3 blend is S 0.27+-0.09 and Bridge/mid-point 0.23+-0.09, which in thin LTE means Tex ~10 K (8-21 K) at S and ~14 K (9-44 K) at the Bridge. This is consistent with Tex<=12 K and excludes Trot=5-6 K. Our blend values N 0.97, S 1.38 (team-a 2.23+-0.08, re-extraction 1.84) agree within the 30% template systematic; the two pixel values overlap in a 9" beam, so summing them double counts.
- What each data set can test: off-nuclear CH3OH can only be tested with the 254 GHz line at 1" and 6". The 9" blend cannot do it. The "Q-branch 1.6x above any model" claim cannot be tested without team-a's model; we give the blend numbers only.

## 3. Other CH3OH lines and the SO 251.826 problem
`scripts/s23_other_ch3oh.py`, `scripts/s24_so_251826.py`. 266.838 GHz (Eu 57 K) is not detected (predicted 0.4-1 mJy/beam at Tex 8-10 K, below the noise), so it gives no constraint beyond Tex<~25 K. The A-type 3(3)-6(3) ladder at 251.74-251.92 GHz overlaps SO 6(5)-5(4) 251.826 (report #69). With Tex 8-21 K from the 254/blend ratio, the ladder could supply 0-100% of the S feature (35% at 10 K, all of it at 15 K). Centroid: SO -28+-12 km/s from HCO+; the ladder is +23 to +47. SO at 251.826 is therefore an upper limit (<=0.28+-0.06 Jy/beam km/s), confidence low. #69 is updated.

## 4. Compact table (#87)
A sentence was added. The table gives the peak pixel at the native beam, including leakage from neighbouring regions. #61-#64 are deblended per-region amplitudes. Peak-pixel/unified ratio is 1.10-1.30 (e.g. HCN 2-1 N 26.9/22.7 = 1.19).

## Report edits
`scripts/s25_update_report.py` (backup out/report_before_s25.json) changed #5, #69, #70 and #87, and the caveats and interpretation text. It added #88-#90 (CH3OH 254 vs HCO+ and location; Tex ratio; other CH3OH lines). The report now has 91 results.

## Uncertain
- The identification rests on one transition plus velocity. The second CH3OH line (266.838) is not constraining.
- The HCO+ 3-2 reference has an unknown line-to-line velocity difference of ~10 km/s (a floor was assumed). The N+S-blended HCO+ at Xb4f is a mixture, so the CH3OH-HCO+ offsets of +7 to +16 km/s are not interpreted.
- The Tex ratio assumes thin, matched filling, the same emitting gas in the 254 and 165 GHz lines, and g=2J+1. It carries a 30% blend-template systematic. N is formally below any LTE ratio (1 sigma).
- Xb4d recovers only ~1/3 of the Xb4f flux at the Bridge (0.13 vs 0.42 Jy/beam km/s, in beams differing 10x in area). This is consistent with extended CH3OH, but flux recovery was not measured.
- The W/E2 CH3OH signals are only 1.7-1.9 sigma and are not claimed.
