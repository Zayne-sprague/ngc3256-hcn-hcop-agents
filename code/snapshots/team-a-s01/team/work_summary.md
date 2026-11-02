# Worker summary — Q0, round 2 (reviewer fixes)

Conventions unchanged: radio velocity, LSRK, pbcor, Jy/beam km/s (K km/s via RJ with each cube's beam), D=40 Mpc. All scripts in /workspace/scripts (driver run_all.sh, section "turn 2"); outputs in /workspace/out; /output/report.json rebuilt by scripts/build_report.py (30 results, valid JSON). `python3 -I scripts/test_pipeline.py` passes (recovers HCO+/HCN/CS/SO and CH3OH).

## 1. Line identification (reviewer fix 1) — reviewer was right
* Bug found: `ident_*` used logA>-3.9 (CH3OH is -4.63) and species labels such as "SO 3Σ v=0" did not match my species list. New `ident.py`: local-velocity matching (HCO+ centre of the same OUS/region, not 2780), all species Eu<200 K, no Aij cut, LTE-ranked (ν·Sμ²·e^{-Eu/30K}), curated-species tier first. `ident_blind.py` re-assigns all 132 blind detections (out/blind_detections_reassigned.csv).
* **163.54 GHz = CH3OH (−1,J)-(0,J) E, J=1-3 blend**: 2.23±0.08 Jy/b km/s (27σ, 1.73 K km/s) in the Xb55 peak pixel. J=4,5,6,7,8 are NOT detected (<~0.3 K km/s vs 0.9 predicted for J=4 at 30 K LTE).
* **251.66 GHz = CH3OH 2(0)-1(-1) E (254.0154)**: centroid test (line251.py) against local HCO+ velocity: CH3OH needs −3.5…+10.7 km/s in the 5 regions with real detections, NS 11/2-9/2 f needs −44…−58 km/s everywhere (and its e partner is absent) → NS rejected. (At vsys the answer looks ambiguous; at S the local velocity matters.)
* Xb4f 2(0)-1(-1): S 0.72±0.07 (10.6σ), Soutflow 0.65, Bridge 0.35, N 0.23 (3σ). Xb4d: Bridge 0.172±0.016 (10.8σ), new W clump 0.113±0.020 (5.8σ), N/S <0.08/<0.05. CH3OH/HCO+(3-2) = 0.051 on Bridge and W clump, <0.013 (N), <0.017 (S): CH3OH is off-nuclear (agrees with Harada+2018 "not dominant at the nuclei").
* Matched beam 9.36×6.28": Xb4d recovers 55±8% of the Xb4f 2(0)-1(-1) flux → extended/resolved out. Q J=1-3 vs 2(0)-1(-1) vs Eu>45 K lines: LTE Trot=9 K (3σ 7-10 K) but χ²=53/6 dof; RADEX (e-ch3oh.dat) grid also χ²min≈55/5 dof → only "low excitation (<~15 K) or non-LTE"; no density/column constraint (N_u(Eu 20 K)=1.7e12 cm⁻² beam-averaged).
* Other CH3OH lines: A-cluster blended with SO 6_5-5_4 (Xb4f S 3.6σ marginal; SO 6-5 limits now carry this blend caveat); 264.36 (Eu 57 K) not detected; 262.83 at spw edge; **264.92 "9(0)-8(1)" 16σ at Xb4d S is the HCO+ red wing, not CH3OH** (removed). No CH3OH in Band 7.
* Other re-assignments: 176.95-176.98 = SO 5(4)-4(3) (178.605); HOC+ stays tentative; B7 338.9/350.65/352.4 GHz = spw-edge artefacts; 176.011 (8σ, B5) and 341.928 (5.6σ) unassigned. New `out/line_table_v2.csv` (old table − wrong rows + 57 CH3OH rows).

## 2. Regions / disk / astrometry (fix 2)
* `regions2.py`→out/peaks_astrometry.csv: peaks of HCO+, HCN, CH3OH, 1.1 mm continuum in Xb4d and Xb4b with Gaussian centroid errors, θ/(2S/N), 10%-beam floor. HCO+/HCN/continuum centroid offsets at N, S, Bridge (Xb4d) ≤0.12" (≤0.13 beam); S at 0.25": ≤0.05". Xb4b splits N into ≥5 sub-clumps; S is one compact peak. New region Wclump (−1.8,−2.4") (HCN-only + CH3OH). out/regions.csv updated (Disk added as a mask).
* **Extended disk**: Xb55 HCO+ 2-1 S/N>5 at r>8" (136 px, 196 arcsec²): HCO+ 2-1 6.3±0.3, HCN 2-1 3.2±0.1 Jy/b km/s (HCN/HCO+=0.51±0.02, matched). Not detected in Xb4f 3-2 / Xb47 4-3 / CH3OH; Xb55 exceeds beam-convolved Xb4f 3-2 by ~3 Jy/b km/s at r=10-12" but the 3-2 maps have a negative zero-level bowl, so extended-2-1 vs missing-flux is undecided.

## 3. Soutflow null test (fix 3) — `wing_null.py`
Wing windows relative to v_S=2722. 583 independent beams. HCO+ 9σ (7σ after empirical scatter 1.29), HCN 12σ; 0 peaks >8σ in the null windows (3-4 line-free windows each). Position (−0.25,−4.46)" HCO+, (−0.17,−4.51)" HCN (error 0.10" adopted), 1.9" = 370 pc from S. Wing flux 0.18±0.02 (HCO+), 0.33±0.03 (HCN) Jy/b km/s; wing/core 0.27/0.67 vs 0.011 at S; no red wing; spectrum is a separate broad component (mom1 2645/2564 vs 2740 km/s at S), NOT a PSF-sidelobe copy of S (predicted 0.008 vs 0.18). HCN/HCO+ there ≈1.8±0.3 (S: 0.40). Neighbouring lines (CH3CHO, CH3NH2, CH2NH) are weak/high-Eu; spectrum is flat-topped.

## 4. H13CN (fix 4) — `h13cn.py`
Explicit CO-wing models fitted together with an HCN-profile template. Zero-CO is rejected (χ²≫). Accepted wing families: H13CN/HCN(4-3) <0.09-0.11 (N), <0.135-0.15 (S), <0.075-0.15 (Bridge) at 3σ → HCN/H13CN >9 (N), >6.7 (S, Bridge). Conditional on the wing family (flexible wings can absorb real H13CN); N best-fit 0.015-0.045±0.015.

## 5. Report (fix 5)
Replaced the "unidentified" and "CH3OH not detected" entries; added CH3OH detections/limits, matched-beam CH3OH, ident rebuild, disk, Soutflow test, astrometry, H13CN, literature cross-check (from abstracts only) and an entry listing which lines are matched in each ratio. HOC+ still tentative.

## Still uncertain
* CH3OH excitation/column (LTE and RADEX both poor; Band-5 baseline ripple ~0.2 K km/s). The 2(0)-1(-1) Eu=20 K vs Q-branch ratio is off by 1.6×; possible unmodelled blend in the Q123 feature.
* Disk detection rests on one tracer/beam (Xb55) with beam-wing and zero-level systematics.
* Soutflow: nature (outflow vs separate cloud) not decided; no quantitative comparison with CO outflow geometry (Sakamoto+2014, Michiyama+2018).
* Literature cross-check not made against line tables; 176.011 and 341.928 GHz features unassigned.
* Sources used for literature: https://iopscience.iop.org/article/10.3847/1538-4357/aaaa70 , https://arxiv.org/pdf/2203.05851
