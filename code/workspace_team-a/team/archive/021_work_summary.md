# Work summary - Q22 (joint core + envelope fit of the 1.7" and 10" HCN/HCO+ ladders)

## What I did
- `scripts/q22_joint.py` (derived from `q20b_eta.py`, same 1.7" ladder data, Tint and nuisances) adds the six team-b p0028 10" intensities per nucleus with the Q20f error model (3% stat, aperture 10/3/5%, common per-band 10% scale shared by HCN/HCO+, analytic marginalisation; line width free 100-300 km/s) and a free envelope zone. Models: M1 core only (Gaussian-core dilution at all apertures = Q20 single phase), M1E core + extended zone with the SAME (T,n,N) and own filling, M2T shared T / different n,N (two-density clumpy variant), M2 independent T, M2H = M2 with Tc 10-25 K and Te 50-150 K (the item's hypothesis). Tempered SMC (6000 particles, 2 seeds), ln evidence + profile chi2 + BIC. Runs: `q22_run.sh` (baseline, Tbg 2.73/5 K, data a = 1.7"+Tint, ab = +10"; variants core size x1.5, free 10" envelope filling, 10" scale 15%, aperture error x2), `q22_run2.sh` (Tk floors from the p0040 CO limits: S 30/40 K, N 20 K).
- `scripts/q22_post.py`: P(T), n, N, tau, abundance/path-length cut for X(HCN)=1e-9..3e-8, H13CN(4-3) limit, core/envelope masses -> `out/q22/q22_post.csv|log`.
- `scripts/build_report_q22.py` appended Q22a-f to /output/report.json (valid against the schema; backup `out/q22/report_before_q22.json`).
- Outputs: `out/q22/runs/*.csv|npz`, `out/q22/run.log`, `run2.log`, `q22_post.*`.

## Results
- Single-phase core cannot fit the 10" ladders (S chi2 177, N 135 vs 1.5-4 for two-zone): the core carries only 2-10% of the 10" flux. Trivial but now quantified.
- S: the hot envelope is NOT required. Same-state gas (M1E: T=24 [20,31] K, eta_core 0.68, extended filling ~1%) has lnZ -25.6, BIC 28.7; M2 (Tc 18 [13,29], Te 68 [38,106] K, P(Te>=50)=0.71 vs prior 0.41) lnZ -28.5, BIC 32.0; M2H lnZ -26.8 (not excluded, Bayes factor ~0.3 vs M1E). Preference for M1E is weak-moderate (dlnZ 0.7-3.2 over variants).
- N: single state is rejected (chi2 61 vs 2.8, dlnZ +23, robust to all variants, +17..+27); extended zone warmer (Te 72 [46,105] K, P(>=50)=0.72, tau32~0.3), core Tc 21 [12,42] K poorly constrained; labels unstable at Tbg 5 K. N-S contrast in either zone not significant.
- CO floor (p0040: S >=30/40 K, N >=20 K): S single-state fit costs dchi2 +0.9/+2.1 (dlnZ -1.1/-2.6): T=34 [31,41] / 44 [41,51] K with eta_core 0.76-0.82. So the Q20 "cool 17-25 K S core" is a mild preference, not a requirement.
- X(HCN): path L_c at X=1e-8 is 2.5 (M2) to 11 (M1E) pc, i.e. below the 50-90 pc core -> clumpy core, volume filling ~0.03-0.2. X=1e-9 loses 40-60% of the posterior (dlnZ -0.5..-0.8; -1.6 with the 40 K floor): disfavoured, not excluded. H13CN: S consistent at 0.5-0.97 (R=30-100); N only 0.08/0.28/0.65, so N prefers R>=50-100 or lower tau43.
- Masses (X=1e-8): S core 1.5e8 [0.3,5.5e8] Msun (x0.3 for Gaussian core ~ Mdyn 0.58e8), envelope within 10" 2.5e9; core fraction 0.06 [0.01,0.22] (M2) / 0.18 (M1E); N 0.09 [0.02,0.45].

## Literature
Brunetti+2021 (2021MNRAS.500.4730B) via p0040/Q20d for the CO floor; Harada+2018 (2018ApJ...855...49H) 1-0 ratios as before; Michiyama+2018 (two-phase outflow) as motivation only. No new external numbers; all inputs are earlier team results (Q12b, Q5, Q20, Q32/p0040) and team-b p0028 intensities.

## Uncertain / not done
10" data are team-b's (not independent); N and S 10" apertures overlap and the envelope is not shared; 1.7" and 10" flux scales independent; uniform-sphere zones, one line width, X equal in both zones; mass fractions depend on envelope area; Tbg>2.73 (IR pumping) only via the 5 K variant; no re-imaging of the 10" ladders at matched uv coverage (a raw-data request would address that and the N/S overlap).
