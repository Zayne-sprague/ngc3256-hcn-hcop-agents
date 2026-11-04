# Q55 round 2 (worker): synthesis of the N/S HCN/HCO+ contrast

## What I did (answers to review fixes, in order)
1. **Matched-beam, core-scale joint test with a model comparison** (`scripts/q55e_joint.py`, `q55d_radex13.py`, `q55h_tables.py`; `out/q55/joint*.log|json`, `hyp_table.csv|md`). One-zone RADEX grids (HCN, HCO+, plus low-column grids for the 13C isotopologues), per-nucleus parameters (Tkin, n, N(HCN), X=N(HCN)/N(HCO+), 12C/13C), flat grid priors, marginal likelihood (ln Bayes factors) for 9 hypotheses (nothing differs, T, n, T+n, X, N(HCN), N(HCN)+X, T+n+N, all) over 12 data sets. Inputs: R10, R32, R43 in ONE common beam (1.0-1.1" core and 2.5"), absolute integrated Tb (beam filling <=1, this kills the thin/cold-S "T-only" solutions of the old test), 13C/12C 1-0 (2015.1.00993, 2.2") and the H13CN 4-3 limit; calibration nuisances shared by N and S. Old 6-datum test (`q55_ns_single_param.py`) mixed 2-6" beams and is withdrawn.
2. **New matched 1-0 at 1"** (`q55f_r10_core.py`, `out/q55/r10_core.csv`): HCN and HCO+ 1-0 from one cube (2016.1.00965 X94a spw27). HCN/HCO+ 1-0: 1.0" N 1.06+-0.06, S 0.73+-0.11 (N/S 1.45+-0.22); 2.5" N/S 1.44+-0.09 (00993: 1.58+-0.04). With 3-2 N/S 2.08 (1.1") and 4-3 N/S 3.7 (1.0", but 1.6 in Xb47) the "75% of the contrast is at J=1-0" becomes **50-75%**; the 4-3 contrast is unresolved (00493 HCN at S is 0.61 of Xb47; HCO+ 0.87).
3. **Normalisation** (`q55g_norm.py`, `out/q55/norm_cont.csv`): line / 1.3 mm continuum in matched beams. N/S continuum = 1.0-1.1 at 1". Per unit dust, HCN, HNC, CN are ~3x brighter at N (3.1, 3.0, 3.2 at 1.1"), HCO+ only 1.4x (3-2), 1.25 (4-3), 2.7 (1-0). So the species that deviates is the HCN/HNC/CN family; "HCO+ excess at S" (the old leading hypothesis) is withdrawn. Td_N/Td_S (0.5-2), S non-dust continuum and optical depth are caveats; no CO/C18O in cached cubes.
4. **Weights removed**; verdict is expressed as exclusions + Bayes factors. **Archive check**: H13CN/H13CO+ 1-0 exist only in 2015.1.00993 (2.2", already used in Q47) and 2016.1.00965 spw25 (0.77", loses 25-45% flux, H13CO+ blends with SiO, H13CN at band edge: not measured); 2018.1.00223 X858 spw29 covers H13CN 4-3 at 0.3-1" but CO 3-2 is 385 km/s away (a polynomial baseline attempt failed, discarded); no cached cube covers H13CO+ 3-2/4-3 or H13CN 3-2. Added predictions table and the observation scaling (Q55e).
5. **Forum**: p0057 (no HCN-vib at N) cited and consistent with our v2=1f limits; p0054-r01 (hot N) vs the old Q55b T-only (cold N) vs Q42/Q49 (unconstrained): P(Tkin_N>Tkin_S) = 0.13-0.76 over data sets (Q55d), sign undecided; the old cold-N result is withdrawn. Q51/Q52 not used.

## Key results (report entries Q55a-g in /output/report.json, 245 results)
- Identical physics at N and S: excluded, ln BF <= -46 in every data set.
- T-only: ln BF -0.9 (core, no R43) to -14; n-only -2 to -8 with 13C; N(HCN)-only -2 to -10; abundance-only (X) -4 to -8 at 1" with 13C but 0.0 at 2.5" without R43 or with Xb47 R43.
- Column+abundance (N thick tau32(HCN) 4-10, S thin-to-moderate 0.1-0.9; logX_N~0, logX_S~-0.25) is best or within 1.1 in all 12 sets; (T,n) excitation within 1.7 at the core. **Not decided** (3 data per nucleus).
- Sensitivity: the thick-N solutions sit at 12C/13C=100 (prior edge); capping at <=56 costs ~3-4 in lnZ and flips the core ranking toward T,n,N(HCN) / excitation.
- Discriminant: N/S of H13CN/HCN (1-0): 2.5 (column+abundance), 1.9 (all), 1.0 (X), 0.9 (T), 0.7 (T,n), 0.6 (n); measured 3+-3.5. Need sigma~0.003 at S, ~7x the 2015.1.00993 integration, tapered ~2" (1" would cost ~90x).
- Verdict list and diagnostic implication in Q55g (S has the LOW ratio; HCN/HCO+>1 is not required of an AGN; not an AGN-specific chemistry).

## Literature
Harada+2018 (Table 13), Michiyama+2018, Izumi+2016, Imanishi+2019 (via Q53), Sakamoto+2010, Aalto+2015 (Tvib), Jimenez-Donaire+2017 and Nakajima+2023 (context only; abstracts/passages seen in lit.py search, not compared quantitatively).

## Uncertain / open
- 4-3 N/S (1.6 vs 3.7) drives the ranking between X-only and column+abundance; needs a single-observation 4-3 pair at 1".
- 12C/13C prior (32-100, grid), one-zone single T/n for both species, 13C data at 2.2" applied to a 1" core, grid edges (T 20/150 K), P(T sign) from independent marginals.
- 00965 1-0 loses flux (absolute 1-0 values are lower limits); Td ratio and S non-dust continuum in the normalisation.
