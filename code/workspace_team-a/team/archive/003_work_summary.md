# Q47 work summary (worker, team-a): verification of team-b p0050 / p0044-r02

**Verdict: CONFIRMED for (a) and (b); (c) only partly testable: our posterior shows the same qualitative thing (13C data are T-neutral, cold branch not excluded), but the specific 0.002-0.07 -> 0.2-0.58 numbers were not reproduced because they belong to team-b's one-zone Harada-intensity fit, which is not in their snapshot.**

## What was done
- (a) `scripts/h13cn.py` (Q4, rerun/reused) -> `out/h13cn_limits.csv`: H13CN 4-3 (Xb47 spw20) on the CO 3-2 wing, HCN 4-3 (spw16) template, 6 wing models, two fit windows.
  N: alpha 3-sigma upper 0.09-0.11 (exp/log wings), HCN/H13CN > 8.3-12 (team-b 0.12); zero-CO model 0.22-0.44; S 0.135-0.375 (>2.7-7.4); Bridge 0.13-0.145 (log wings). Limit is wing-model dependent; reproduces team-b.
  `scripts/q47_alpha43.py`: our RADEX grids give alpha(4-3)=0.02-0.03 at n<=3e5 (team-b: 0.016-0.035) but 0.05-0.07 at n=1e6 and 0.10-0.12 for X=32, n>=1e6, Tk>=30 K: the datum is uninformative below 1e6 but starts to bite at the highest densities (small revision of "does not constrain").
- (b) `scripts/q47_ratios10.py` -> `out/q47/ratios10.json`: from archive 2015.1.00993.S 1-0 cubes (X961 spw0 for HCN/HCO+, X95d spw0+spw1 for 13C), every channel convolved to a circular 2.2", r=2.2" apertures, estimators T (template, free linear baseline), A (+-150 km/s), P (peak pixel); errors = max(formal, scatter in 40 empty off-source apertures).
  H13CN/HCN: N 0.0217+-0.008 (T) / 0.0204+-0.0058 (A); S 0.007+-0.007 / 0.006+-0.009; C 0.026+-0.010 / 0.021+-0.006.
  H13CO+/HCO+: N 0.0192+-0.0039 / 0.0276+-0.0049; S 0.0114+-0.0092 / 0.0144+-0.0054; C 0.0111+-0.0060 / 0.0152+-0.0049.
  A and T agree with team-b (e.g. N HCO+ 0.0276 vs 0.0271, S 0.0144 vs 0.0150, N HCN 0.0204 vs 0.0207). Our peak-pixel S values (0.022, 0.033) disagree with theirs (-0.0006, 0.014): unstable pixel choice at S, not used. Our errors are 1.5-2x larger. N-S: 1.8 sigma (H13CO+, A), 1.3-1.4 sigma (H13CN). So X_S/X_N > 1 is a <=2 sigma hint; S H13CN/HCN is consistent with 0, so X_S/X_N has no upper bound in our errors (team-b's "<=1.8" is not an upper bound).
- (c) `scripts/q47_grid13.py` (RADEX H13CN 1-0 grid, `out/q47/grid_h13cn_10.npz`) and `scripts/q47_reweight.py` (-> `out/q47/reweight_1_0.csv`): importance-reweight our four N-only nested-sampling runs (q50, own envelope, CMB grids, band-scale 0.10) with the 1-0 ratios; X ~ LN(45, 0.2 dex)<=100 or logU(30-100); 10% model error.
  Before: P(Tk<=30)=0.47 (0.28-0.55). After: HCO+ only: mass kept 18-31%, P=0.27-0.32; HCN only: 2-6%, P=0.35-0.51; both: 1-3%, P=0.39-0.60 (T50 26-51 K). Low surviving mass = our thick-HCN solutions over-predict H13CN/HCN (as team-b notes, ~2x). T is not constrained (16-84% >~12 to >~90 K). The cold branch is not excluded and not favoured; strong exclusion depends on band-scale priors and on mixing 1.7" 12C with 2.2" 13C ratios, which supports team-b's softened wording.

## Literature
Harada et al. 2018 (1.7" H13CO+/HCO+ peak 0.0458+-0.0049 at N, via team-b p0050 and our out/q4c/harada_table13_rows.csv; not re-read this turn); 12C/13C 21-67 (Tang+2019, Martin+2019) as quoted by team-b, not independently checked.

## Uncertain / remaining
- (c) numbers of team-b not rerun (their s235 machinery absent); our post-hoc reweighting treats the core component only and a one-zone ratio.
- 2.2" circular kernel approximated (axes differ <15%); H13CO+ lies at the spw edge (two spws averaged); SiO 2-1 contaminates H13CO+ outside +-150 km/s (A is the clean estimator); no H13CO+ LAMDA file (HCO+ file at N/X).
- A genuinely matched fit (all 12C and 13C at one 2.2" aperture) still missing.
- Report: three new results Q47 (a),(b),(c), one recipe step and one caveat in /output/report.json.
