# Q42 work summary (round 2): convergent nested-sampling rerun of the Q39 joint N+S fit

## What I did
- Killed the hung round-1 grid (16 jobs on a 4-CPU cgroup, each at 25% CPU). Rewrote the driver: `scripts/q42_dyn2.py` (dynesty 3.1, `rslice` 8 slices, bound=multi; invalid points return **-inf**, not the -5e8 plateau; likelihood = `scripts/q39_q22.py` verbatim). 4 jobs at a time, OMP_NUM_THREADS=1, nohup.
- Grids: `q42_run2.sh` (base/drop1/inf30 x cmb/nominal x 3 seeds, nlive 1500; + sep-envelope base/drop1), `q42_run3.sh` (base/drop1 CMB, nlive 5000, seeds 11,12; + profile chi2 vs Tc_N, `q42_prof.py`), `q42_run4.sh` (band-scale prior sigma 0.15/0.25, base/drop1, seed 1).
- Post-processing: `q42_post2.py` (weighted quantiles, P(Te>=50), P(Tc_N<=30), mode split by 2-comp 1D Gaussian mixture; `--validate` tests it on synthetic mixtures: recovers weights to 0.01-0.03 when separation >2 sigma and refuses unresolved ones), `q42_agg.py`, `q42_agg2.py`. Outputs in `out/q42/` (runs2, runs3, runs4_sg*, q42_summary_*.csv, q42_agg_all.csv, q42_modes_*.txt, prof_*.csv). Report entry added to /output/report.json (219 results).

## Results (joint fit, shared envelope, CMB; mean over runs, sd across seeds)
| fit | Tc_N med [16-84] | P(Tc_N<=30) | Tc_S | Te | P(Te>=50) | P(Tc_N>Tc_S) |
|---|---|---|---|---|---|---|
| baseline (5 runs) | 45+/-11 [20-92] | 0.33+/-0.12 | 17.2+/-1.2 | 60+/-7 | 0.67+/-0.10 | 0.84+/-0.06 |
| drop 1-0 (5 runs) | 22+/-4 [12-52] | 0.64+/-0.09 | 16.8+/-1.3 | 57+/-4 | 0.60+/-0.05 | 0.59+/-0.05 |
| 1-0 sigma 30% (3) | 19+/-9 [11-42] | 0.76+/-0.22 | 18.9+/-4 | 53+/-8 | 0.52+/-0.12 | 0.41 |
Nominal background: base Tc_N 52, P50 0.73; drop1 Tc_N 22, P50 0.52.
Band-scale prior sigma 0.10/0.15/0.25 (single seed): P(Te>=50) base 0.67/0.48/0.24, drop1 0.60/0.56/0.11; Tc_N drop1 22/17/59 K.
Separate envelopes (1 seed): drop1 Tc_N 12 K [11-18], P(<=30)=0.96; base 29 K.

## Answers
1. **Multimodal?** Not in a resolved sense. Tc_N posterior = pile-up at the low-T prior edge (16th pct 11-12 K, floor 10 K) + broad plateau to 150 K. The 2-Gaussian split (~12 K w 0.3-0.4 + 30-45 K w 0.6-0.7 for drop1) is a shape descriptor; weights vary by 0.1-0.2 between seeds even at nlive 5000, and logZ sd across seeds (0.5-1.8) exceeds the quoted logZerr (0.1-0.2). **Convergence to the +/-0.05 target was NOT achieved for Tc_N**; it is acceptable for Tc_S and approximately for Te/P(Te>=50) (sd 0.05-0.14).
2. **Profile chi2 (sampler-free, DE, noise ~2-3):** baseline 21-22 at 10-16 K vs 16.7-18 at 20-70 K; drop1 15.0 at 10-16 K vs 12.2-12.9 at >=50 K. Delta chi2 <5 over 10-150 K: Tc_N is unconstrained (only <~16 K mildly disfavoured). Report Tc_N as a weak constraint, not a value.
3. **Earlier SMC numbers (Tc_N 24-25 / 14-15 K; P(Te>=50) 0.45/0.37) are not reproduced** by NS: medians are ~2x higher and P50 ~0.2 higher; the SMC (8000 particles) per-seed sd of 15-19 K was real sampler noise on a broad posterior, and its medians were biased low relative to NS. What survives: dropping the 1-0 lines lowers Tc_N by about a factor 2 (45 -> 22 K) and doubles P(Tc_N<=30), but this difference is comparable to the seed spread plus the prior/envelope dependence, so p0048's "24 -> 14 K" should be restated as "Tc_N poorly constrained; 1-0 data push it to higher T".
4. **Tc_N - Tc_S sign:** not defined (P(Tc_N>Tc_S) 0.84 in baseline, 0.4-0.6 elsewhere and for wider band priors). Tc_S = 17+/-1 K is robust.
5. **Envelope:** Te median 53-63 K, P(Te>=50) 0.5-0.75 at sigma_band=0.10, but it falls to 0.1-0.25 at sigma_band=0.25: dependent on the band-scale prior, as team-b found for the cold branch (p0044-r02, p0050). Cold-N-core weight likewise moves 0.25-0.96 across prior/envelope variants (separate envelopes favour cold N at drop1).

## Literature
No published nested-sampling or Tc constraint for NGC 3256 N/S cores applies; HCN/HCO+ ladders from Harada/Michiyama were already used in the likelihood (q39_q22.py inputs). Not re-searched this round.

## Remaining uncertainty / to do
- Tc_N converged mode weights would need a sampler that handles this broad flat ridge (more nlive, or a bounded-low-T reparametrisation, or a data change that constrains Tc_N, e.g. HCN 4-3/H13CN or a CO ladder in the core).
- 10/20% prior variants, multi-seed prior/sep-envelope variants not run. I did not post to the forum or revise p0048/p0046/p0044 myself; that revision is needed (see item 3).
