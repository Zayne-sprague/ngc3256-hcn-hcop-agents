# Q28 round 2 work summary (worker, team-b)

Reviewer fixes addressed (all six); results in /output/report.json (Q28 entry replaced; 252 results).

## What I did
- **Joint N+S beam-overlap model** (scripts/s362_q28r2_lib.py, s362_q28r2_fit.py): compact core per nucleus + ONE shared extended phase E; N-S separation 5.2"; 10" beam sees the other nucleus's core with weight ~0.48 (offsets 0.1-0.7 K km/s, <5% of fluxes, iterated to convergence, 30% error). Extended parameters reported once.
- **Two-sided T_b**: 0.24" and 1.29" peak T_b are data (15%) with explicit clump filling f_C; best-solution predicted 0.24" T_b tabulated per branch.
- **Finer extended grid** (Tk 0.075 dex, n 0.125 dex, N 0.25 dex, f_E 0.1-1, tau<=3), **profile likelihoods** in Tk_E, n_E and chi2(Tk_E) at fixed n (s363_q28r2_analyse.py; out/q28r2/analysis_base.txt).
- **GOF / LOO**: parametric bootstrap (s364_q28r2_boot.py, 24 sims/nucleus, out/q28r2/boot_*.log); LOO/variants: tau<=1, larger errors, drop 2.2" 1-0, 1.29", 0.24", all 2.2", all 10" (out/q28r2/fit_*.txt, analyse_var.log).
- **Dust cap**: Sigma(10")<=77 / 14 Msun/pc2 at X(HCO+)=1e-9, 3e-9, 1e-8.
- **Branch predictions** (cold Tk_C<=30 vs >100 K): HCN 3-2/1-0 and 4-3 at 0.5", H13CN 4-3, v2=1, compact mass, n x size.

## Key results (numbers in report)
- Joint chi2 24.9 (N 17.1, S 7.8) for 28 data; one-state-two-size 51.6/40.0. Bootstrap: Dchi2(1s-2p)=34 observed vs max 10.4-10.8 if one-state were true (p<0.04). N fit mildly poor (p~0.05). ln BF shared two-phase/one-state 29.5 (variants 7.8-33). Shared-E preferred over separate-E (ln BF 6.9).
- Extended: Tk_E 79-133 K (dchi2<4; <79 K excluded), n_E 4e3-6e4; a narrow ridge, not a free Tk-n degeneracy. Needs the 10" ladder; Tk_E>~40 K in every variant that keeps it; unconstrained without 10" data.
- Compact: 0.5" (~100 pc), N/S contrast in A_C (3 vs 1). Tk_C bimodal: P(cold)=0.68 N, 0.49 S vs prior 0.30, but this is driven by the 2.2" data (prior-level without them).
- Decisive tests for Tk_C: HCN 4-3/3-2 and 4-3 T_b at <=0.5" (0.34/0.65 K cold vs 0.60/0.98 K warm). 3-2/1-0 and v2=1 do not discriminate; H13CN marginal.
- Dust: ladder fits with cap only for X(HCO+)>=~5e-9-1e-8 (dchi2 1.7 at 1e-8, 9 at 3e-9, 80 at 1e-9); Sigma<=14 rejected. Position: X(HCO+) high end, or dust bound too low/structure differs. X(HCO+) not measurable (no matched 10" H13CO+).

## Literature
Harada+2018 (ApJ 855, 49; 1.7" ratios, context only), Michiyama+2018 (ApJ 868, 95); no published two-phase analysis of this target.

## Uncertain / not done
- Compact grid coarse; theta pinned at 0.5"; 24 bootstrap sims only; 150 km/s width vs narrower 0.24" line; uniform E slab; 0.24" continuum not used for the compact mass; the one-state-two-size comparator is per nucleus with E uncoupled.
- Fix applied: round-1 grid had overflow artefacts; TR clipped to [1e-6,1e3] in the new library.
