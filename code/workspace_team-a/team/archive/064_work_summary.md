# Q49 work summary (round 2, reviewer fixes)
Backups: out/report_before_q49.json, out/report_before_q49b.json.

## Done
1. Forum posts read: p0042 (old: Tc_N 16 [12,32], Tc_S 20.5, P(Te>=50)=0.55/0.63, "cores 16-21 K", Q22e core fractions 0.009/0.08, title "P=0.2-0.6") ; p0044 (team-b's own; only our reply r01 quotes P=0.43+-0.20, Tc 18-31 K - r03 by us already retracts 24-25/14-15 K and cites Q42); p0046 (title/findings: Te~47 K, P=0.43+-0.20, Tc_N 18-31 K, Tc_S 18-25 K; radiation-field result itself stays valid); p0048 (Tc_N 24->14 K, P 0.45->0.37, "N-core Tc depends on 1-0", "unconstrained above ~15 K" - headline retracted).
   **Need retraction replies (editor, pointing to p0051):** p0048 (full retraction of the Tc_N 24->14 and P 0.45->0.37 headline; sampler noise, SMC Tc_N seed sd 4-19 K), p0042 (Tc_N 16 [12,32] K, P(Te>=50)=0.55/0.63, Tc_N-Tc_S p=0.41, Q22e mass fractions -> recomputed), p0046 (P=0.43+-0.20 and Tc_N 18-31 K superseded; Tbg/continuum part stands). p0044: no reply needed (team-b post; r03 already carries the correction), optionally note r01's 0.43/18-31 K is superseded. I did not post replies (editor task).
2. Dependency audit: results derived from SMC posteriors flagged "DERIVED FROM SUPERSEDED SMC FIT; not re-evaluated" with confidence low and superseded_by: 109-117 (Q12c), 123, 164, 165, 200 (Q22d), 217 (Q40e), plus earlier 201-203, 205, 210, 212.
3. Q22e (204) recomputed from Q42 dynesty samples (scripts/q49_mass.py -> out/q49_mass.log; shared env, baseline, CMB, 3 seeds): core mass fraction S 0.019 [0.004,0.083] (SMC 0.009), N 0.031 [0.012,0.096] (SMC 0.084); drop-1-0: S 0.010, N 0.067. M_env(10") 3.9-5.5e10 Msun at X=1e-9 (x2-3 above SMC). Core path S median 4-22 pc (X=1e-9), N 19-30 pc; seeds scatter x2-5. SMC N>S core-fraction contrast is not reproduced. Conclusions on clumpiness / X>~4e-9 hold qualitatively but P(M_core<=Mdyn) not recomputed.
4. Result 219 (Q49) wording: "1-0 data push to higher T" replaced by "may favour, not significant": 45 vs 22 K shift +23+/-12 K vs seed sd 11 K, separate-envelope drop1 12 K, profile dchi2<5 vs noise 2-3. P(Te>=50): 0.5-0.7 at sigma 0.10 is a spread over variants/priors; 0.11-0.24 at 0.25 single seed each (seed sd 0.05-0.14). Same replacement applied in all older flagged entries.

## Literature
No new literature needed (Harada+2018 ApJ 855, 49 provides the 1.7" ratios used by the fits; unchanged).

## Uncertain
P(M_core<=Mdyn) and the H13CN test were not redone with NS samples; items 109-117 etc. are flagged, not recomputed (an NS refit of the Q12c model would be needed). N-only/S-only fits and multi-seed sigma 0.10/0.20 remain open (as per review notes).
