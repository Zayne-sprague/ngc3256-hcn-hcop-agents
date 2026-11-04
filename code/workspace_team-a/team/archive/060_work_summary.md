# Q49 work summary: retraction/restatement of the SMC-based temperature numbers

**What I did.** Searched /output/report.json (all string fields) and agenda plan_coverage for the old values (Tc_N 24->14 K, 24-25/14-15 K, P(Te>=50) 0.45/0.37, Tc_N 16 [12,32] K, 'cool 12-32 K cores', P(Tc_N>=Tc_S)=0.41, P(Te>=50)=0.43/0.55). No new computation; the numbers to restate come from Q42 (results[218], out/q42/q42_agg_all.csv, forum p0051). No literature needed (an internal retraction).

**report.json edits** (backup out/report_before_q49.json):
- results[201, 202, 203, 205, 210, 212] (Q22a/b/c/f, Q36c, Q39): value prefixed "[SUPERSEDED/RESTATED by Q42 and Q49 ...]" with the restatement, confidence set to low, field superseded_by added. Original text kept after the prefix for traceability.
- results[212] (p0048): explicitly retracts Tc_N 24->14 K and P(Te>=50) 0.45->0.37 as sampler-noise artefacts (SMC per-seed Tc_N medians scatter 4-19 K).
- new results[219] (Q49): the restated statement and a list of what is retracted vs not retracted.
- Not edited: results 109-117 (Q12 single-nucleus fits, a different model, already 'unconstrained').

**Restated (use these only):** Tc_N poorly constrained (profile Delta chi2<5 over 10-150 K; only T<~16 K mildly disfavoured); 1-0 data push it to higher T (NS median 45 K baseline vs 22 K without 1-0, 16-84% ~12-90 K); Tc_S = 17 +/- 1 K robust; P(Te>=50 K) = 0.5-0.7 at band-scale sigma 0.10, 0.1-0.25 at sigma 0.25 (Te median 53-63 K, 16-84% 35-95 K); Tc_N-Tc_S sign undefined.
**Not retracted:** Tc_S 17-21 K, ratios, the need for an extended second zone, 'hot envelope allowed, not required', p0046 radiation-field result (envelope background ~CMB; it does not depend on the SMC temperatures, but its quoted P(Te>=50)=0.43+-0.20 is marked superseded).

**agenda.json:** plan_coverage text for Q39/p0048 rewritten (retracted); backup out/agenda_before_q49.json. The Q49/Q42 question text still quotes the old values by necessity.

**Still to do (editor):** post replies to p0048, p0046, p0044 (and p0042 if it quotes Tc_N 16 K / P=0.55) saying: old 24/14 K and 0.45/0.37 retracted, see p0051, with the restatement above. I cannot post. Residual: the p0044 report revision (r03) already exists per the notebook; I did not re-read the forum posts, so check them for the other old numbers.
**Uncertain:** Tc_N and P(Te>=50) remain prior (sigma_band) and shared-envelope dependent; N-only/S-only fits not run (Q-open); needs new data (H13CN 4-3, CO ladder) rather than more sampling.
