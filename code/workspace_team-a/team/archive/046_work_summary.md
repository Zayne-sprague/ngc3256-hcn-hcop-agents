# Q39 work summary
Question: does the Q22/Q36 envelope result (P(Te>=50), core Tc) depend on the 1-0 data or on IR pumping?
Method: scripts/q39_q22.py (copy of q36_q22.py, shared envelope 'sh' and independent 'sep'; Q39_MODE=drop1 removes the three 1-0 ratio rows, inf30 sets 1-0 sigma=30%; Q39_IR=tag loads HCN core grids with 14 um v2=1 pumping from scripts/q39_irgrid.py, grid points with N(v2=1)/N(v=0)>0.07 excluded). 16 seeds x 12000 SMC particles per configuration, 288 runs (scripts/q39_jobs.txt, q39_run.sh). Post: scripts/q39_post.py -> out/q39/q39_summary.csv, q39_post.log.
Note: option (c) (2-1+3-2+4-3 only) is identical to (a), since 1-0 enters only those three data rows.
Results (16 seeds; cmb / nominal background): P(Te>=50) 0.45+-0.26 / 0.46+-0.20 baseline; 0.37/0.40 drop 1-0; 0.36/0.38 30% 1-0 errors (sem 0.05, seed sd 0.2): no significant change. Tc_S unchanged (~21 K). Tc_N median falls 24-25 -> 14-15 K (large per-seed scatter at baseline, sd 15-19 K), i.e. N core's hot-tail solution depended on 1-0. ln Z baseline -52.1, drop1 -43.6/-44.3 (not comparable).
IR pumping: T_IR<=100 K (W 0.3-1) no effect; T_IR=150 K W>=0.3 and 200 K disfavoured (dlnZ -4..-10, chi2 +8..+19) and would lower P to 0.02-0.19.
Literature: none new used; builds on p0042/p0046 (Q22/Q36) and Q24 v2=1 limits (Q4c).
Uncertain: IR model assumptions (A=1.7/s, collision rates copied, W/T_IR free), 1-0 error model, multimodality (seed sd 0.2).
Suggested follow-ups: independent N-core temperature tracer (e.g. CH3CN/HC3N), Band 3 matched 1-0 with uv-matching.
