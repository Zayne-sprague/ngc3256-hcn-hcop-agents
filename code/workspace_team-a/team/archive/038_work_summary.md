# Q36 work summary: radiation field / Tbg for the Q22 fit

## Result
Question: does the continuum-implied radiation field favour Tbg 2.73 K or >=5 K in Q22, i.e. hot (>=50 K) or warm (~35 K) envelope?
**Answer: neither extreme. A uniform Tbg=5 K is not physical for the envelope; the realistic field is core-only and frequency dependent, and leaves Q22's envelope at Te ~47 K, P(Te>=50)=0.43 +-0.20 (same as CMB-only, 0.43).**
- Continuum (archive images): Tb_RJ at matched 10.5" is only 5-8 mK (N, S) in Bands 5/6/7, <=10% of J_CMB at 258 GHz. Compact cores: S 0.59 K at 0.24" (intrinsic ~0.80 K, deconvolved 0.31x0.26", 3.5 mJy); N 0.18 K (resolved 1.6x1.1", mean 0.11 K). Agrees with Sakamoto+2014 (S 860um Tb 0.24 K at 0.43").
- Dust/ff split (10.5", 3 bands): ~50% dust at 258 GHz, ff dominated <~200 GHz; Michiyama+2020 give ff 76% (N)/90% (S) at 99 GHz (their Fig. 7, 2.5-3.5" apertures). Global FIR Td 29-34 K from their 100/160(/70) um flux densities (S250 printed equal to S160, unused).
- Equivalent Tbg (nominal, per line 1-0/2-1/3-2/4-3): S core 5.7/4.0/4.3/5.2 K, N core 3.4/3.1/3.2/3.7 K; bracket (dust-only ... x2 clumping) S 3-9 K at 3-2/4-3. Envelope 2.7-3.6 K even with 10x sub-beam concentration.
- Q22 rerun (6 seeds x 12000 particles per case): P(Te>=50) 0.42 nominal / 0.35 dust-only / 0.56 high / 0.43 CMB-only; uniform 5 K 0.13; envelope-only 5 K (cores CMB) 0.11 -> the 0.55->0.18 flip of the Q22 reviewer is entirely an envelope-background effect, excluded by the 10" continuum. Cores stay cool (Tc_N 18-31, Tc_S 18-25 K). ln Z differences <=1.7 (weak).
- Seed scatter of P(Te>=50) is 0.1-0.4 for identical setups (multimodal posterior): scenario differences 0.35-0.56 are not significant.

## Method / scripts (all in /workspace/scripts, outputs /workspace/out/q36)
- q36_cont.py: matched-beam smoothing, Tb at N/S (q36_cont.json/.log). q36_sed.py: ff+dust fits, FIR Td (q36_sed.log).
- q36_bg.py: J_bg(nu) scenarios -> equivalent Tbg (q36_bg.log). q36_rt.py: RADEX-equivalent solver with frequency-dependent background (validated by q36_validate.py: median 0.02-0.1% vs team RADEX grids, 98.5-99.7% of points within 2%; fails only at N~1e18 where RADEX itself is unconverged).
- q36_grids.py: grids with the same axes as grid4_*.npz (grid_{mol}_{N,S,E}_{scen}.npz). q36_q22.py (= q22b_shared.py with per-zone grids), q36_run.sh, q36_post.py -> q36_post.log, q36_q22_summary.csv, runs/. jobs.txt: the "" (default grid) jobs were mangled by xargs and not used; q36cmb is the CMB control.
- Report entries Q36a/b/c in /output/report.json (backup out/q36/report_before_q36.json).

## Literature
Sakamoto+2014 (ApJ 797, 90, Sec 6.1); Michiyama+2020 (ApJ 895, 85: ff fractions, SED Fig. 7, FIR flux densities, Table 3 sizes); Michiyama+2018, Harada+2018 context.

## Uncertain
- Geometry factor eta (0.5-1), dust fraction in the core (0.5-1: if ff is extended the 0.3" core is dust-dominated), sub-beam clumping; ff extrapolation to 88 GHz (affects 1-0 only).
- Xb4d/Xb4f absolute scale differs ~35% (Q26): core Tb at 0.24" from Xb4b only, not cross-checked.
- No Michiyama 99 GHz nuclear flux densities were extracted (tables not parseable), so the ff share is from our 3-band fit, not from H40a directly.
- IR (14um) pumping not included (Q24), would add excitation of the cores, not the 10" envelope.
- Envelope P(Te>=50) remains prior- and seed-sensitive; decisive data would be a separate envelope temperature tracer (e.g. CH3CCH/H13CN ratios or CO high-J at 10").
