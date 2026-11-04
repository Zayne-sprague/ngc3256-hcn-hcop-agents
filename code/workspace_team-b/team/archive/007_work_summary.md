# Q44 work summary (round 2)
Scripts: s260*, s261_q44_matched_fit.py (matched fits), s262_q44_analysis.py, s263_run_m32.sh (M3/M32), s264-s265 (core/ring vectors, 13C, null tests), s266_run.sh/s266 fit, s267 analysis, s268_q44_freeband_m32.py (free 40/40% band scale, N M32). Outputs: out/q44/analysis_q44*.txt/json, s264.txt, s265.txt, s260f.txt.
Report: entries Q44a-d (indices 221-224); backup report_backup_q44r2.json.
## Results
- Q44a headline rewritten: P(cold|N)=0.01-0.2 across priors (BF 2-60): prior-dependent, not excluded. X_N 66-77 is an aperture mean, above Tang/Martin 21-67 and in tension with Q40 (31-35). MC (CS-scaled) is non-matched, labelled supplement only. Two-zone case (team-a p0042, p0044-r03) acknowledged, not modelled.
- Q44b: M32 (12C 3-2 + 13C 1-0): N P(cold) 0.14-0.28, BF 1-4; S 0.4-0.65; C 0.2-0.27; free 40/40% band scale N 0.22-0.31. dchi2(cold-hot) at X=45, N: M +22 (10/15%) / +7.7 (20/25%), M32 -1.5/-1.2. M3 control: no discrimination. So the 12C 1-0 lines + band-scale prior (<=25%) carry the hot preference; 13C + 3-2 do not discriminate.
- Q44c downgraded to low: core-ring H13CO+ 1.9-2.1 sigma, H13CN 1.4; null test (39 noise apertures): N H13CO+ above all nulls, two-line P=0.026; S and Bridge flat. r<0.3" is a beam-dominated pixel.
- Q44d (new): core (r<1.1") and ring fits separately: X_N,core=67 [48,90], P(cold)=0.06-0.10 (0.3 for X<=178); ring X 77-81. Core X is not lower than the aperture mean, so Q40's 31-35 is not reproduced in a matched core aperture.
## Uncertain
Two-zone matched fit not done (reviewer note 2); 12C 3-2 at 1.1" from Xb4d not attempted; 12C 4-3 at 2.2" not possible (beam 4.1x2.7"). No new literature beyond that already cited (Harada+2018 Table 13; Tang/Martin range as quoted in the review).
