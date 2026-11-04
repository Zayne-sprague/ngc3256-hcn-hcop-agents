# Q57 work summary
**Question**: how much of the N/S HCN/HCO+ contrast (1.45 at 1-0, 1.50 at 3-2, 2.2" beam) can excitation/density/temperature/opacity alone produce vs an abundance-ratio change?

**Done**: joint N+S one-zone RADEX LVG evidence calculation (scripts/s353_q57_joint_bracket.py, s354_q57_bestfits.py, s355_update_report_q57.py; outputs out/q57/s353_*.txt, s354_*.txt, cube_N/S.npy). Data: Q44 matched 2.2" vector (HCN, HCO+ 1-0 and 3-2) + 13C 1-0 (Q29 s185_box, consistent with Q59), X=12C/13C marginalised (prior 40-100 and 20-200), filling factor and a shared 10% 3-2 band scale marginalised. Hypotheses Hexc / Hcol / HcolX / Habun / Hfull (see report entry).

**Results** (X 40-100 / 20-200):
- Equal-abundance fits reproduce contrast 1.09 (1-0), 1.23 (3-2) (Hexc); 1.13/1.35 (Hcol, opacity free); observed 1.45/1.50.
- Delta chi2 vs full: Hexc 7.7/7.6, Hcol 5.0/6.0, HcolX 3.2/1.2, Habun 3.7/5.2. Bayes factors vs full: Hexc 1.25/0.63, Hcol 0.59/0.29, HcolX 0.64/1.06 -> not excluded by evidence.
- Full fit: N(HCN)/N(HCO+) N 3.2-4.0, S 2.0-2.5; A_N/A_S = 1.6 [1.26-2.5], P(>1)=0.95/0.92. Agrees with Q44 (3.6 vs 2.2).
- A free 12C/13C per nucleus (X_S >~ 130 vs X_N ~80) can replace the abundance change (HcolX); that X_S is above Q40's 65-85 and the Tang+19 range.
- No hypothesis needs thick HCN; best fits are low-n (1e3-1e4), warm/hot (Tk>100) for X<=100.

**Literature**: Harada+18 (ApJ 855, 49; Tab. 13 ratios), Tang+19 (A&A 629, A6; 12C/13C), Jimenez-Donaire+17; earlier Q40/Q44/Q56/Q59 results of the team.

**Uncertain**: one-zone, equal filling and Tex for 13C; 13C data are one N detection and S H13CN consistent with 0 (needs ~10x deeper data, Q59); grid edge at log n = 3 for the best fits; A grid 0.1 dex; HNC and 4-3 not included (4-3 not matched); evidence is prior/Occam dependent, so only the max-likelihood statements (excitation-only cannot reach 1.45) are prior-robust.
