# Work summary: Q27 round 2 (beam filling in the RADEX ladder fits)

## Done (reviewer fixes, in order)
1. Floor sensitivity: floor T_R,peak = I/(1.0645 w), w 100-250 km/s, f 0.7-1.0 (scripts/s108_q27r2.py, out/q27/s108.txt part A). The 250-pc floor is insensitive (N Tk 112-126 K, P(Tk>50)=0.84-0.86). The 49-pc floor is not: N Tk 100 -> 71 -> 35 K, P(Tk>50) 0.81 -> 0.68 -> 0.29 from (w=100,f=0.7) through (150,0.85) to (250,1.0). The 0.24" T_b is beam-limited, so the floor gives an UPPER limit on eta (<=0.01-0.02) and a lower limit on tau.
2. Columns/mass (part D): N(H2)=eta N(HCO+)/X, X 1e-9..1e-8, 10" beam, compared with 10"-smoothed Xb47 continuum read as all dust (Sigma <= 14-77 Msun/pc2, M <= 6e7-3e8 Msun) and the 3e3 Msun/pc2 CO nuclear peak. The thin eta=1 branch (Sigma 28-360) fits; the thick eta~0.01-0.02 branch (Sigma 90-3600 median, M 4e8-1.5e10) exceeds the dust bound unless X(HCO+)>~1e-8. Weak discriminator (X, kappa, Td unknown).
3. Separate eta_HCN, eta_HCO+ (flat box and |ln q|<=ln10) and X(HCN/H13CN)=30/50/100 (part B): Tk>50 K survives (P 0.88-0.96 N, S, Bridge); separate eta has lower ln Z by 1.2, so no evidence for differential filling; Tk 158-178 K.
4. p0012/p0009 cold-gas test (part C): P(Tk<=30 K) 0.02-0.07, Bayes factor cold/warm 0.03-0.11; chi2(25 K) = 3.1-3.5 vs 0.1 if eta may be 0.01 (thick), 36-41 for eta>=0.05, ~49 for eta>=0.3. A cool core is allowed only dilute and thick (eta~0.01, tau(3-2) 10-12); a cool gas filling the 10" beam is excluded. The comparison is not independent of p0012 (different aperture/deblend): Delta lnZ 2-3, not a rejection.

## Correction found
tau labels in s62/s64/s66/s104: F.tauH[:,1] is tau(2-1), not 3-2 (grid J index 0 = 1-0). The s104 H13CN limit was applied to tau(3-2). Fixed in s108 (tau32=[:,2], tau43=[:,3]); the Tk conclusion is unchanged (shared eta + H13CN: Tk 158 [79,223] vs 178 [112,251] before), but previously quoted 'tau(HCN 3-2)' values in Q3r4/Q27 are tau(2-1).

## Files
scripts/s104_q27_filling.py (round 1), scripts/s108_q27r2.py (round 2, 3 min), scripts/s109_update_report_q27r2.py (report update); out/q27/s104.txt, s108.txt, s108.json; report_backup_q27r2.json (pre-update report). /output/report.json: Q27 entry replaced (193 results).

## Literature/forum
p0012 (team-a S 1" Tkin 20-30 K), p0009 (team-a), p0011/p0026 (team-b), Sakamoto+2014 nuclear Sigma peaks as quoted in the report.

## Uncertain
Floor needs the true line width and a like-for-like peak T_b; eta here is a 10" area filling, not the 50 pc filling; X(HCO+) and dust kappa/Td (x3) unmeasured; one-zone RADEX, equal Tex in H13CN limit; the dust-based Sigma is not the CO mass; excitation vs dilution remains degenerate (Bayes factors 8-15, prior dependent).
