# Q55 work summary (synthesis of the cause of the N/S HCN/HCO+ contrast)

## What I did
- Collected the surviving results (Q4/Q4c v2=1f, Q47a-c H13CN/13C, Q42/Q49/Q50 Tc, Q8/Q54 ladders & chemistry, Q14m CH3OH/CS, Q53) from /output/report.json; excluded SMC-derived entries (Q12c/Q20/Q22/Q40e posteriors) and Q51 (out/q51 holds only unreviewed single-seed runs; item still in progress, so no core-mass/X_crit claim).
- New computation: scripts/q55_ns_single_param.py (Q40i RADEX grids, N vs S pair, which ONE parameter must differ) -> out/q55/ns_single_param.{log,json}.
- scripts/q55_report.py appends Q55a (test table), Q55b (single-parameter RADEX), Q55c (verdict) to /output/report.json (241 results, validated structure; backup out/report.pre_q55.json).

## Results
- v2=1f: Tb ratio <0.024-0.06 (N), <0.049-0.087 (S) -> Tvib<290-430 K; strong IR pumping excluded, Tvib 200-300 K not.
- H13CN alpha<0.10 N / <0.14 S: informative only at n>=1e6 (predicted 0.05-0.12 there vs 0.02-0.03 at n<=3e5); tau10(HCN) 4-11 both nuclei.
- X_S/X_N = 0.98 (0.65-1.58), P(<=1.8)=0.90: prior-filled, abundance/opacity difference neither seen nor excluded.
- Ladders: HCN 3-2/1-0 N/S 1.2; contrast at 2.5": 1-0 1.58, 3-2 1.81, so ~75% (log) of the contrast exists at 1-0 -> J-dependent excitation <=25-30%. Note 4-3/3-2 N/S is 1.07 at 4.6x3.1" but 1.6 at 1.4": beam-dependent.
- RADEX N vs S: identical physics excluded (chi2 43 for 6 data); any single parameter suffices (chi2 1.1-3.2); X needs only +0.25 dex (x1.8).
- Tc_N unconstrained 12-95 K; Tc_S 17+-1 K; P(Tc_N<=30)=0.33/0.64 (with/without 1-0).
- Chemistry (Q54): HCN, HNC, CN, CS all low vs HCO+ at S; no XDR signature; CH3OH/H2 S/N 7.7 (4.3-18.7).
- Verdict: not IR pumping, not mainly J-dependent excitation, not HCN-specific AGN chemistry; HCO+(+CH3OH) excess/composition at S weighted ~0.45, excitation/density/T ~0.30, opacity ~0.15, pumping <=0.05, XDR chemistry ~0.05 (expert-judgement weights, not computed posteriors).
- AGN-diagnostic implication from Q53: S falls below SB values at 1-0/3-2/4-3; N SB-like at >=188 pc, non-discriminating at its 56 pc core.

## Literature
Harada+2018, Ledger+2021, Izumi+2013/2016, Aalto+2015, Sakamoto+2010, Meijerink+2007 (as already used in the cited items); forum p0050-r01, p0053-p0056.

## Uncertain / decisive
Abundance vs opacity vs temperature stays open. Needed: H13CN + H13CO+ (and C34S) 1-0 or 3-2 at N and S in uv-matched 0.3-1" beams (alpha~0.03 depth), a matched-beam thermometer at the N core, uv-matched 1-0 map. Q51/Q52 outputs, when delivered, should be added as table rows. The single-parameter test mixes ratios from different beams (consistency only).
