# Q22 work summary (round 2, after reviewer "fix")

## What I did
Wrote `scripts/q22b_shared.py`: a JOINT N+S non-LTE fit (32 data, same likelihoods as round-1 `q22_joint.py`: 1.7" ladders, 0.24" Tint, six team-b p0028 10" intensities per nucleus) with separate cores per nucleus and ONE shared extended zone (same T,n,N,filling), and each 10" aperture = own core + 0.5x the other core (edge of aperture) + envelope. Envelope variants: sh (shared), shf (shared physics, separate fillings), sep (independent envelopes), same (single state per nucleus), shH (Tc 10-25, Te 50-150 K), none (cores only). Tempered SMC (12000 particles, 4 seeds base / 3 seeds variants, 153 runs).
Sensitivity: 10" scale sigma 30%, other-core weight 0/1, drop 10" 4-3 / 2-1 / all HCN / all HCO+, leave-one-line-out, Tbg 5 K. `scripts/q22b_post.py`: evidence tables, P(T) under log-T and linear-T priors (importance reweighting), masses/X.
Run: `scripts/q22b_jobs.txt` via `scripts/q22b_runone.sh` (xargs -P 15); report: `scripts/build_report_q22b.py` (Q22a,b,c,e,f replaced; Q22d from round 1 kept).
Outputs: `out/q22b/` (runs/*.csv|npz, q22b_summary.csv, q22b_PT.csv, q22b_mass.csv, post.log, run2.log). Previous report in out/q22b/report_before_round2.json.

## Results (Tbg 2.73 K unless stated)
- Models (chi2 / k / lnZ): cores only 133/12/-118; single state per nucleus 62.6/14/-77.3 (N 52, S 10); **shared envelope 14.1/17/-52.5 (BIC 73)**; independent envelopes 15.1/22/-60.9 (BIC 91); hypothesis shH 14.0/17/-49.9.
- Shared vs independent: dlnZ +8.5 (double counting removed). Shared vs single-state: dlnZ +25; survives 30% 10" scale, other-core weight 0/1 (+27/+18), Tbg 5 K (+26), dropping 4-3 or 2-1 (+24), leave-one-line-out (+18..+25). It falls to ~+5 if a whole species is dropped from the 10" data, so it rests on the HCN-vs-HCO+ pattern, not on any single line. It is a second (lower-n, lower-N, tau~1, filling ~0.26) zone that is needed, not a second temperature per se.
- Te(shared) = 53 [34,78] K (log prior), 56 [39,78] (linear prior); P(Te>=50)=0.55 (prior 0.41) / 0.63 (prior 0.71); Tbg 5 K: 35 [25,51], P=0.18/0.21; shf 0.19. Variants span 0.1-0.8: the hot-envelope claim is not supported as a requirement; the hypothesis (Tc 10-25, Te>=50) fits equally well (not excluded). P(Te>Tc)=0.90 for both nuclei at 2.73 K (0.6-0.8 at 5 K).
- Cores: Tc_N 16 [12,32], Tc_S 20.5 [13.5,32] K (Tbg 5 K: N 30 [19,66]); N-S difference not significant (P(Tc_N>=Tc_S)=0.41). Consistent with Q20 S core 17-25 K and with p0041's bimodal/poorly constrained N temperature and ~2 sigma abundance contrast.
- Masses (mu=2.8, X(HCN) equal everywhere, envelope = 10" circle): core fraction of HCN-traced gas within 10" S 0.009 [0.003,0.05], N 0.08 [0.02,0.39] (independent-envelope fit S 0.023, N 0.076). M_core(S)=2.3e7 (X=1e-8), 2.3e8 (1e-9); M_core(S)<=Mdyn=0.58e8 needs X>~4e-9 (shared) / ~1e-8 (independent). Core path lengths 0.1-1 pc at X=1e-8 << 50-66 pc: clumpy cores. Round-1 M_core 1.5e8 at X=1e-8 is superseded.

## Literature
Harada+2018 (ApJ 855, 49; 1-0/3-2 ratios used as in round 1), Sakamoto+2014 (ApJ 797, 90; outflow context; I could not extract a CO mass from it, so the "total gas" comparison is qualitative), team-b p0041.

## Uncertain
10" data are team-b p0028 values, not uv-matched to 1.7" (re-imaging/tapering not done); the 0.5 cross-weight and uniform shared envelope are assumptions; lnZ seed scatter 0.4-1.1 (2.9 at k=22); profile chi2 = best final particle; P(Te>=50) depends on Tbg, T prior and variant; mass numbers scale with 1/X and the envelope area; no independent CO-derived mass.
