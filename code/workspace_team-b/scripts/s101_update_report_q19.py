"""Q19: add the Q19 result, recipe step and amend the Q16 contrast wording in /output/report.json (idempotent)."""
import json
p='/output/report.json'; d=json.load(open(p)); R=d['results']
R[:]=[x for x in R if not x['quantity'].startswith('Q19')]
d['recipe']=[x for x in d['recipe'] if not x['description'].startswith('Q19')]
for x in R:
    if x['quantity'].startswith('Q16: matched-beam 4-3/3-2') and 'Q19 amendment' not in x['value']:
        x['value']+=(' [Q19 amendment: the HCN N/S 1.26 (2.9 sigma) is the UNCORRECTED ratio. After the Q23 ramp (x0.98) and differential Xb47/Xb4f recovery it is 1.22 (68% 1.14-1.30) with the Q12 priors alone, 1.02-1.09 anchored on the continuum, 0.97 [0.87,1.09] with the team-a physical uv model; pooled 1.10 [0.93,1.25], P(>1)=0.71, i.e. consistent with unity at <1 sigma. Quote it as "N/S T43/T32 1.26 observed, 1.1+-0.15 after recovery: not established", see Q19.]')
        x['confidence']='low'
A=lambda n,v:{'name':n,'value':v}
R.append({'quantity':'Q19: N/S HCN and HCO+ T43/T32 ratio of ratios after differential Xb47/Xb4f short-spacing recovery, Q23 edge-ramp bias and a continuum anchor (check of team-a p0022)',
 'target':'N and S nuclei, Xb4f native beam 5.72x3.54" (J=3,4); continuum at the B5 beam 8.95x5.98"; HCN and HCO+',
 'value':('Observed N/S T43/T32 (Xb4f beam, s82): HCN 1.262+-0.074 (stat 0.042, 0.2" registration 0.061), HCO+ 1.164+-0.064; team-a p0022 reproduces it (1.257 / 1.112 at the Xb55 beam; 1.23 [1.10,1.35] with their MC). Q23 ramp raises HCN 4-3 by +4+-4% (N) / +6+-5% (S): N/S x0.98+-0.025. '
  'Model: observed = true x D, D=(R47/R4f)_N/(R47/R4f)_S, D>1 if S (more extended) loses more 4-3 flux. HCN N/S true, median [68%]: (no recovery) 1.24 [1.16,1.32]; A, our Q12 priors one-sided (R47/R4f 0.92-0.97 N, 0.89-0.96 S; D<=1.034) 1.22 [1.14,1.30]; '
  'B, continuum anchor: continuum (Xb47/Xb4f)_N/(..)_S at the B5 beam = 1.14 (peak), 1.21 (r<3"), 1.01 (r<6"), 0.99 (r<10") -> with equal continuum spectral index 1.09 [1.00,1.18] (peak), 1.02 [0.94,1.10] (r3"), allowing d_alpha_N-S sd 0.3: 1.09 [0.96,1.22] / 1.02 [0.90,1.15]; '
  'C, team-a physical-Bmin uv model (D=1.28+-0.12, a prior-dependent model; their unbounded-L variant leaves 1.18-1.29) 0.97 [0.87,1.09]; pools: A+B 50/50 1.17 [1.02,1.27] (P>1 0.87); A+B+C equal weights 1.10 [0.93,1.25], 95% [0.82,1.35], P(>1)=0.71. '
  'HCO+: observed 1.164, A 1.15 [1.08,1.22], B 1.02 [0.94,1.11] (peak) / 0.96 [0.89,1.04] (r3"), C 0.98 [0.90,1.07]; pool A+B+C 1.06 [0.93,1.18]. '
  'RESULT: p0022 is CONFIRMED in its reproduction (the raw ratio) and in its statement that the contrast is not robust: the differential recovery needed to remove it (D~1.1-1.3) is not constrained by our relative-recovery priors alone (they give D<=1.03, and then N/S stays 1.22, 2.5-3 sigma) but is required by the continuum (N/S of the Xb47/Xb4f continuum = 1.14-1.21 at the matched beam, which is either recovery or a real alpha difference of 0.25-0.44 consistent between the 170-258 and 258-347 GHz pairs) and by the physical uv model. '
  'Beyond 1 sigma? Uncorrected yes (2.9 sigma HCN, 2 sigma HCO+); corrected (pool) no: HCN 1.10 [0.93,1.25] contains 1 (0.6 sigma), HCO+ 1.06 [0.93,1.18]. Our tentative N>S excitation contrast (p0015) is therefore downgraded to "not established; N/S T43/T32 = 1.1+-0.15, upper bound ~1.3-1.35 (95%)". '
  'Tbg / IR pumping: the T43/T32 ratio of ratios is a data ratio and does not depend on Tbg; Tbg only enters the Tk conversion. Existing S/N fits (s73/s74, Q16): Tbg 5 K S aperture ladders Tk 20-80 K (dchi2(80)=3.1, (100)=8.9), N 65-150 K; Tbg 10 K S 20-30 K only (dchi2(80)=16) and N 30-100 K, i.e. a higher Tbg widens, not closes, an N-S Tk gap in the p0012 8-datum fit, but Tbg>5 K is not supported here (beam-averaged dust Tbg,eff at 0.24" 3.0-3.7 K at N/S, Q16) and IR pumping is weakly constrained (HCN v2=1 / HCN <0.01-0.09, Tvib<270-360 K, moderate pumping not excluded). The N/S Tk ratio from the four ladders is 1.3-1.6 [0.7-2.7], P(Tk_N>Tk_S)=0.6-0.67, unchanged: not significant. Team-a extended-gas T_N/T_S 1.2 [0.96,1.43] is also consistent with 1.'),
 'uncertainty':'ratio of ratios: stat+registration 0.074 (HCN), 0.064 (HCO+); recovery term D model-dependent (1.0-1.3); continuum anchor error 5%+alpha prior; pooled 68% HCN [0.93,1.25]; 95% [0.82,1.35]',
 'unit':'ratio of RJ brightness temperature ratios (N/S)',
 'assumptions':[A('distance','D=40 Mpc'),A('beam','Xb4f native 5.72x3.54 PA82 for lines; continuum at B5 beam 8.95x5.98 (out/s14_common_beam.csv), peak/r3/r6/r10 apertures; the continuum aperture-dependence of D (1.14-1.21 vs 0.99-1.01) is itself a recovery/size signature'),
  A('recovery priors','Q12: R47/R4f 0.92-0.97 (N), 0.89-0.96 (S) HCN, one-sided D~U[1,0.92/0.89]; this assumes S loses at least as much as N'),
  A('continuum anchor','equal continuum spectral index N and S (sd 0, 0.3, 0.5 tested); line and continuum assumed filtered alike; continuum may contain free-free at S'),
  A('team-a model','Bmin>=15 m physical uv model from p0022 (/snapshots/team-a-s09, q8c): prior-dependent; used as scenario C, not verified (MCMC tau 400-1700)'),
  A('ramp','Q23 HCN 4-3 low-edge ramp, N/S x0.98+-0.025'),A('pooling','equal weights over scenarios A, B (peak, sd 0.3), C: a judgement, as in p0022')],
 'confidence':'low',
 'evidence':'s100_q19_NS_recovery.py -> out/q19/s100.txt/json; inputs out/q16b/s82.txt (observed), out/s14_common_beam.csv (continuum), out/q16/s74.txt (Tbg fits), /snapshots/team-a-s09/out/q8c (p0022 reproduction); s101_update_report_q19.py'})
d['recipe'].append({'order':len(d['recipe'])+1,'stage':'modelling','description':'Q19: check of p0022: observed N/S T43/T32 (s82) x Q23 ramp x one-sided differential recovery (Q12 priors, continuum anchor, team-a uv model), pooled; Tbg/IR-pumping from existing s73/s74 fits; report patch (s101)',
 'software':'python3 numpy','parameters':[A('MC draws','2e5 per scenario'),A('scenario pool','A,B,C equal weights')],'scripts':['/workspace/scripts/s100_q19_NS_recovery.py','/workspace/scripts/s101_update_report_q19.py']})
d["interpretation"]+=' Q19: the N>S ladder contrast (HCN T43/T32 N/S 1.26 raw) is not robust (p0022 confirmed): after ramp and differential Xb47/Xb4f recovery it is 1.10 [0.93,1.25] (HCO+ 1.06 [0.93,1.18]), consistent with unity; no N-S excitation difference is claimed.'
json.dump(d,open(p,'w'),indent=1)
