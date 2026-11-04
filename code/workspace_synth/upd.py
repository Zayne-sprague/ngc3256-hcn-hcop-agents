import json
d=json.load(open('/output/report.json'))
t=json.load(open('/snapshots/team-a-s07/TEAM_REPORT.json'))
R=d['results']
new=t['results'][99:]
assert len(new)==10
for r in new:
    r['evidence']=("p0012 (team-a) "+r['evidence'].replace('scripts/','/snapshots/team-a-s07/scripts/').replace('out/','/snapshots/team-a-s07/out/')
      +"; p0011-r01 (team-a reply) cross-checks against team-b's 10\" fit")
    R.append(r)
# S/N temperature headline: mark
for r in new:
    if r['quantity'].startswith('Q2: Credible ranges of the same') or r['quantity'].startswith('Q2: Credible ranges of Tkin'):
        r['evidence']+="; STANDS as the Q2 temperature result: supersedes the physT30-prior Tkin 30-50 K at S quoted in the p0009 entries (floor artefact, floors 0/15/20 K agree)"
R[94]['evidence']+="; p0011-r01 (team-a): COMPATIBLE. Team-a's 1\" matched-beam fit gives N Tkin 50-100 K (median 80, weak) and S 20-30 K with >=80 K excluded (p0012), so team-b's centre Tk>~56 K at 10\" is dominated by the warm N; densities log n 4.5-5.5 at both nuclei agree with n<10^5.25. N-vs-S T contrast remains untestable at 10\" and only suggestive at 1\""
for i in (98,99,101):
    R[i]['evidence']+="; NOTE (p0012): the physT30 (Tkin>=30 K) prior used here was shown to pile S on the floor; the floor-free fit gives S Tkin 20-30 K. Opacity/abundance conclusions are unchanged (thick classes)"
d['interpretation']+="""

(9) Q2 non-LTE ranges, floor-free (team-a p0012, confirmed compatible with team-b p0011 by p0011-r01). At 1\" matched beams (N-S separated), one-zone RADEX with optically thick, sub-thermal HCN and HCO+ gives log n(H2) 4.5-5.5 at both nuclei (not distinguishable, agreeing with team-b's n<10^5.25 at 10\"), N Tkin ~25-100 K (pooled median 80 K, poorly fit, chi2 19-23/8, class dependent) and S Tkin ~15-50 K (median 25 K, >=80 K excluded at >3 sigma; T>~19 K one-sided 95%). The earlier S 30-40 K was an artefact of a 30 K prior floor and is withdrawn. Team-b's 10\" centre Tk>~56 K is dominated by the warm N and is consistent. The 'warm N, cool S' contrast is only suggestive (N class B fits at 29 K; allowed ranges overlap at 25-50 K); the robust difference is the excitation pattern Tex(HCN)/Tex(HCO+) 0.86 (N) vs 0.63 (S), valid only for Tbg<=5 K (continuum-based Tbg,eff 3.0-3.7 K; for Tbg>=10 K the S equal-Tex branch wins). Bridge, disk, wings and the S outflow lobe have one or two ratios only: n>=1e5 cm-3 from R32, no temperature constraint. X(HCN) is not measured; with a beam-filling HCN column the CO surface density needs X(HCN)>~1e-7 (Sigma 470 Msun/pc2) to >~1e-8 (Sigma 1e4); the earlier 'fills <1-3% of the cloud' statement is withdrawn. H13CN and CS add no constraint on the solutions."""
d['caveats']+=[
 "Q2 (p0012): one-zone RADEX, equal-weight pooling of classes A/B/F2 is a choice; N is fitted poorly (HCN 3-2/2-1 0.35 is 2.8 sigma below every class); S result conditional on Tbg<=5 K; grid has a 10 K edge and 0.25 dex steps; Tbg modelled as blackbody and an unresolved hot dust core is not excluded; no IR-pumping model; hyperfine structure not modelled.",
 "Q2 N-vs-S temperature contrast is suggestive only; Bridge, disk, wing and outflow temperatures are unconstrained (2 data vs 3 parameters). Harada+2018 1-0 inputs are literature values at 1.7\" (assumption)."]
json.dump(d,open('/output/report.json','w'),indent=1)
