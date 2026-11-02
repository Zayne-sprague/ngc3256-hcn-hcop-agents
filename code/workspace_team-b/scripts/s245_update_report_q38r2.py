"""Q38 round 2: make the report self-consistent. Revises entries (disk annuli, matched-J 10", observed 10" ladders, 10" region comparison, disk contrast, disk positions) and rewrites the Q38 entry.
Inputs: out/q38/s242 (NNLS injection), s243 (RADEX E1 vs N at matched 2.2"), s244 (clump vs nucleus ladders). usage: python3 -I s245_update_report_q38r2.py"""
import json
R=json.load(open('/output/report.json')); res=R['results']
def find(p):
    for x in res:
        if x['quantity'].startswith(p): return x
    raise KeyError(p)
def rep(x,old,new):
    assert old in x['value'],(x['quantity'][:40],old[:50]); x['value']=x['value'].replace(old,new)
HDR=('[Q38r2 REVISED: the E1/E2/W "disk" numbers are 10" clump POSITION values blended with nuclear wings (R10 at E1 is 0.63 at 10" vs 0.49 at uv-complete 2.2"). '
     'Same-method ladder (de-contaminated totals): nucleus (N+S+bridge+halos) R21/R32/R43 = 0.58/0.45/0.38 (R32/R21 0.78+-0.12, R43/R21 0.66+-0.07); E1 0.35/0.24/<0.06 (R32/R21 0.69+-0.17); W 0.49/0.33/<0.33 (0.68+-0.15). '
     'Clump/nucleus double ratio of R32/R21 = 0.88+-0.25 (E1), 0.87+-0.23 (W) (0.5-0.6 sigma from 1; vs N alone 0.78+-0.23, 1.0 sigma): the 2-1 -> 3-2 decline is COMMON to nuclei and clumps. NO "faster decline of the disk" and NO disk-vs-nucleus excitation contrast is supported at 2-1/3-2. '
     'The 4-3 clump value is an upper limit (E1 R43/R21 < 0.17 vs nucleus 0.66+-0.07: formally lower, but it rests on Xb47 zero level and PB 0.5 at E1: suggestive only). See the Q38 entry.] ')
# 99
x=find('Disk (off-nuclear) HCN/HCO+ and its radial trend')
rep(x,'so R falls from J=2-1 to 3-2 faster than at the nuclei (0.52->0.30 vs 0.58->0.42 at the midpoint). Consequence: disk gas has HCN/HCO+ ~10-15% lower than the nuclear complex at 2-1 (modest, not a factor 2); consistent with Harada+2018 that nuclear positions are not special; the 3-2 disk deficit cannot be separated from missing short spacings.',
 'which is NOT a faster decline than at the nuclei: the clump 2-1 -> 3-2 step (R32/R21 0.69+-0.17 E1, 0.68+-0.15 W from de-contaminated totals) equals the nuclear one (0.78+-0.12), double ratio 0.88+-0.25. Consequence: only the 2-1 annulus value (0.49+-0.03+-0.1 at 8-12") is a disk measurement; it is 10-15% below the nuclear complex (modest, not a factor 2; the +-0.1 zero-level systematic makes it ~1 sigma); consistent with Harada+2018 that nuclear positions are not special; for 3-2 and 4-3 no disk ratio exists (NO-GO) and the 3-2 deficit cannot be separated from missing short spacings.')
x['value']=HDR+x['value']; x['confidence']='low'
# 96
x=find('Matched-J HCN/HCO+ at a common circular 10"')
rep(x,'Where measurable the clump ratios decline with J faster than the nuclear ones (E1 0.51->0.33->0.15), consistent with lower excitation of the off-nuclear gas, but Band 7 primary-beam correction at r~8-10" (PB~0.5) is not in the error budget',
 'The clump position ratios fall with J (E1 0.51->0.33->0.15) but this is NOT shown to be faster than the nuclear decline once nuclear-wing leakage is removed (R32/R21 clump/nucleus 0.88+-0.25; 4-3 an upper limit); the Band 7 primary-beam correction at r~8-10" (PB~0.5) is not in the error budget')
x['value']=HDR+x['value']
# 111
x=find('Observed HCN and HCO+ J=2-1, 3-2, 4-3 ladders at the common 10"')
rep(x,'T(2-1)/T(3-2) is 2.6 (HCN) and 2.0 (HCO+) at N, 3.8 and 2.4 in the disk; T(4-3)/T(3-2) 0.34/0.45 at N, 0.14/0.33 in the disk.','T(2-1)/T(3-2) is 2.6 (HCN) and 2.0 (HCO+) at N, 3.8 and 2.4 at the clump positions; T(4-3)/T(3-2) 0.34/0.45 at N, 0.14/0.33 at the clump positions (10" blended values; the clump-vs-N differences are not established after de-contamination, see header).')
x['value']=HDR+x['value']
# 112
x=find('Does the 10" ladder differ between regions?')
rep(x,'Disk (E1,E2,W) vs Bridge/N/S: DIFFERENT.','Disk clump positions (E1,E2,W; 10" blended) vs Bridge/N/S: DIFFERENT at 10" in the 4-3 ratio and (weakly) 3-2, but [Q38r2] the same-method de-contaminated 2-1/3-2 ladder shows no clump-vs-nucleus difference (double ratio 0.88+-0.25) and the 4-3 clump value is an upper limit, so this 10" difference is carried by 4-3 (Xb47 PB~0.5, zero level) and nuclear-wing leakage, not established as excitation.')
# 133
x=find('Disk (E1) vs centre excitation contrast')
old=x['value']; i=old.index('[Q38 REVISION')
body=old[:i]
new=('[Q38r2 DOWNGRADED] The E1-vs-centre Tk contrast below applies ONLY to the 10" blended ladder and is NOT confirmed with uv-complete matched 2.2" data. '
 'RADEX (s243, 1-0 X961 12m + 3-2 Q34 12m+7m, HCN and HCO+, K km/s at 2.2", 5% per line, shared 10%/15% band scales, Tk<=300 K): N Tk 40 K [14,158]; E1 Tk 178 K [79,251], P(Tk_E1<Tk_N)=0.17 (W: 0.11); P(n_E1<n_N)=0.51; log A = log N(HCN)/N(HCO+): N 0.70 [0.50,0.70], E1 0.30 [0.20,0.40], P(A_E1<A_N)=1.00. Bayes factors vs free (N vs E1): same-Tk&n (abundance/column free) +1.00, same-A (excitation only) -2.10, same-A&N -1.31, same-Tk +0.01, same-n +0.94, all-equal -14.2. '
 'So at matched resolution E1 differs from N in HCN/HCO+ abundance/opacity (R10 0.49 vs 0.84 is already present at J=1-0) and not by lower Tk; the 10" posterior P(Tk_E1<Tk_centre)=0.91-0.99 is withdrawn as evidence of cooler disk gas (only 4 data per position at 2.2", Tk weakly constrained: the statement is "no lower Tk at E1", not "E1 is warm"). '
 'Matched 2.2" r31=T(3-2)/T(1-0) clump/N: E1 HCN 0.80+-0.07, HCO+ 0.81+-0.03 (a 20% lower r31, not the factor 2 of 10" T43/T32); E2 0.34/0.40 and W 0.55/0.46 depend on the 12m-only 3-2 recovery (E2 HCO+ 1.9 K km/s 12m vs 4.7 comb) and are lower limits on r31 recovery. '
 'The 10" analysis follows (superseded in interpretation): ')
x['value']=new+body
x['uncertainty']='10" analysis: ratio-of-ratios sigma 0.057-0.21 (ln), P(Tk_E1<Tk_centre)=0.91-0.99 NOT supported at matched 2.2" (P=0.17); matched 2.2": stat 2-8% per datum plus 5% per line, Tk E1 [79,251], N [14,158]'
x['confidence']='low'
x['evidence']+=' ; Q38r2: scripts/s243_q38_radex_E1.py -> out/q38/s243.txt/json'
# 136
x=find('Disk positions E1, E2, W: n and Tk')
x['value']='[Q38r2: the E1 "Tk<~56 K" statement is from the 10" blended ladder (HCN 4-3 weak) and is NOT reproduced with uv-complete matched 2.2" 1-0+3-2 data (E1 Tk 178 K [79,251], P(Tk_E1<Tk_N)=0.17, see the disk-contrast entry): withdrawn as a physical temperature.] '+x['value']
x['confidence']='low'
# Q38 entry
q=find('Q38: is the disk (E1/E2/W clump) HCN/HCO+ decline')
v=q['value']; j=v.index('CONCLUSION:')
q['value']=v[:j]+('(5) NEW, same-method nucleus (s244; s58d NNLS native-beam totals): nuc_total R21 0.58+-0.06, R32 0.45+-0.03, R43 0.38; R32/R21 0.78+-0.12, R43/R21 0.66+-0.07 (N alone 0.90/0.79/0.53: 0.88+-0.15, 0.59+-0.07; 10" position N 0.78, 0.59). Clump/nucleus double ratio of R32/R21: E1 0.88+-0.25, W 0.87+-0.23 (vs N alone 0.78+-0.23, 0.77+-0.21), i.e. 0.5-1.1 sigma from 1: the 2-1->3-2 decline is common to nuclei and clumps. 4-3: E1 R43/R21<0.17 vs 0.66+-0.07 nucleus (formally lower, but E1 HCN 4-3 total is 0.00+-0.16 Jy km/s at PB~0.5, zero level moves the 10" value 0.07-0.22: suggestive only). '
 '(6) NNLS injection test (s242; simulated skies = fitted nuclear+clump model + optional unmodelled 25"/40" extended component, ACA-like uv coverage 8.7-45 m, Hogbom clean, PB, native-beam maps refitted with the same design): recovered/true R(J) for E1 and W = 1.00-1.01 (R21), 1.01-1.03 (R32), 0.92-0.93 (R43), R32/R21 recovered 1.00-1.03 for a flat truth (0.70 -> 0.70-0.73 for a true decline 0.5/0.35), nuclear R ladder exact (+-1%); E2 is biased (R32/R21 1.02-1.12; R43 recovered 0.2-0.5 of truth) and is excluded. Caveat: identical morphology in both lines, no noise, assumed uv coverage, so this validates the fit machinery, not the real uv sampling. '
 '(7) E2: 12m-only vs comb 3-2 differ by a factor 2.5 (HCO+ 1.9 vs 4.7 K km/s; R32/R10 0.18+-0.23 vs 0.63+-0.11) and the NNLS E2 is biased: E2 is DROPPED from the ladder/double-ratio claims; E1 and W 12m-only and comb agree (R32/R10 0.72 vs 0.74; 0.89 vs 0.89). '
 '(8) RADEX with matched 2.2" 1-0 + 3-2: no lower Tk at E1 (P(Tk_E1<Tk_N)=0.17), E1 differs from N by lower HCN/HCO+ abundance/column ratio (see the disk-contrast entry). '
 'CONCLUSION (CONFIRMED in part / REVISED for p0032 and p0033): p0032 is right that 8-12" 3-2 and 4-3 disk annulus ratios are not measurable, and our 0.52/0.34/0.14 were clump-position values at 10", not an annulus ladder (REVISED). Robust: HCN/HCO+ at the E1 and W clumps falls by ~30% from 2-1 (or 1-0) to 3-2 (R32/R21 0.69+-0.17 / 0.68+-0.15; R32/R10 0.74+-0.07 / 0.89+-0.06 at 2.2"), the SAME decline as in the nuclei (0.78+-0.12; R32/R10 0.75+-0.01 at N). NOT supported: a faster disk decline, or disk excitation lower than the nucleus (matched 2.2" r31 clump/N 0.80-0.81, RADEX P(Tk_E1<Tk_N)=0.17). The 4-3 clump ratio is an upper limit (<~0.2). The decline with J is therefore a property of HCN/HCO+ in all regions (a downgrade of the "disk decline with J" claim to "common decline with J", with the 1-0 disk ratio 0.49 (p0033: only 6-10% decline at 1-0 in the disk) being the real clump-vs-nucleus difference: 0.49 vs 0.84 at N).')
q['uncertainty']=q['uncertainty']+'; double ratios include 5% per line calibration on each ratio (10% on the ratio of ratios), sign-flip and model-variant spread; NNLS injection bias <=3% (R21,R32), 8% (R43) for E1/W under the assumed uv coverage'
q['evidence']+='; scripts/s242_q38_injection.py, s243_q38_radex_E1.py, s244_q38_clump_vs_nucleus.py, s245_update_report_q38r2.py'
json.dump(R,open('/output/report.json','w'),indent=1)
