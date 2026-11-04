import json
p='/output/report.json'
d=json.load(open(p))
R=d['results']
def find(s):
    for r in R:
        if r['quantity'].startswith(s): return r
    raise KeyError(s)
B3='/snapshots/team-b-s03/'
# 5: N/S ratio
r=find('Integrated-intensity ratio HCN/HCO+ at matched J')
r['evidence']+=(" || p0005 + p0004-r01 (team-b, "+B3+"out/s30_ratios.txt, s34_native_ratios.txt): independent same-pixel native-beam measurement CONFIRMS: 0.25\" N 1.03+-0.03, S 0.48+-0.02, Bridge 0.56; 0.94x0.74\" N 0.88-0.90, S 0.39-0.40, Bridge 0.48-0.50; 6.5\" (Xb4d/Xb4f) N 0.54/0.54, S 0.39/0.42 (N/S 1.35-1.4 at 6.5\", slightly above team-a's 1.2 because of position/aperture choice; both show the contrast shrinking with beam). Standing matched-beam values unchanged; native-beam N rises 0.60->0.88->1.03 from 5.7\" to 0.25\" while S stays 0.38->0.40->0.48 (see radial-trend result). Synthesiser spot-check of s30_ratios.txt Xb4d native line agrees.")
# disk
r=find('Extended disk region definition')
r['evidence']+=(" || p0005/p0004-r01 (team-b, "+B3+"out/s36_disk.txt): independent annulus analysis at a common 10\" beam (PB>=0.5 mask, far-field zero level) finds R(2-1)=0.575 (0-4\"), 0.554 (4-8\"), 0.489+-0.004 stat (8-12\"; +-0.03 incl. bowl, up to +-0.1 for zero-level non-flatness), i.e. a 12+-5% decline, in agreement with team-a 0.47-0.50 vs 0.56-0.62. Beyond 12\" Band 5 is NO-GO (negative bowl); Band 6/7 8-12\" annuli are NO-GO. Team-b withdrew its earlier 'R=0.30 at 8-20\"' as an artefact of the negative 17-24\" ring. Two teams now agree on the 2-1 disk decline (modest, 10-20%, not a factor 2); 3-2/4-3 disk ratios are unconstrained.")
r['confidence']='low'
r['value']+=" || SYNTHESIS: standing disk result = HCN/HCO+ 2-1 at r=8-12\" is 0.49 (+-0.03 stat+zero level; systematic up to +-0.1) vs 0.55-0.58 inside r<8\" (two teams agree), a ~10-20% lower ratio; 3-2 and 4-3 disk ratios not measurable (go/no-go failed). Within r<8\" all arrays show the ratio falling outward by 10-20% (2-1 0.575->0.554; 3-2 0.44->0.40; 4-3 0.32->0.27, team-b)."
# Sblue
for key in ('Southern blue-shifted dense-gas component','Q1: Soutflow candidate'):
    r=find(key)
    r['evidence']+=(" || p0005 (team-b, "+B3+"out/s38_sblue2.txt) CALIBRATED UPDATE: Xb4d wing -400..-150 km/s, r=0.75\": HCN 0.449+-0.062 (7.3 sigma), HCO+ 0.320+-0.060 (5.3 sigma); R_wing=1.40 (68% 1.12-1.79, 95% 0.90-2.35), P(R>1)=0.93; R_wing>R_core(0.55+-0.13) robust (P=0.999); null 0/400 random positions with both wings >3 sigma; peak-pixel injection shows no bias. Replication weak: Xb4b convolved to Xb4d beam gives R=0.86+-0.25 (difference 1.3 sigma), Xb4f 0.83+-0.25; R rises with blueshift (0.55 core, ~1 at -200..-150, 1.3-2.6 beyond, bins 1-3 sigma). STANDING VALUE: wing HCN/HCO+ 3-2 = 1.2+-0.4 (range 1.0-1.5; combines Xb4d 1.40 and Xb4b-smoothed 0.86; team-a two-Gaussian 1.20 (1.01-1.39) agrees); the earlier box-window values 1.7-1.8 are dropped. Robust: wing ratio > core ratio (~2-3x the S-core 0.40-0.55); NOT established: wing ratio > 1. Team-a full-window 0.75 (Xb4d) is core-dominated and consistent.")
find('Southern blue-shifted dense-gas component')['confidence']='low'
r=find('Q1: Soutflow candidate')
r['value']+=" || SYNTHESIS (after p0005): standing wing ratio 1.2+-0.4 (not 1.7-1.8); R_wing>R_core robust, R_wing>1 unconfirmed."
# deblend
r=find('Q1b: Deblended HCN/HCO+ per component')
r['evidence']+=(" || p0005 (team-b, "+B3+"out/s37_deblend2.txt): independent signed-template deblend at 10\" gives 2-1 N 0.64+-0.12, S 0.46+-0.08; 4-3 (6.5\") N 0.46+-0.07, S 0.35+-0.07; S/N ratio 0.71+-0.17 (2-1), 0.75+-0.17 (4-3), template 3-2 0.67. Agrees with team-a (N 0.69/0.47, S 0.49/0.25..0.35 at 2-1/4-3); S is lower than N at every J, but injection tests give 15-25% (up to +50% at 2-1 for S) template bias, so 2-1 contrast stays unconstrained (1-3 allowed) and 4-3 contrast ~1.3-1.9 (team-a 1.9+-0.5, team-b 1.3+-0.4): treat as 1.3-1.9, consistent with the 3-2 value 2.0 at 1\".")
# ladders: wider systematic
for key in ('Q1: Ladder ratios at 10','Q1b: Deblended ladder ratios'):
    r=find(key)
    r['evidence']+=(" || p0005 (team-b, "+B3+"TEAM_REPORT cross-band systematics entry, out/s37_deblend2.txt): cross-band flux ladders cross Bands 5/6/7 and different arrays (recovery differences 10-35%); team-b recommends +-25-35% on ladder slopes instead of 14%, and finds 10\" Tb ladders T(4-3)/T(3-2) 0.41-0.57 and T(2-1)/T(3-2) 1.6-2.2, equal at N and S within errors. Adopt +-25-35% uncertainty; 'N=S ladders' conclusion unchanged (no difference detectable). 'HCN/HCO+ falls ~25% per J step' is a within-band ratio trend (R(4-3)/R(3-2)=0.74-0.80 at 6.5\"/10\") with ~10% extra array systematic.")
    r['uncertainty']='+-25-35% (cross-band flux scale and array recovery; was +-14%)'
# RADEX
r=find('Q1b: RADEX')
r['evidence']+=(" || p0005 (team-b, "+B3+"out/s33_radex_q1.txt) independent RADEX on the 10\" midpoint ladders (R=0.588,0.445,0.335): 31/56 grid cells reach chi2<4; the accepted N(HCN)/N(HCO+) varies inversely with n (A~3-4.5 at n~1e4-1e5, A~1 at 1e6, A~0.6-0.9 at n>=3e6). So 'X>=3 and n<1e5' (team-a) is NOT robust: it holds only for the low-density branch; n, T and abundance ratio are degenerate with ladder data alone. H13CN opacity or other tracers needed.")
r['value']+=" || SYNTHESIS: only the degeneracy is established: ladders are sub-thermal; n<~1e5 with HCN-rich columns OR n>=1e6 with comparable columns are both allowed (team-b grid)."
r=find('Origin of the N/S HCN/HCO+ contrast')
r['value']='qualitative, UNRESOLVED between abundance/opacity and excitation: N and S HCN and HCO+ ladders are equal within errors (no density/temperature difference detected, although a +-25-35% ladder systematic and 2-1 blending limit this), and CS 7-6 tracks HCO+ not HCN; but RADEX (two teams) cannot separate n, T and abundance ratio (team-a n<1e5/X>=3 vs team-b degenerate branch A~0.6-4.5). An abundance/opacity-like difference is favoured at ~2 sigma only. The contrast is spatially concentrated: N compact core (r<0.5\") has R~1, S is flat at 0.4 at all radii.'
r['confidence']='low'
r['evidence']+=' || p0005 (team-b '+B3+'out/s33_radex_q1.txt, s39_ratio_maps.txt) added degeneracy and radial structure.'
# new result: radial trend
R.append({
 'quantity':'Radial trend of HCN/HCO+ 3-2 about the N and S nuclei (cumulative apertures and annuli; Xb4b 0.25", Xb4d 0.94", Xb4f 5.7x3.5")',
 'target':'N and S nuclei, r<0.25-3"',
 'value':'N (Xb4b): annuli 0-0.25" 0.97+-0.02, 0.25-0.5" 0.87, 0.5-1" 0.65, 1-2" 0.41+-0.02, 2-3" 0.39+-0.06; cumulative r<1" 0.75, <2" 0.61, <3" 0.57 (Xb4d 0.86/0.75/0.63/0.58; Xb4f 0.60/0.58/0.56). S: flat 0.43 (r<0.25") / 0.40 / 0.37 / 0.35 / 0.39 at all r<3". The N rise inward is a compact (r<~0.5", ~100 pc) HCN-bright core diluted by larger beams (aggregate fluxes agree between Xb4b and Xb4d to 3-10%, so not filtering). Peaks of HCN, HCO+ coincide with continuum to <=0.07". Outflow-side region about Sblue (r<1.5"): 0.51+-0.03 vs S aperture 0.37+-0.01.',
 'uncertainty':'statistical 0.01-0.03; 5% HCN-vs-HCO+ calibration; HCN v2=1f excluded; vibrationally excited HCN and opacity not tested',
 'unit':'ratio (Jy km/s / Jy km/s, same pixel; Tb ratio = x1.0126)',
 'assumptions':[{'name':'window','value':'|v|<350 km/s about region velocity, same pixel for both lines'}],
 'confidence':'medium',
 'evidence':'p0005 (team-b) '+B3+'scripts/s39_ratio_maps.py -> out/s39_ratio_maps.txt, .png. Consistent with team-a matched-beam N 0.83 / S 0.40 at 1.05" (p0004) and its beam-dependence result; single-team radial profile, medium confidence.'})
# recipe
d['recipe'].append({'order':15,'stage':'verification','description':'Synthesis update 3: reconciled team-b p0005/p0004-r01 (native-beam ratios, radial trend, 10" disk annuli, calibrated Sblue wing, signed-template deblend, RADEX degeneracy) with team-a p0004; spot-checked values against s30_ratios.txt and s38_sblue2.txt.','software':'python json','parameters':[],'scripts':['/workspace/synth/update3.py']})
# interpretation
t=d['interpretation']
t=t.replace("Its wing HCN/HCO+ 3-2 is 1.2-1.8 depending on estimator (two-Gaussian 1.2, box windows up to 1.8; HCO+ wing only 5-7 sigma) versus 0.40 at S: HCN-enhanced, qualitatively",
 "Its wing HCN/HCO+ 3-2 is 1.2+-0.4 (Xb4d 1.4, Xb4b-smoothed 0.86+-0.25, two-Gaussian 1.2; earlier box-window 1.7-1.8 dropped; HCO+ wing only 5-7 sigma) versus 0.40-0.55 in the S core: the wing ratio exceeds the core ratio robustly (>3 sigma), but a ratio >1 is not confirmed by the independent array. Qualitatively")
t=t.replace("the enhancement factor is uncertain","the enhancement factor is uncertain (~2-3x the S core)")
t=t.replace("Because the ladders do not differ, the factor-2 N/S contrast is most plausibly an HCN/HCO+ abundance or opacity difference, not density or temperature (equal abundance disfavoured at only ~2 sigma).",
 "Because the ladders do not differ within a +-25-35% cross-band systematic, an abundance/opacity-like difference is the simplest reading of the factor-2 N/S contrast (equal abundance disfavoured at only ~2 sigma), but team-b's RADEX grid shows n, T and the abundance ratio are degenerate, so excitation is not excluded. The contrast is localized: N has a compact (r<0.5\") HCN-bright core (R 0.97 at r<0.25\", 0.4 at 2-3\") while S is flat at ~0.4 at all radii.")
t=t.replace("RADEX grid (low confidence, 1 dof) gives n(H2)<~1e5 cm-3, Tkin>~30 K at N, N(HCN)/N(HCO+)>=3 at both nuclei, with n and T degenerate.","RADEX grids (low confidence, 1 dof) give n(H2)<~1e5 cm-3, Tkin>~30 K, N(HCN)/N(HCO+)>=3 on one branch (team-a) while team-b shows an equally acceptable high-density branch (n>=1e6, abundance ratio ~1); n, T and X are degenerate.")
t=t.replace("(5) Disk. Only a candidate Band 5 HCO+/HCN 2-1 stack at r>8\" exists (HCN/HCO+ 2-1 ~0.47-0.50, 15-25% below the nuclei at 10\")",
 "(5) Disk. Only Band 5 HCO+/HCN 2-1 at r=8-12\" gives a disk ratio, now measured by both teams (0.47-0.50, i.e. 10-25% below 0.55-0.62 at r<8\"; zero-level systematic up to +-0.1); beyond 12\" and in Bands 6/7 it fails go/no-go")
d['interpretation']=t
d['caveats']+=["p0005 (team-b): the 8-12\" disk ratio carries a zero-level systematic up to +-0.1 because the far-field bowl is not flat (17-24\" ring negative); the 12% decline is therefore ~1-2 sigma if the full systematic is used.",
"Sblue wing ratio is window-dependent (-350..-100 km/s gives 1.02+-0.18 as core HCO+ leaks in); Xb4b HCN wing flux is only 53% of Xb4d (partly resolved out), so the blue-wing R relies on one array for R>1.",
"Team-b injection tests show deblended S/N ratios carry 15-25% template bias (S coefficient biased high); 2-1 N/S contrast remains unconstrained."]
json.dump(d,open(p,'w'),indent=1)
