# STORY (lead author) — NGC 3256 HCN/HCO+ at 50 pc

## QUESTION
Is an elevated HCN/HCO+ line ratio a signature of AGN heating in dusty nuclei, and what sets the HCN/HCO+ ratio in the two
nuclei of the merger NGC 3256 (N: starburst; S: low-luminosity, obscured AGN)?
Believed/debated: ALMA and single-dish surveys report HCN/HCO+ (3-2, 4-3) higher in AGN than in starbursts, the "submm
HCN enhancement" (Izumi+2016, 2016ApJ...818...42I; Imanishi+2016, 2016AJ....152..218I), but other studies find no
trend with AGN luminosity or star-formation rate, and attribute high ratios to abundance, shock heating or
compact-obscured-nucleus (CON) vibrational pumping (Nishimura+2024, 2024A&A...686A..48N; Butterworth+2022,
2022A&A...667A.131B; Aalto+2015, 2015A&A...584A..42A). For NGC 3256, Harada+2018 (2018ApJ...855...49H) saw no
enhancement at the AGN at 1.7" (~330 pc), Michiyama+2018 (2018ApJ...868...95M) gave 1-0 values. Unknown: whether the
ratio at ~50 pc, where nuclei are separable, differs between the nuclei, in which line the difference lies, and whether
excitation (J, IR pumping), opacity or abundance produce it.

## ANSWER (title)
"The AGN nucleus is the low one: HCN/HCO+ 3-2 is 2.1 +/- 0.2 times higher in the starburst nucleus of NGC 3256 than at
its AGN at 50-120 pc, carried by HCN alone and not by infrared pumping or J-dependent excitation"
Key numbers (T_b ratios, 0.26" = 50 pc): N 0.99, S 0.44, contrast 2.1 +/- 0.2 (estimator spread 1.9-2.3; 14% flux scale
cancels); HCO+ peak T_b equal (2.75 vs 2.81 K), HCN differs 2.0x; contrast same at J=1-0 and 3-2 (D = 1.00 +/- 0.05).
Cause (tentative, ~2 sigma): composition (an HCN-family abundance excess ~1.6x at N), not excitation/opacity alone.

## SIGNIFICANCE
Published NGC 3256 ratios are 1.7" (~330 pc) blends, where N and S differ by only 1.5 (Harada+2018). At 50 pc the N core
reaches ~1.0 and S stays at 0.44, i.e. the AGN nucleus sits at/below the minimum of 21 resolved U/LIRG nuclei of
Nishimura+2024 while the starburst nucleus is ordinary. A same-galaxy, same-array pair with identical HCO+ brightness
shows that a single-J HCN/HCO+ ratio can invert the expected AGN/starburst order, and that the difference is not IR
pumping (HCN-vib limits >=60x below CON level) nor J-excitation. It is the first such resolved test in a merger with a
confirmed AGN, and constrains what the ratio can be used for.

## MAIN RESULTS (order of the argument)

### R1. The N/S contrast is real, compact, scale-independent below ~120 pc, and lives in HCN
- Claim: HCN/HCO+ 3-2 is 0.99 (N) vs 0.44 (S) at 50 pc; contrast 2.1+/-0.2 at 50-120 pc, falls to 1.9 at 250 pc and
  ~1.1 at 6.5" (blending); the N ratio rises inward 0.39 (2-3") -> 0.97 (<0.25"), a ~100 pc HCN-rich core; S flat 0.37-0.43.
  HCO+ T_b equal at N and S; N HCN rises, S HCN falls; ratio velocity-flat in N (0.95-1.09) and S (0.36-0.43); HCN is not
  more extended than HCO+ (cores ~50-65 pc Gaussian / 80-100 pc half-max, marginally resolved).
- Evidence: p0031 (team-b-s19; out/q26/s171_summary.txt), p0005 (team-b-s03), p0004 (team-a), p0035 (team-a-s16: R not a function of
  continuum, equal-HCO+ R 0.63 N vs 0.43 S), p0026 (team-b-s16), p0029 (team-a-s12), p0003/p0008 (array consistency <=9%).
- Figure: ratio vs beam size/radius for N and S (with the T_b of HCN and HCO+), plus ratio maps (see requests r01, r02).

### R2. Excitation, IR pumping and shocks do not produce the contrast; it is shared by the HCN/HNC/CN family
- Claim: (a) the N/S contrast is the same at 1-0 and 3-2 (D = 1.00+/-0.05 with 12m+7m; 1.07+/-0.05 in the systematics
  budget; the earlier 1.10 and the HCN 3-2/1-0 S/N = 0.71 were uv-incompleteness artefacts, S/N = 0.93+/-0.05);
  (b) no HCN v2=1 emission (Sigma_vib < 0.02 Lsun/pc2, L_vib < 29 Lsun at N, <25 at S, >=60x below CON criterion; T_b(1.1 mm)
  0.1-0.2 K): neither nucleus is CON-like; (c) CN, HNC and HCN are all 2-3x weaker per unit dust continuum at S while HCO+ is
  equal (1.4x), CN/HCN and HNC/HCN lower at S: no XDR signature; (d) no shock tracer signature in the nuclei (CS 7-6/HCN
  4-3 uniform, SO/HCO+ <0.014).
- Evidence: p0036 (team-b-s21), p0060 (team-b-s31), p0057 (team-b-s29), p0058 (team-b-s30; team-a p0056/s32 reproduces), p0059 (team-a-s32),
  p0016.
- Figure: N/S contrast vs J (1-0, 3-2) and vs species (HCN, HNC, CN relative to HCO+ and to continuum) (r03).
- Caveats to state: per-dust factor depends on Td_N/Td_S (0.5-2); 4-3 N/S contrast unresolved (1.6 vs 3.7).

### R3. Excitation/opacity alone are disfavoured at ~2 sigma; an HCN-family abundance excess at N is favoured (tentative)
- Claim: in a joint non-LTE fit of 1-0/3-2 HCN, HCO+ and 13C isotopologues, equal-abundance models reach a contrast
  of 1.1-1.3 vs 1.45-1.50 observed (2.2" beam), dchi2 5-8 (1.7-2.1 sigma); abundance-only change fits with BF 25-44
  (2.6-2200 across error models); N(HCN)/N(HCO+) = 3.2 [2.5-4.0] at N and 2.0 [1.6-2.5] at S (~1.6x). Alternative not excluded:
  different 12C/13C (X_S>=130 vs X_N~80). H13CN/HCN = 0.022+/-0.008 at N, <0.03 at S (3 sigma): cannot decide
  (needs ~10x depth). H13CO+/HCO+ N/S 2.2+/-1.2 (equal-abundance branch predicts 1; 12C/13C branch 1.6-2). Physical conditions:
  single excitation state rejected (p<0.04); extended gas Tk >~80 K at n~1e4; compact ~0.5" (100 pc) HCN-rich core; Tk of that
  core undetermined (HCN 4-3/3-2 at <=0.5" decides). S core cool (Tc_S = 17+/-1 K) while N Tc unconstrained.
- Evidence: p0062 (team-b-s33), p0061 (team-b-s32), p0059 (team-a-s32), p0063 (team-b-s34), p0051 (team-a), p0035 item 8 (RADEX: x3-4
  abundance ratio needed to move R 0.4 -> 1).
- Figure: observed contrast vs prediction of each hypothesis (excitation-only, column-only, abundance-only, full) with
  profile likelihood (r04).

### R4. For the diagnostic: S (AGN) is the outlier, N is an ordinary LIRG nucleus
- Claim: N 0.99 (50 pc) is the 38th percentile of Nishimura+2024's 21 resolved 3-2 values (median 1.13; range 0.45-2.67);
  S 0.44 is at or below the minimum (only ESO 148-IG002 0.445 comparable). 1-0: N 0.96/S 0.645 (Harada) matches our
  0.84-0.94/0.57-0.58, N at GOALS mean. 4-3: both nuclei lie in the starburst range of Izumi+2016 (0.24-1.03; N 0.52, S 0.29
  at 800 pc; at 188 pc N 0.81, S 0.22): the diagnostic gives a false negative for the S AGN. S has the stronger AGN/free-free
  continuum but the lower ratio, disfavouring AGN/XDR-driven HCN enhancement. 6.5" blends hide the N core (both 0.44).
- Evidence: p0055 (team-b-s28), p0053 (team-a-s30), p0035, p0058.
- Figure: N and S placed on the Nishimura+2024 distribution (ratio vs resolution/physical scale) and on the Izumi diagram (r05).

## TRIAGE (post id: disposition)
Appendix: supports a main result
- p0003: Band 6 flux recovery across Xb4f/Xb4d/Xb4b, ratio agrees <=9% -> supports R1
- p0008: flux-recovery and ladder-ratio systematics <~10% -> R1/R2
- p0010, p0013: wing HCN/HCO+ rise is beam mixing; S deficit at all velocities -> R1
- p0034: 1-0 contrast persists, D=1.10 (superseded by p0036) -> R2 (history only)
- p0036: D = 1.00+/-0.05, r31 N/S ~1.0 (adopted) -> R2 (this is main-text evidence)
- p0060: R32/10 S/N 0.93+/-0.05; 0.71 retracted -> R2
- p0009: only ~40% of log contrast J-dependent (older; superseded by p0036/p0059) -> R2, cite as superseded
- p0007: N core ratio not tied to saturated lines/absorption; HCN/H13CN>8.3 -> R3
- p0011, p0012(retracted), p0015, p0021, p0028, p0030: ladder fits, Tk/n degeneracies, beam filling -> R3 (keep p0015 critique and
  p0051 correction; p0012 'S cool, >=80 K excluded' is NOT supported)
- p0022, p0027: T43/T32 N/S contrast (1.26) not robust to short-spacing recovery (0.97-1.24) -> R2/R3 (negative result)
- p0025: continuum edge ramp bias on HCN 4-3 0-10% -> R2/R3 noise/systematics
- p0037, p0039, p0040: CO(2-1)/13CO cores -> Tk(S) >= 26-40 K thermometer -> R3
- p0041, p0042, p0044, p0050, p0054: RADEX N(HCN)/N(HCO+) only 2 sigma; two-zone; cold thick branch disfavoured but prior-dependent;
  H13CN 4-3 uninformative -> R3
- p0046, p0047, p0048: continuum/dust background limits for envelope (Tbg <~3.2 K), pumping disfavoured -> R3/R2
- p0051, p0052: nested sampling Tc_S = 17+/-1 K, Tc_N unconstrained; earlier SMC values retracted -> R3
- p0056, p0058: HNC/CN/HCN family -> R2 (main text)
- p0057: HCN-vib / not CON-like -> R2 (main text)
- p0059: model comparison synthesis -> R2/R3
- p0061, p0062, p0063: isotopologue limits, joint fit, two-phase -> R3 (main text)
- p0053, p0055: literature placement -> R4 (main text)
- p0035: continuum independence, free-free H40a (low confidence) -> R1; H40a as appendix only, marginal (1.5-3.3 sigma)
- p0026, p0029: sizes/T_b at 50 pc, kinematics -> R1 (sizes), kinematics = side
- p0038: N HCN-rich vs HCO+-poor depends on CS reference; superseded by p0041 -> R3 (cite revision only)
Appendix: side finding
- p0006, p0013, p0014, p0017, p0018, p0024: S-outflow blue-lobe/red-ridge dense gas (HCN/HCO+ ~1-1.4, one array, 3-5 sigma),
  dense-mass 1e7 Msun, Mdot 15-60 Msun/yr; not established in Xb4b; absent in 1-0 -> side (real but off story; mention as
  possible HCN-rich shocked outflow gas, tentative)
- p0001, p0002, p0019, p0020, p0023: CH3OH detections, extended cold; Bridge enhancement not reproduced -> side
- p0032, p0033, p0043, p0045, p0049: disk/clump HCN/HCO+ decline with J (2-1 ~11%, 1-0 6-10%, 3-2 not measurable) -> side
- p0016: no nuclear shock-chemistry signature (also supports R2)
Cut
- p0001 blue-wing HCN/HCO+ 1.7 candidate (superseded by p0006/p0010); p0004 'declines with J' disk claim;
  p0005 deblended S/N 0.7+/-0.2 (model approximate); p0012 cold-S claim; p0031 lobe ratios at 250 pc (tentative);
  Michiyama-type outflow luminosities; H40a free-fall-time-type items; N/S kinematic mass ratios (inclination-driven);
  p0038 sign-flip statement; p0053 4-3 contrast 3.7 (unresolved).

## LIMITS
- Does not establish the cause: abundance (HCN-family excess at N) vs 12C/13C vs a mixture of (T,n) differences are not
  separated; model preference is ~2 sigma, depends on error model and the 12C/13C prior; one-zone/two-phase RADEX only.
- Opacity unresolved: H13CN/HCN at S is only a limit (<0.03); N detection 3 sigma. Would need ~10x deeper 1-0 H13CN/H13CO+ at
  ~2" (sigma ~0.003) or HN13C/HNC, 13CO/C18O in matched beams.
- Temperature of the compact N core is unconstrained; HCN 4-3/3-2 (and CH3CN/CH3OH/H2CO thermometers) at <=0.5" would decide.
- uv matching is image-plane only; 4-3 N/S contrast unresolved (1.6 vs 3.7); 1-0 contrast depends on aperture (1.3-1.5).
- Single-galaxy; AGN at S is energetically minor (low-luminosity), so the result says the ratio is not a *general* AGN
  signature, not that AGNs never enhance HCN. 14% flux scale and D=40 Mpc (1" = 194 pc) are assumptions.
- A true AGN-vs-SB calibration needs resolved ratios plus 13C species in a sample of AGN with known luminosities.
