# STORY (lead author, revision 3 after referee report; story_ok was false) — NGC 3256 HCN/HCO+ at 50 pc

## QUESTION
Is an elevated HCN/HCO+ J=3-2 ratio a signature of AGN heating in a dusty nucleus? Tested in NGC 3256 (N: starburst; S: obscured,
low-luminosity AGN) at 50-250 pc, where two nuclei of one galaxy can be compared with the same array and tuning.
Debated: submm "HCN enhancement" in AGN (Izumi+2016; Imanishi+2016) vs. no trend with AGN luminosity and shock, abundance or
vibrational-pumping alternatives (Nishimura+2024; Butterworth+2022; Aalto+2015). Harada+2018 give N/S of about 1.5 (1-0) and 1.9 (3-2)
at 1.7" (330 pc) (their Table 13, as transcribed by the team; to be verified, q10).

## ANSWER (title: a description, not a cause)
"The HCN/HCO+ J=3-2 ratio rises inward at the starburst nucleus of NGC 3256 but stays low at its low-luminosity AGN."
N rises from 0.38 (2-3") to 1.00 (<0.25", 49 pc) over ~200 pc while S stays at 0.33-0.44 (0.45 at 2-3", where extended gas
dominates); Xb4b and Xb4d agree within 0.03. The N/S contrast is ~2 at <=120 pc in 12-m-only data (2.1, an UPPER estimate: missing
short spacings are untested below 1" and can only lower it), 1.44 flux-complete at 430 pc (same as the published 1.5). No CON-level
vibrational HCN at either nucleus (beam-filling Sigma_vib < 0.02 Lsun/pc2). Cause (abundance, opacity, 12C/13C, excitation) NOT
separated. So in this galaxy the ratio does not flag the AGN (one galaxy, energetically minor AGN).

## SIGNIFICANCE
New relative to Harada+2018: the resolved profile (N rising inward, S flat) and the CON-level pumping limit. The contrast itself at
1.7" is already published (1.5 at 1-0; ~1.9 at 3-2); our flux-complete 1.44-1.5 reproduces it, our 12-m-only 2.1 at <=120 pc is
unconfirmed for flux. Within-galaxy counter-example to a general AGN signature, consistent with no trend in Nishimura+2024 (definitions
differ; comparison indicative only). Not a statement about all AGN.

## MAIN RESULTS (one claim each)
R1 (Sec 3.1): the N ratio rises inward (0.38 -> 1.00 over ~200 pc), S stays low (0.33-0.44, 0.45 at 2-3"); arrays agree to 0.03.
R2 (Sec 3.2): the N/S contrast is ~2 at <=120 pc (12-m-only, upper estimate), 1.44 flux-complete at 430 pc; HCN N/S 2.0-2.5, HCO+ N/S
   0.72-1.37 by estimator (no "mostly HCN").
R3 (Sec 3.3): HCN, HNC and CN are all weak relative to HCO+ at S in flux-complete 1-0 (1.6-2.1 at 2.5").
R4 (Sec 3.4): no CON-level vibrational HCN (alpha < 0.05 N, < 0.10 S at 48 pc); excludes the CON range only at N, marginally.
R5 (Sec 3.5): cause not separated: fit does not reproduce the contrast; thin-line contrast H13CN/H13CO+ N/S ~2.5 (0.7-7 at 1 sigma)
   is a hint of opacity at most.
Discussion: placement among Nishimura+2024 / Izumi+2016 (indicative; matched aperture), Ledger+2021 comparison.
Moved to appendices: J-ladder, 12-m-only species, fit grids, flux recovery. Deleted: Appendix B side findings, unsupported claims.

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
- p0022, p0027: T43/T32 N/S contrast (1.26) not robust to short-spacing recovery (0.97-1.24) -> R2/R3 (negative result; MOVED to main text as a LIMIT on J-independence, with p0053)
- p0025: continuum edge ramp bias on HCN 4-3 0-10% -> R2/R3 noise/systematics
- p0037, p0039, p0040: CO(2-1)/13CO cores -> Tk(S) >= 26-40 K thermometer -> R3
- p0041, p0042, p0044, p0050, p0054: RADEX N(HCN)/N(HCO+) only 2 sigma; two-zone; cold thick branch disfavoured but prior-dependent;
  H13CN 4-3 uninformative -> R3
- p0046, p0047, p0048: continuum/dust background limits for envelope (Tbg <~3.2 K), pumping disfavoured -> R3/R2
- p0051, p0052: nested sampling Tc_S = 17+/-1 K, Tc_N unconstrained; earlier SMC values retracted -> R3
- p0056, p0058: HNC/CN/HCN family -> R2 (main text)
- p0057: HCN-vib / not CON-like -> R2 (main text)
- p0059: model comparison synthesis -> R2/R3
- p0061, p0062, p0063: isotopologue limits, joint fit, two-phase -> R3 (one main-text paragraph; BF tables, error-model and prior grids, two-phase details -> appendix)
- p0053, p0055: literature placement -> R4 (main text)
- p0035: continuum independence, free-free H40a (low confidence) -> R1; H40a as appendix only, marginal (1.5-3.3 sigma)
- p0026, p0029: sizes/T_b at 50 pc, kinematics -> R1 (sizes), kinematics = side
- p0038: N HCN-rich vs HCO+-poor depends on CS reference; superseded by p0041 -> R3 (cite revision only)
Cut (off-question; Appendix B deleted)
- p0006, p0013, p0014, p0017, p0018, p0024: S-outflow blue-lobe/red-ridge dense gas (HCN/HCO+ ~1-1.4, one array, 3-5 sigma),
  dense-mass 1e7 Msun, Mdot 15-60 Msun/yr; not established in Xb4b; absent in 1-0 -> side (real but off story; mention as
  possible HCN-rich shocked outflow gas, tentative)
- p0001, p0002, p0019, p0020, p0023: CH3OH detections, extended cold; Bridge enhancement not reproduced -> side
- p0032, p0033, p0043, p0045, p0049: disk/clump HCN/HCO+ decline with J (2-1 ~11%, 1-0 6-10%, 3-2 not measurable) -> side
- p0016: no nuclear shock-chemistry signature (also supports R2)
Cut
(EDITOR: p0053 4-3 contrast 3.7 is moved OUT of 'cut' into main text R2/LIMITS as an unresolved tension with the 3-2 contrast.)
- p0001 blue-wing HCN/HCO+ 1.7 candidate (superseded by p0006/p0010); p0004 'declines with J' disk claim;
  p0005 deblended S/N 0.7+/-0.2 (model approximate); p0012 cold-S claim; p0031 lobe ratios at 250 pc (tentative);
  Michiyama-type outflow luminosities; H40a free-fall-time-type items; N/S kinematic mass ratios (inclination-driven);
  p0038 sign-flip statement; (p0053 retained, see above).

## LIMITS
- The 2.1 core contrast is 12-m only at <=0.6"; the uv-matched 12m+7m test failed three times (q01, q03, q06 unanswered with numbers);
  the Xb4b/Xb4d comparison (-0.05) is NOT a short-spacing test (at 1.5" it gives -0.02 whereas adding 7-m gives -0.2 to -0.3). The -0.3
  is a plausible allowance, not a bound.
- Pumping: excludes CON-level only; moderate pumping (Tvib 150-300 K) or a <47 pc buried source open; S limit overlaps the CON range.
- Cause not established; fit does not reproduce 1.45; Delta chi2 5-8 has no agreed sigma (dof 1: 2.2-2.8, dof 2: 1.7-2.4).
- Opacity: H13CN/H13CO+ N/S contrast 2.5 with 1 sigma range 0.7-7 (T, A) or 1.8-3.7 (P, noise-biased).
- Temperatures: T_k(S) 17 K (nested sampling, one zone) vs >=26-40 K (CO thermometers) unreconciled; N core unconstrained.
- Annulus errors statistical + 0.03-0.05 systematic; 0.26" beam cannot resolve structure inside ~50 pc.
- Single galaxy; 14% flux scale and D=40 Mpc assumptions.

## EDITOR NOTES (round 3)
- Story re-centred on the profile; contrast demoted to a flux-uncertain upper estimate; Harada comparison made quantitative.
- Unsupported claims removed (team-b gave no numbers: shock tracers, velocity flatness, observing blocks, ratio vs continuum, Falstad limit).
