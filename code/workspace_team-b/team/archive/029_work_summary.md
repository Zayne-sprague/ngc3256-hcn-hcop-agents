# Q17 round 2 work summary (worker, team-b)

Report: /output/report.json now has 174 results. The Q17 kinematics entry is kept and caveated. The old hypothesis and mass entries are replaced by 4 new ones (null test, hypotheses, r31, mass). Re-run `scripts/s107_update_report_q17r2.py`; it is idempotent from `out/report_before_q17r2.json`.

## Scripts (all in /workspace/scripts, run with `python3 -I`, outputs in /workspace/out/q17r2)
- s106a_q17_core_null.py -> s106a.txt/json, maps npz. Core-brightness-matched null, three baseline variants, 1" grid of r=1" apertures within 10" of the S-N midpoint.
- s106b_q17_channels_res.py -> s106b.txt/json, s106b_channel_maps.png. 40 km/s bins, ridge vs lobe, products/resolutions.
- s106c_q17_hypotheses.py -> s106c.txt/json. chi2/BIC comparison of S counter-lobe, N flow and Bridge strip on the red-wing map.
- s106d_q17_r31_matched.py -> s106d.txt/json. 1.9"-beam, 2.1"-aperture 3-2/1-0 ratio, RADEX n/Tk intervals.
- s106e_q17_mass_fig.py -> s106e.txt, s106e_wing_vs_core.png. Revised mass, Mdot and momentum.

## Results
1. **Core null.** The red wing is NOT just core leakage, but the significance is lower than 9.7/6.1 sigma.
   - HCN: ridge 0.246 Jy km/s against a core-scaled prediction of 0.027. The excess is 0.22+-0.05 (4.7 sigma_scatter). It is 0.18+-0.03 (7 sigma) if the N and S nuclear apertures are also excluded from the fit.
   - HCO+: ridge 0.196, prediction 0.04, excess 0.15+-0.03 (5.7 sigma; 5.0 sigma excluding N and S).
   - Baseline order and window choice change the ridge flux by only ~10%.
   - Mirrored blue window: the lobe excess is 0.33 HCN (8.8 sigma) and 0.21 HCO+ (9 sigma). The blue wing at the ridge is -0.036/-0.002, so the red skew is one-sided.
   - The reviewer's mean red wing of 0.07-0.11 at core>2 comes mostly from the N and S nuclear apertures. Away from them, at core 0.3-2, the red wing is 0.005-0.03.
   - Caveat: only 3-10 control apertures match the ridge core brightness (it is the brightest non-nuclear Bridge gas), so the scatter error is itself uncertain.
   - Null-corrected ridge/lobe flux ratio: 0.67 (0.55-0.78), down from 0.83 raw.
2. **Channels and resolution.**
   - HCN ridge flux in 40 km/s bins: 0.046, 0.081, 0.077, 0.069, 0.022 (160-200, 200-240, 240-280, 280-320, 320-360 km/s). Control sd is ~0.015, so only the three bins at 200-320 are >3 sigma. HCO+ is similar, with 3.7-4.6 sigma at 200-320 and nothing beyond 320.
   - The lobe is >3 sigma in 6 of 7 HCN bins and is a monotonically falling wing, whereas the ridge profile is a narrow 200-320 km/s feature. The profile shapes differ.
   - Products: Xb4d, tapered Xb4d and Xb4b (convolved to 1.05") agree within errors (HCN 0.246/0.261/0.221). The ridge is not resolved out.
   - Xb4f (5.7x3.5") blends ridge, N and S, so it gives no independent coarse-beam test.
3. **Hypotheses** (chi2/BIC on the 291-point red-wing map).
   - A single component on the S axis fits at r=3.1-3.25" (about 600-630 pc), with dBIC -96 (HCN) and -46 (HCO+).
   - A free-position N-flank component converges to the same clump, 2" S of N. It is not N-centred.
   - The counter-direction of the N blue lobe (PA 111 from N) gains nothing (dBIC +9/+15). There is no red wing on N's far side.
   - A Bridge strip is disfavoured in HCN (dBIC +33 relative to the localized clump) and is ambiguous in HCO+ (dBIC ~+1). The between-nuclei strip collapses to a ~0.1-0.2" clump.
   - Non-discriminating observables: collinearity, equal |v|, HCN/HCO+.
   - Verdict: undecided. The S counter-lobe is the most probable origin (~55-65%) but not established.
   - Pereira-Santaella+2024 supports bipolarity only within r<250 pc, so it says nothing about the ridge at 640 pc.
4. **r31 at matched beam** (1.9" beam, 2.1" aperture, M18 windows).
   - HCN: ridge 0.16+-0.05 vs lobe 0.37+-0.12 (2.0 sigma).
   - HCO+: ridge 0.13+-0.05 vs lobe 0.19+-0.04 (0.9 sigma).
   - This is not a significant difference, so the earlier 0.23 vs 0.6 is superseded.
   - RADEX gives n ~1e4-1e5 cm-3 for both. r31 does not discriminate a cooling outflow from heated Bridge gas.
5. **Mass.** Lobe-scaled and tentative.
   - Ridge dense mass in 180-360 km/s: 1.4e6 (X=1e-8) to 4.7e6 Msun (X=3e-9) at lobe excitation. The full range with the excitation factor (0.8-2.4) is 9e5 to 1.3e7 Msun.
   - Mdot at i=80: 2.6-8 (X=1e-8) to 9-26 Msun/yr (X=3e-9), scaling as tan i.
   - The "bipolar total 1.5x lobe" claim is removed.
   - Confidence: low-moderate.

## Remaining uncertainty
- Few matched-core controls, and the Bridge core/continuum there is complex.
- No coarse-beam or independent line test: the Xb4f beam is too large, and the H13CN/CS/SO windows at the ridge were not checked (reviewer's new question 4).
- The M18 1-0 flux and aperture matching is unverified, since there is no 1-0 cube here.
- JWST velocity-sign structure beyond 250 pc is unknown.
- No kinematic test separates a Bridge streaming flow from an outflow beyond the +250 km/s one-sided offset.
