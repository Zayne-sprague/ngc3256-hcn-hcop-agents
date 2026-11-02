# Q6 work summary (round 2, worker team-b)

Item: sizes, line widths, surface brightness of the 0.24" dense-gas structures (HCO+ / HCN 3-2, continuum); beam filling and T_b.
Round 1 (s150-s157, unchanged, see out/q6/) gave half-max sizes, widths, peak T_b, missing flux, RADEX filling. This round answers the four review fixes.
Report: /output/report.json (191 results; Q6 entries patched, 3 new "Q6r2:" entries). Run order: s157 then s164 (s157 rebuilds the round-1 entries, s164 patches/adds).

## Scripts (/workspace/scripts) and outputs (/workspace/out/q6r2)
- s160_q6r2_tb_scale.py -> s160.json: Xb4b cubes -> T_b maps -> flux-conserving elliptical-kernel smoothing to 55/120/250 pc circular beams and the Xb4d beam; 3-ch boxcar peak T_b; compared with native Xb4d.
- s161_q6r2_core_pedestal.py -> s161.json: beam-convolved core (elliptical) + fixed-width circular pedestal fits (1.2/1.8/2.5"), radial profiles.
- s162_q6r2_profile_manual_cont.py -> s162.json: radial half-power sizes, curve of growth, regular vs manual 26.6 MHz cube, native continuum peaks, T_b, dust Sigma.
- s163_q6r2_fit_virial.py -> s163.json, tb_vs_scale.png: T_b(theta) slopes and core+pedestal fits; virial variants.
- s164_update_report_q6r2.py: report patch.

## Results
1. T_b(scale) (K RJ, native 0.24" / 55 / 120 / 250 pc): HCO+ N 2.75/2.32/1.14/0.69; S 2.81/2.49/1.20/0.49; Bridge 1.87/1.75/1.13/0.56; E1 2.02/1.92/1.44/0.52; E2 2.04/1.67/0.64/0.27. HCN N 2.65/2.10/1.04/0.50; S 1.18/1.04/0.54/0.19; Bridge 1.13/1.01/0.60/0.28; E1 0.75/0.68/0.40/0.13; E2 0.74/0.63/0.30/0.12.
   Fall of 3.3-7.6x from 49 to 250 pc, slope T_b ~ theta^-0.7..-1.2 (point source -2, smooth medium 0). Xb4b smoothed to the Xb4d beam reproduces native Xb4d to <=15% (9/10 cases). The dense gas is therefore not a "smooth medium" as claimed for CO 2-1 by Brunetti+2021 (no trend 55-120 pc); here already 1.5-2.6x drop between 55 and 120 pc.
   Core+pedestal model (w=1.8" fixed): S HCO+ d=0.35", Tc 3.5 K, core filling in the 0.24" beam 0.68; HCN cores d 0.34-0.55", Tc 0.9-1.4 K; N is degenerate (a plateau profile, not core+Gaussian). The earlier "dilution only 1.2-1.35" holds for the ~100 pc core only; the report text was corrected. Filling at <50 pc is still not measurable (no flattening seen).
2. Size robustness: S HCO+ core 0.54" (pedestal 1.8-2.5") vs half-max 0.50" vs radial 0.48" -> adopted 0.50+-0.05" (97 pc); core share of the flux ~20-30%. N is unconstrained by the core+pedestal fit; half-max 0.52", radial 0.40", manual-cube half-max 0.36" -> adopted 0.45 (+0.07/-0.15)" and flagged as imaging dependent (difference >0.05"). Bridge/E1/E2 are 0.4-0.8" (elongation, blending). Peak K km/s agrees between the regular and manual cubes to 3-6%.
3. Continuum with native pixels: S 1.85 mJy/beam (rms 0.023), T_b 0.59 K, dust-only Sigma 5.0e3 (Td 50 K) / 2.3e3 (100 K) Msun/pc2 (upper limit); N 0.56, SE 0.56, E1 0.39, E2 0.28, Bridge 0.19 mJy/beam. The earlier smoothed 1.64 was 10% low.
4. Virial: M_vir=1040 sigma^2 R with six conventions (R=1.91 sigma_r vs half-max radius, sigma 1G vs mom2 vs gradient-corrected, radial-profile D, 35 Mpc): nuclei 0.9-3e8, clumps E1/E2 3-8e7 Msun; the CO median (5.4e7) is matched by E1/E2 to a factor 2 only; upper limits.

## Literature
Brunetti+2021 (2021MNRAS.500.4730B, abstract: smooth medium, no trend with scale 55-120 pc, T_b up to 37 K); Brunetti+2022 (2022MNRAS.515.2928B, catalogue medians); Sakamoto+2014; Rosolowsky & Leroy 2006 (virial convention).

## Uncertain
- Xb4b lacks scales >2-3": the 120/250 pc points are lower limits by ~<=15% (the true slope is ~0.1 shallower).
- Core+pedestal decomposition for N; the fitted pedestal width was fixed (1.2-2.5"), not derived from Xb4d (a free width was unconstrained).
- N size imaging-dependent (0.36-0.52"); weak E1/E2/HCN entries (S/N 6-12).
- Per-beam chi2 of the core+pedestal fits is 60-700: the structures are not Gaussians, the quoted fit errors are statistical only.
- Brunetti's CO values were taken from the abstract; no CO map was compared pixel-to-pixel with ours.
