# Q8 work summary: cross-band / cross-array flux consistency (worker, round 1)

## What I did
Convolved the Xb4b (0.25") and Xb4d (0.94x0.74") HCO+ 3-2 and HCN 3-2 pbcor cubes to the Xb4f header beam (5.73x3.54" / 5.72x3.54", PA 82) with an analytic Fourier Gaussian (flux rescaled by beam-area ratio), removed a per-pixel line-free baseline, and integrated (|v-VSYS|<350 core, -600..+550 all) at region peaks and in r=4/8/12" apertures. Also compared native-beam Xb4b vs Xb4d (existing s15), tapered/manual products vs standard cubes, and CH3OH 254 GHz Xb4d-smoothed vs Xb4f.

## Scripts and outputs (all re-runnable, `python3 -I scripts/<name>`)
- scripts/s26_flux_recovery.py -> out/s26_recovery.json/.txt, out/s26_spectra_{HCO+,HCN}.npz (aperture & peak spectra)
- scripts/s27_ch3oh_recovery.py -> out/s27_ch3oh_recovery.json
- scripts/s28_manual_tapered.py -> out/s28_manual_tapered.json/.txt
- scripts/s15_crossous.py (existing) -> out/s15_crossous.txt
- scripts/s29_update_report_q8.py -> /output/report.json (now 95 results, schema-valid; 4 new Q8 entries)

## Results
1. **Native-beam consistency**: at matched (0.94") beam Xb4b/Xb4d agree to 0-8% (HCO+ N 1.01, S 0.95, Bridge 0.92; continuum 0.92-1.02), consistent with 5% calibration. The differences below are therefore spatial filtering, not calibration.
2. **Recovery at the Xb4f beam** (core moment-0 at the peaks, vs Xb4f): Xb4d 0.87 (N), 0.81 (S), 0.87 (Bridge) for HCO+; HCN 0.87/0.77/0.91. Xb4b 0.78/0.67/0.76 (HCO+), 0.78/0.61/0.80 (HCN). Faint clumps (E1,E2,W,Sblue): Xb4d 0.3-0.85, Xb4b 0.15-0.7 (noisy, max-in-1" bias).
3. **Apertures** (HCO+ core Jy km/s): r<8": Xb4f 94, Xb4d 61 (0.65), Xb4b 44 (0.47); r<12": 94, 47, ~0 (interferometer bowl; zero-spacing missing). So ~35% (Xb4d) to ~55% (Xb4b) of the Xb4f flux within r<8" is resolved out. HCN: 0.58 / 0.44 at r<8". Xb4f itself has no short spacings; its own missing flux is unknown (aperture flux saturates at r~8").
4. **Line and continuum are filtered alike** (continuum peak Xb4d/Xb4f 0.84 N, 0.79 S; Xb4b 0.73, 0.68). **HCN/HCO+ 3-2 is robust**: N 0.600/0.598/0.599 (Xb4f/Xb4d/Xb4b), S 0.380/0.360/0.347, Bridge 0.513/0.535/0.544 -> cross-array spread <=2% N, <=9% S, <=6% Bridge. For faint clumps (E2, Sblue) arrays differ up to 2x: do not mix arrays there.
5. **Tapered/manual vs standard**: Xb4d tapered 2-12% above standard-convolved (HCO+ 1.02-1.10, HCN 0.99-1.12); Xb4b manual vs standard 0.97-1.00 at bright regions, 1.23 at faint Sblue. Usable with a 10-20% systematic for faint wings.
6. **CH3OH 254 GHz**: after smoothing Xb4d to the Xb4f beam, Bridge 0.46 vs 0.60 (0.77+-0.3), S 0.67, W 0.73, Sblue 0.82, similar to HCO+ (0.83-0.94). The earlier "Xb4d gets 1/3 of Xb4f" was beam dilution of extended CH3OH, not resolved-out flux; at most 20-30% (~1.5 sigma) is resolved out. Closes the Q0 residual.

## Implications
Matched-beam ratios from one array are reliable at N/S/Bridge; ratios across arrays or Xb4f-vs-high-res comparisons must carry a 10-35% recovery systematic for extended gas. For Q1, use Xb4d (or Xb4b smoothed) for resolved N/S/Bridge ratios; extended-disk fluxes are only in Xb4f/Xb47/Xb55.

## Uncertain
- Faint-region peak ratios are positively biased (max within 1") and noise-limited.
- Per-channel beam variation of Xb4f (header beam used) <~2% not corrected.
- Cannot say how much total flux Xb4f misses (no total-power/short-spacing data); 12" aperture totals of Xb4b/Xb4d are bowl-dominated and not physical fluxes.
- Band 5/7 vs Band 6 flux scale (needed for Q3) not examined here.
