# q05: literature checks for Sec. 3.5

Sources were read from the papers themselves (full text via `lit.py` and arXiv PDFs). Our numbers come from `scripts/q53*.py` and `out/q53/`. I could not see Table C4 or the draft. Item (c) rests on what I could reconstruct from our outputs.

## (a) Nishimura et al. 2024 (A&A 686, A48; 2024A&A...686A..48N; arXiv 2402.15436)
- **Quantity.** The ratio is L'_HCN/L'_HCO+ for J=3-2. L' is in K km/s pc², so it is the ratio of velocity-integrated Rayleigh-Jeans brightness (flux × λ²), not a peak ratio. The Rayleigh-Jeans form is the same T_b scale we use. Only the tiny frequency factor (267.558/265.886)² = 1.013 separates L' ratios from integrated-T_b ratios.
- **Velocity window.** Flux densities are integrated over a window about the systemic velocity, chosen per galaxy: ±150, ±200, ±250, ±300 or ±500 km/s (Table 3). Ours is a fixed ±350 km/s, the same for both lines and both nuclei (`q53e/g`).
- **Apertures.**
  - "Resolved": an aperture the size of the synthesized beam (about 0.3", 50–100 pc for LIRGs, 200–400 pc for ULIRGs).
  - "Global": an elliptical aperture 10×10 times the beam, i.e. 0.5–1 kpc for LIRGs and 2–4 kpc for ULIRGs.
  - Beams are natural-weighting, about 0.3". There is no common-beam taper across galaxies. The ratio ranges are 0.4–2.3 (global, mean 1.1) and 0.4–2.7 (resolved, mean 1.2), from Table 4.
- **NGC 3256.** It is not in their sample. They cite it (Harada+2018) as a case where an elevated ratio is unrelated to an AGN.
- **Comparison with ours (HCN/HCO+ 3-2, T_b, ±350 km/s, `out/q53/lit_table.csv`, `q53g`).** Their "resolved" scale matches our 0.25–0.3" core point. Their "global" scale of 0.5–1 kpc is closest to our 1.2 kpc aperture (about 6").

| scale | N | S |
|---|---|---|
| core, r<0.25" (56 pc) | 0.98 ± 0.02 | 0.42 ± 0.02 |
| 200 pc | 0.81 ± 0.03 | 0.40 ± 0.02 |
| 1.2 kpc | 0.52 ± 0.03 | 0.42 ± 0.03 |

Caveats:
- Our core aperture is a radius of 0.25", i.e. a diameter of about 2 beams. Theirs is one beam diameter.
- Our values are fixed-position or aperture values in the pbcor cubes, convolved to a common circular beam.
- At N the ratio is beam dependent (1.0 → 0.5). The core is the right comparison for their "resolved" ratio, and the 1.2 kpc value for their "global" ratio.
- Both N and S fall inside the 0.4–2.7 range. The N core is at the sample mean (about 1.0–1.2), and S is at the low end.
- The Nishimura sample has no AGN/SB split of the ratio. They report no systematic trend with L_bol,AGN.

## (b) Izumi et al. 2016 (ApJ 818, 42; 2016ApJ...818...42I; arXiv 1512.03438), Table 1
- **What the ratio is.** R_HCN/HCO+ is the HCN(4-3)/HCO+(4-3) integrated-intensity ratio on the brightness-temperature scale. It is NOT J=3-2. They add a systematic error (15% if not given).
- **Starburst range 0.24–1.03 is confirmed, with a caveat.**
  - It is the range of the high-resolution (<500 pc) SB entries other than the two LMC clumps: NGC 1614 0.24 ± 0.06, M82 0.41, NGC 7469 ring C/D 0.45/0.48, B 0.75, NGC 1097 ring 0.82, NGC 253 1.03 ± 0.22.
  - M82 and NGC 253 are JCMT/APEX single-dish points (353 and 262 pc).
  - The LMC clumps N113 (0.21) and N159 (0.11) are also in that table as SB. Including them gives 0.11–1.03.
  - Our median of 0.48 (`q53i`) uses the 0.24–1.03 set.
- **AGN range (high-resolution ALMA).** 1.11–3.19: NGC 7469 1.11, NGC 1068 AGN position 1.53, CND-S 2.58, E-knot 2.84, CND-N 3.14, W-knot 3.19, NGC 1097 2.01. The median is 2.6, and the text states that AGNs show higher ratios than SB.
- **Low-resolution entries (>500 pc).** These mix the two classes. Examples are NGC 1365 0.59, NGC 4945 0.78 (AGN, SB-ring contaminated), NGC 3256 whole galaxy 0.39 ± 0.11, and IRAS 13242 SB 1.23. Our classification used only the high-resolution ranges.
- **Use for NGC 3256.** Compare only our 4-3 values (not 3-2) to these ranges: N 0.80 (0.69–0.93) at 188 pc and about 0.95 (0.81–1.12) at the 56 pc core. S is 0.23 (0.18–0.29) at 188 pc and 0.29 (0.25–0.35) at the core (`out/q53/mc.json`, report entry Q53).

## (c) 4-3 ratio of S: 0.29 versus 0.317
**I believe the two numbers are different quantities, and 0.29 is not a 4.1×2.7" measurement.**
- **0.317 ± 0.011** is the S HCN 4-3/HCO+ 4-3 T_b ratio from the single Xb47 tuning, convolved to a circular 4.5" beam. The window is ±350 km/s with a linear baseline (`out/q53/fixedwin43.csv`, row HCN47/HCO47, 4.5, S). N is 0.515 ± 0.008.
- **0.29 in our outputs** is the extrapolated S core value (0.3", 56 pc): 0.29 (0.25–0.35). It equals the 3-2 core ratio times R43/R32 with a flat-to-linear extrapolation (`scripts/q53h_mc.py`, entry Q53, `out/q53/mc.json`). It is not a 4.1×2.7" value. Please check whether Table C4 mislabelled this value.
- **Native-beam Xb47 S values I can reproduce**
  - Peak-pixel T_b ratio at 4.13×2.70": 0.039/0.130 = 0.30 (`out/q4/peak_tb.csv`; N 0.54).
  - Fixed-position ratio at the common 6.5" and 10" beams: 0.302 and 0.305 (`out/q1/q1_ratios.csv`).
  - The 0.317 is at 4.5".
- **Reason for 0.30 vs 0.317.** The estimators differ (peak pixel vs a ±350 km/s integral with a baseline). The beams also differ (native elliptical vs circular 4.5"). At S the HCN profile is narrower than HCO+ (about 120 vs 192 km/s FWHM) and offset by about 25 km/s (`out/q53/spec43_overlay.png`). The ratio of peaks therefore differs from the integrated ratio at the 5% level. Treat the difference as a method systematic of about 0.02, not a physical beam effect.
- **Separate point worth stating.** The Xb47 ratio at 4.5" (S 0.32) is higher than the ratio from the high-resolution 4-3 cubes at the same beam (S 0.22). The high-resolution cubes recover only 61% of the Xb47 S HCN flux (HCO+ 87%). The 4-3 S ratio at ≤1" is therefore uncertain from 0.22 up to about 0.32 (see `paper_answers/r08/answer.md` §2).

## (d) Bibliographic data and numbers to compare
All bibcodes were checked against the paper library headers.

**Falstad et al. 2021.** A&A 649, A105 (article number), bibcode 2021A&A...649A.105F, arXiv 2102.13563. Title: "CON-quest: Searching for the most obscured galaxy nuclei".
- **Not an NGC 3256 paper.** NGC 3256 does not appear in its text in the library. It is the CON-quest survey of HCN-vib J=3-2, with a CON criterion of Σ_HCN-vib > 1 L⊙ pc⁻² for r > 5 pc. If the draft takes N/S numbers from Falstad, the source is probably wrong, so please check.
- **Our HCN-vib limits.** N L_vib < 29 L⊙ (48 pc), S < 25 L⊙; Σ_vib < 0.02 L⊙ pc⁻², at least 60× below the CON criterion (r08 §3). These are our limits, not Falstad numbers.

**Ledger et al. 2021.** MNRAS 504, 5863, bibcode 2021MNRAS.504.5863L, arXiv 2104.08224. NGC 3256 data are at a common 2.2" beam (469 pc), 1.1" pixels, with a uv cut at 15 kλ (Table 3). Nuclei are 3-pixel-type apertures centred on the 93 GHz peak. Intensities are K km/s.

| ratio | Ledger NGC 3256 | ours |
|---|---|---|
| HCN 3-2/1-0 | N 0.39 ± 0.04, non-nuclear 0.27 ± 0.03; S not given in the text I read (the abstract gives ~0.3–0.4); slope 1.53 ± 0.07 | S 0.335 ± 0.028; 0.31–0.36 at both nuclei (1-0 from Harada+2018, 3-2 from Xb4b/Xb4d); N equal to S within errors |
| HCN 4-3/3-2 | N 0.48 ± 0.07, S 0.46 ± 0.07, non-nuclear 0.45 ± 0.06 (slope 0.99 ± 0.04) | not at matched beam; our R43/R32 depends on the missing flux (r10) |
| CN 1-0/HCN 1-0 | N 1.7 ± 0.1, non-nuclear 1.14 ± 0.09, S about 1.1 | not measured by us |
| CN 2-1/1-0 | 0.95 ± 0.03 (slope) | not measured |

- **Ledger's inference.** They attribute the within-galaxy variation to excitation at the two nuclei. Our HCN 3-2/1-0 is equal at N and S, so it does not need an N/S excitation difference.
- **Beam.** At 2.2", Ledger's N and S are separate (5" apart) but each is beam-diluted. Ours is quoted at 1.7–2.5".

**Pereira-Santaella et al. 2024.** A&A 681, A117 (article number), bibcode 2024A&A...681A.117P, arXiv 2309.06486. JWST/NIRSpec CO v=1–0 at 4.67 µm around NGC 3256-S.
- **Obscured nucleus.** The bright 4.6–4.7 µm continuum is dust-embedded. CO v=1–0 is in absorption in both P and R branches toward the nucleus. The disk shows P-R asymmetry and the outflow shows emission.
- **Excitation by the AGN.** The 4.7 µm continuum vibrationally excites CO in the outflow out to r ~ 250 pc. Gas temperature from P(2)/R(0): above 40 K in the central 100 pc and below 15 K at 250 pc, in 100 pc (0.5") apertures. The paper gives no numbers for HCN.
- **Comparison with our results.** The statement "S is dust-embedded and has an obscured 4.7 µm continuum" is theirs. Our result is the HCN-vib limit: Tvib below about 300–430 K, no CON-level pumping at 48 pc. It does not exclude a buried source smaller than 47 pc that is diluted in the beam.
- **Cautions.** The 40 K and 15 K temperatures are for the CO outflow, not for dense gas at the nucleus. They do not constrain HCN excitation at N or S.

## Assumptions and checks
- Distances and scale: 1" = 188 pc (our adopted scale), whereas Pereira-Santaella adopt d = 40 Mpc (about 194 pc/").
- The Nishimura beam, aperture and window definitions are from their Sects. 4.1 and 4.3 and Table 3.
- Izumi values are transcribed from their Table 1 (`scripts/q53_context.py`).
