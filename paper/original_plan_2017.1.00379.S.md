# Dense-Gas Conditions in the NGC 3256 Merger and Its AGN-Driven Outflow

**ALMA project:** 2017.1.00379.S — Physical properties of dense gas in an AGN-driven outflow  
**PI:** Harada, Nanase | **Targets:** ngc_3256 | **Band(s):** 5, 6, 7 | **Public since:** 2019-04-16  
**Publication status:** unpublished (see `input_files/literature.md` for the evidence)

## Abstract

Galactic outflows are thought to quench star formation by expelling molecular gas, yet star formation has been detected within some AGN-driven outflows, plausibly triggered by gas compression; how density, temperature, and chemistry respond in dense outflowing gas remains poorly constrained. NGC 3256, a nearby luminous infrared merger with two nuclei and a known molecular outflow, is an ideal test case, and archival ALMA project 2017.1.00379.S provides Band 5/6/7 cubes spanning the HCN and HCO+ J=2-1, 3-2, and 4-3 ladders plus CS(7-6), H13CN(4-3), and SO(6_5-5_4). We will re-image these lines in CASA, smooth them to a common 7.01-arcsec beam, and extract fluxes in apertures pre-committed via a resolution-honest Step-0 diagnostic comparing the nucleus-outflow separation to the beam. Fixed-temperature (30/60/100 K) RADEX LVG modeling will constrain gas density and dense-gas fraction per aperture, with opacity and shock cross-checks from H13CN, CS, and SO. We expect to determine whether the active region differs from an outer-disk control and, regardless of spatial resolvability, to deliver the first published catalog of these diagnostic lines for NGC 3256, informing the star-forming potential of AGN-driven outflows.

## Proposed work and key results

## Stage 1: Data acquisition (Weeks 1–2, ~15 hrs)
Download all five member OUS products for project 2017.1.00379.S from the ALMA Science Archive: `uid://A001/X1273/Xb47` (Band 7), `uid://A001/X1273/Xb4b` and `uid://A001/X1273/Xb4d` and `uid://A001/X1273/Xb4f` (Band 6, three resolutions), and `uid://A001/X1273/Xb55` (Band 5). Retrieve both the pipeline-calibrated measurement sets and the QA2 image products. Organize by band/resolution in a consistent directory tree.

## Stage 2: Inspection and QA (Weeks 2–3, ~15 hrs)
Load QA2 cubes and continuum images in **CARTA** to check image fidelity, field coverage, and obvious artifacts. Use `listobs` and `plotms` in **CASA** to confirm spectral windows against the expected line list (verify each transition's rest frequency in **Splatalogue**, since the provided list is unverified). Note synthesized beams per member OUS — this fixes the 7.01″ common-beam ceiling.

## Stage 3: Step-0 aperture diagnostic (Weeks 3–4, ~15 hrs)
**Figure 1:** overlay the dual-nucleus positions (Sakamoto et al. 2006), the outflow position (Michiyama et al. 2018), and the CO-outflow morphology (Brunetti et al. 2024) on the Band 5 (7.01″-beam) continuum or line map, with the beam FWHM drawn to scale, using **astropy** (WCS, `SkyCoord`) and **matplotlib**. Measure nucleus–outflow separation in beams and pre-commit to the aperture scheme this dictates (either two apertures — "active region" and "disk control" — or beam-filling-factor-corrected native-resolution apertures). Document the decision before any flux is extracted.

## Stage 4: Re-imaging (Weeks 5–9, ~50 hrs)
In **CASA**, run `tclean` on all six HCN/HCO+ transitions plus CS(7–6), H¹³CN(4–3), and SO(6₅–5₄) from calibrated measurement sets, then `imsmooth` each cube to the common 7.01″ beam. Image the Band 5/6/7 continuum similarly. This is the largest time sink; budget extra weeks for iterative cleaning and self-checks against QA2 products.

## Stage 5: Spectral extraction and moment maps (Weeks 10–12, ~30 hrs)
For each pre-committed aperture, extract integrated spectra using **spectral-cube** and **astropy**; compute moment-0/1 maps in **CASA** (`immoments`) or spectral-cube for visual context. **Figure 2:** integrated spectra of the HCN/HCO+ ladder per aperture, stacked for comparison.

## Stage 6: Flux and line-ratio measurement (Weeks 12–13, ~15 hrs)
Measure integrated fluxes (Gaussian or velocity-window integration) with uncertainty propagation from the archive's line sensitivities; report 3σ upper limits for non-detections. Compute HCN/HCO+ ratios at matched J per aperture, extending **Figure 2** with a ratio bar chart.

## Stage 7: RADEX LVG modeling (Weeks 1–4 of Sem. 2, ~35 hrs)
Run a fixed-temperature RADEX grid (T = 30/60/100 K, density 10²–10⁷ cm⁻³) in **Python** (`pyradex` or direct RADEX calls) against the HCN and HCO+ ladders per aperture. **Figure 3:** density-vs-temperature LVG solution plot with χ² contours per aperture.

## Stage 8: Opacity cross-check (Weeks 5–6, ~15 hrs)
Compare HCN(4–3)/H¹³CN(4–3) and HCN/CS(7–6) ratios per aperture to test the optically-thin assumption underlying Stage 7; flag apertures where HCN(4–3) opacity biases the LVG solution.

## Stage 9: SO shock-tracer analysis (Weeks 6–7, ~15 hrs)
Measure SO(6₅–5₄) flux or upper limit per aperture; compare SO/HCN as a shock-chemistry indicator distinguishing compression-driven from purely gravitational gas.

## Stage 10: Summary table and literature benchmarking (Weeks 8–9, ~15 hrs)
**Figure 4:** combined table/figure of fluxes, ratios, LVG densities/temperatures, opacity diagnostics, and SO/HCN per aperture, plotted against literature dense-gas values (e.g., other LIRG outflow studies) for context.

## Stage 11: Interpretation (Weeks 10–11, ~15 hrs)
Synthesize whether density/temperature/chemistry differ between active-region and disk-control apertures, or whether the beam-limited result requires a single-region interpretation, tying back explicitly to the Stage 3 decision.

## Stage 12: Thesis writing (Weeks 12–15, ~35 hrs)
Draft chapters around Figures 1–4 and the summary table, emphasizing the catalog of fluxes/ratios/densities/upper limits as the guaranteed deliverable regardless of spatial resolvability.

## Background reading

1. Michiyama et al. (2018), *ALMA Observations of HCN and HCO<SUP>+</SUP> Outflows in the Merging Galaxy NGC 3256* — [2018ApJ...868...95M](https://ui.adsabs.harvard.edu/abs/2018ApJ...868...95M)
   — The direct predecessor study of HCN/HCO+ outflow tracers in NGC 3256 whose outflow position and morphology this project's Step-0 aperture diagnostic is built on — the closest work not to duplicate.
1. Pereira-Santaella et al. (2024), *The CO-to-H<SUB>2</SUB> conversion factor of molecular outflows. Rovibrational CO emission in NGC 3256-S resolved by JWST/NIRSpec* — [2024A&A...681A.117P](https://ui.adsabs.harvard.edu/abs/2024A&A...681A.117P)
   — JWST-resolved molecular outflow analysis of NGC 3256-S itself, the most recent and directly competing characterization of dense/warm gas in this exact outflow.
1. Ledger et al. (2021), *Observed CN and HCN intensity ratios exhibit subtle variations in extreme galaxy environments* — [2021MNRAS.504.5863L](https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.5863L)
   — Reports HCN/CN and HCN(3-2)/(1-0) intensity ratio trends specifically in NGC 3256, the closest existing line-ratio dataset this project's HCN/HCO+ ladder measurements must be checked against.
1. Harada et al. (2018), *ALMA Astrochemical Observations of the Infrared-luminous Merger NGC 3256* — [2018ApJ...855...49H](https://ui.adsabs.harvard.edu/abs/2018ApJ...855...49H)
   — Multi-transition, multi-species ALMA astrochemical survey of NGC 3256 that establishes the baseline molecular composition and line list the project's re-imaged transitions build on.
1. Brunetti et al. (2021), *Highly turbulent gas on GMC scales in NGC 3256, the nearest luminous infrared galaxy* — [2021MNRAS.500.4730B](https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.4730B)
   — Establishes the GMC-scale physical conditions (turbulence, surface density) of NGC 3256's molecular ISM that provides the disk-control-region context for interpreting the outflow apertures.
1. Michiyama et al. (2020), *Star Formation Traced by Optical and Millimeter Hydrogen Recombination Lines and Free-Free Emissions in the Dusty Merging Galaxy NGC 3256—MUSE/VLT and ALMA Synergy* — [2020ApJ...895...85M](https://ui.adsabs.harvard.edu/abs/2020ApJ...895...85M)
   — Characterizes nuclear vs. disk-wide star formation in NGC 3256 via recombination lines, giving the star-formation context needed to interpret whether outflow gas is star-forming.
1. Nishimura et al. (2024), *CON-quest. II. Spatially and spectrally resolved HCN/HCO<SUP>+</SUP> line ratios in local luminous and ultraluminous infrared galaxies* — [2024A&A...686A..48N](https://ui.adsabs.harvard.edu/abs/2024A&A...686A..48N)
   — Demonstrates the HCN/HCO+ line-ratio method as a dense-gas and excitation diagnostic across a sample of local U/LIRGs, the core measurement technique this project applies per aperture.
1. Sato et al. (2022), *APEX and NOEMA observations of H<SUB>2</SUB>S in nearby luminous galaxies and the ULIRG Mrk 231. A possible relation between dense gas properties and molecular outflows* — [2022A&A...660A..82S](https://ui.adsabs.harvard.edu/abs/2022A&A...660A..82S)
   — Uses a dense-gas/shock tracer (H2S, analogous to this project's SO) to link dense gas properties to molecular outflows in LIRGs, motivating the SO/HCN shock-chemistry diagnostic in Stage 9.
1. Otter et al. (2026), *Clumpy, Dense Gas in the Outflow of NGC 1266* — [2026ApJ...997..361O](https://ui.adsabs.harvard.edu/abs/2026ApJ...997..361O)
   — Recent study of clumpy dense gas in an AGN-driven outflow (NGC 1266) directly bearing on the science question of whether outflow compression can trigger star formation.
1. Harada (2018), *High-Temperature Chemistry in External Galaxies* — [2018IAUS..332...25H](https://ui.adsabs.harvard.edu/abs/2018IAUS..332...25H)
   — Review of high-temperature/shock/UV/X-ray/cosmic-ray chemical heating mechanisms in external galaxies, giving the chemistry background needed to interpret HCN, HCO+, CS, and SO abundance responses.

---
*Draft thesis project proposal generated by almathesis on 2026-07-15; project idea and plan derived from the ALMA Science Archive and NASA ADS. All references were selected from ADS query results.*
