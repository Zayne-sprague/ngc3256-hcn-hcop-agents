# The AGN nucleus is the low one: HCN/HCO$^+$ J=3–2 is about twice as high in the starburst nucleus of NGC 3256 as at its AGN on 50–200 pc scales, because HCN, HNC and CN are weak relative to HCO$^+$ at the AGN

**Agent teams A and B (forum), written by the paper team**

## Abstract

An elevated HCN/HCO$^+$ line ratio in the submillimetre lines is widely used to flag accreting black holes in dusty galaxy nuclei, but other studies attribute high ratios to chemistry, shocks or infrared pumping, and the diagnostic has rarely been tested where two nuclei of one galaxy can be separated. We use archival ALMA Band 6 data of the merger NGC 3256, which hosts a starburst-dominated northern nucleus (N) and an obscured, low-luminosity AGN in the southern nucleus (S), to compare the two at matched resolution of 50–250 pc. The brightness-temperature ratio HCN/HCO$^+$ J=3–2 is $0.99\pm0.14$ at N and $0.44\pm0.06$ at S, a contrast of $2.1\pm0.2$ at 50–120 pc, with one-sided, untested uncertainty from missing short spacings (up to $-0.3$). The contrast sits in HCN, not HCO$^+$: the HCO$^+$ peak brightness is the same at both nuclei, whereas HCN is twice as bright at N. HNC and CN show the same deficit at S relative to HCO$^+$, N is HCN-rich inside about 100 pc, and no vibrationally excited HCN emission is detected at either nucleus, excluding CON-level pumping. Non-linear radiative-transfer fits mildly (about $2\sigma$) disfavour excitation alone, but cannot separate abundance from opacity. The ratio at S sits at the bottom of the distribution of resolved luminous-infrared-galaxy nuclei, while N is typical. A single-line HCN/HCO$^+$ ratio therefore does not flag, and here inverts, the AGN/starburst order.

## 1 Introduction

The HCN/HCO$^+$ ratio of rotational transitions is a long-used diagnostic of the physical and chemical state of dense gas in galaxies. Observations of the J=3–2 and 4–3 transitions led several groups to report a "submillimetre HCN enhancement": higher HCN/HCO$^+$ in nuclei hosting active galactic nuclei (AGN) than in starbursts, which was interpreted as the signature of X-ray or mechanical heating, or of chemistry driven by the AGN radiation field (Izumi et al. 2016, ApJ 818, 42; Imanishi et al. 2016, AJ 152, 218). The claim, however, is contested. Resolved ALMA observations of luminous and ultraluminous infrared galaxies find a wide spread of ratios with no trend with AGN luminosity or star-formation rate, and associate the highest ratios with shock-heated outflows (Nishimura et al. 2024, A&A 686, A48). High ratios can also arise from abundance differences (Butterworth et al. 2022, A&A 667, A131), or from infrared pumping of HCN through its bending mode in compact, heavily obscured nuclei (CONs; Aalto et al. 2015, A&A 584, A42), in which case vibrationally excited HCN lines should be visible and the J=3–2 ratio would reflect radiative excitation rather than gas conditions.

NGC 3256 is a good laboratory for this question. It is a nearby ($D\approx40$ Mpc, $1''\approx194$ pc) merger with two nuclei separated by about 5$''$ (1 kpc): a northern nucleus that is starburst-dominated (N) and a southern nucleus that hosts a low-luminosity, obscured AGN (S). Published ratios of the two nuclei come from beams of about 1.7$''$ (330 pc), in which the nuclei differ by only a factor of about 1.5 and no enhancement at the AGN was seen (Harada et al. 2018, ApJ 855, 49); J=1–0 values at about 2$''$ were given by Michiyama et al. (2018, ApJ 868, 95). What was missing is the comparison at 50 pc, where each nucleus is resolved from its surroundings; the information on whether the difference between the nuclei lies in HCN or in HCO$^+$; and tests of whether it is produced by excitation (rotational or vibrational), opacity or abundance. The ALMA archive lists no publication of the Band 6 programme used here.

In this paper we ask: is an elevated HCN/HCO$^+$ ratio a signature of AGN heating in dusty nuclei, and what sets the ratio in the two nuclei of NGC 3256? Our answer is that the ratio is about twice as high at the starburst nucleus as at the AGN at 50–200 pc scales, that this is because HCN (together with HNC and CN) is weak relative to HCO$^+$ at the AGN, that no vibrational-pumping signature is present, and that the single-line ratio therefore does not mark the AGN; the cause of the difference is not yet established.

## 2 Data and methods

**Data.** The main data are ALMA Band 6 observations of NGC 3256 (project 2017.1.00379.S) covering HCN and HCO$^+$ J=3–2 (265.9 and 267.6 GHz rest) in three array configurations: a long-baseline 12-m configuration giving a native beam of $0.25''$ (about 50 pc), a more compact 12-m configuration, and the 7-m array. Complementary data are HCN, HCO$^+$, H$^{13}$CN, H$^{13}$CO$^+$, HNC and CN J=1–0 from projects 2015.1.00993.S and 2015.1.00412.S, HNC and CN J=3–2 from 2015.1.00412.S and 2018.1.00493.S, and HCN and HCO$^+$ J=4–3 from 2018.1.00493.S and 2018.1.00223.S. The 1.3 mm continuum comes from the Band 6 data. We adopt $v_\mathrm{sys}=2780$ km s$^{-1}$, $D=40$ Mpc, and a flux-scale uncertainty of 10% per line (14% per ratio); the flux scale multiplies N and S alike and cancels in all N/S comparisons.

**Measurements.** Primary-beam-corrected cubes are converted to Rayleigh–Jeans brightness temperature $T_b$, convolved with a flux-conserving kernel to a common circular beam, and baseline-subtracted using line-free channels. Velocity-integrated values are taken over $\pm350$ km s$^{-1}$ about the systemic velocity (HCN regridded onto the HCO$^+$ velocity grid). The line ratio $R=T_b(\mathrm{HCN})/T_b(\mathrm{HCO^+})$ is the ratio of these integrals in an aperture of radius equal to the beam FWHM centred on the HCO$^+$ peak of each nucleus; the N/S contrast is $C=R_N/R_S$. Uncertainties are derived from line-free channels. We vary the estimator (aperture, single pixel, peak channel), window ($\pm150$ versus $\pm350$ km s$^{-1}$) and baseline, and treat the spread as a systematic (Appendix C).

**Missing flux.** Interferometric images lose extended emission. This matters here because S is more extended and fainter than N and loses more flux, which inflates $C$. Comparing the two 12-m configurations shows a change of $-0.06$ in $C$ at 1.0$''$. For beams of 1.5–2.2$''$ we have flux-complete 12-m+7-m images, and these lower $C$ by 0.14–0.30. No flux-complete test exists at beams $\le0.6''$, which is the main untested systematic of the core contrast.

**Modelling.** Physical conditions are explored with one-zone non-LTE radiative-transfer (RADEX) models fitted jointly to the 1–0 and 3–2 lines of HCN, HCO$^+$ and their $^{13}$C isotopologues in a matched 2.2$''$ beam, and compared with Bayes factors. Details are in Appendix A.

## 3 Results

### 3.1 The N/S contrast is real, compact, and lives in HCN

**HCN/HCO$^+$ J=3–2 is $0.99$ at N and $0.44$ at S at 50 pc, a contrast of $2.1\pm0.2$ that is carried by HCN alone.** The ratio map at 0.26$''$ shows a compact high-ratio region centred on N and uniformly low values at S and along the bridge between the nuclei (Fig. 1). The two nuclei have the same HCO$^+$ peak brightness ($2.5$ versus $2.7$ K, N/S $=0.95$), but the HCN peak is $2.0\times$ lower at S ($1.2$ versus $2.5$ K); the difference is therefore in HCN, not in HCO$^+$.

![Figure 1. Peak brightness temperature $T_b$ (RJ, K) of HCO$^+$ J=3–2 (left) and HCN J=3–2 (middle), and the HCN/HCO$^+$ $T_b$ ratio (right), NGC 3256, at a common 0.26$''$ (50 pc) circular beam. Pixels below 3σ are black (left, middle) or grey (right). Crosses mark the northern (N) and southern (S) nuclei and the bridge between them; offsets are in arcsec (east to the left) and pc (1$''$ = 194 pc). The circle is the beam. Ring medians within 0.13$''$ of the markers give ratios of 0.99 (N), 0.48 (S) and 0.51 (bridge); the peak-channel values quoted in the text are maxima and are higher.](figs/fig1_overview.png)

The contrast depends on beam size, and we quote it as such (Fig. 2). In equal circular beams it is $2.26\pm0.06$ at 0.26$''$ (50 pc), $2.32$ at 55 pc, $2.12\pm0.04$ at 120 pc, $1.99\pm0.04$ at 1.0$''$ (194 pc), $1.90$ at 250 pc and $1.64$ at 2.2$''$ (430 pc) with 12-m data only; with flux-complete 12-m+7-m data at 2.2$''$ it is $1.44$–$1.50$. At the largest beams (6.5$''$, 1.3 kpc) N, the bridge and S blend and the ratio is about 0.4 everywhere. We adopt $C=2.1\pm0.2$ at 50–120 pc, where the statistical error is about 0.05, the spread between estimators is $\pm0.1$, and the error is widened to include the untested missing-flux term. That term is one-sided: it can only lower the contrast, by up to about 0.3 if the bias at 0.6$''$ resembled that measured at 1.5–2.2$''$. We therefore never quote the contrast as scale-independent.

![Figure 2. (a) HCN/HCO$^+$ J=3–2 $T_b$ ratio at N (orange), S (blue) and the bridge (green) as a function of beam FWHM, in an aperture of radius equal to the beam. Filled circles: long-baseline 12-m data; open squares: compact 12-m configuration at 1.29$''$; triangles at 6.5$''$: compact-array blends. Grey bands: ±14% flux scale, which cancels in N/S. (b) Peak $T_b$ of HCO$^+$ (solid) and HCN (dashed, open) at N and S. At 0.26$''$ the HCO$^+$ peak is equal at N and S while HCN is 2.0 times lower at S. Colours are colour-blind safe.](figs/fig2_ratio_vs_beam.png)

The beam dependence arises in the N core. Between 55 and 250 pc the N ratio falls from $0.99$ to $0.69$ (power-law index $-0.23\pm0.03$), whereas the S ratio is nearly flat ($0.43$ to $0.37$, index $-0.10\pm0.03$). In annuli around each nucleus, measured at 0.26$''$ resolution, the N ratio is $1.00\pm0.02$ inside 0.25$''$ (radius $<49$ pc), $0.88\pm0.02$ at 0.25–0.5$''$, $0.65\pm0.01$ at 0.5–1$''$, $0.43\pm0.02$ at 1–2$''$ and $0.38\pm0.05$ at 2–3$''$, while at S it is $0.44\pm0.01$, $0.39\pm0.02$, $0.33\pm0.02$ and $0.34\pm0.03$ over the same annuli (Fig. 3). The N nucleus therefore contains an HCN-rich core of about 100 pc diameter embedded in an envelope with ratio 0.4–0.6, and the contrast exists only inside about 1$''$ (190 pc). The outermost annulus at S (0.45) is a low-surface-brightness region dominated by extended gas and missing-flux effects and is not part of the compact profile. The two 12-m configurations agree within $0.03$ for the inner 1$''$ at N and S (ratio at N $0.84$ in both at 1.0$''$); at 0.5–1.5$''$ the compact configuration gives ratios $4$–$6\%$ higher, mainly because the long-baseline data lose about 15–20% of the flux at S. This lowers the contrast by only 3%, from $1.79$ to $1.74$.

![Figure 3. HCN/HCO$^+$ J=3–2 $T_b$ ratio (left) and annulus-mean integrated $T_b$ of HCN (solid) and HCO$^+$ (dotted) at N (middle) and S (right) as a function of radius from the HCO$^+$ peak. Blue: long-baseline data at 0.26$''$; orange: the same data convolved to 1.0$''$; green: the compact 12-m configuration at 1.0$''$ (flux cross-check). Error bars are statistical. Grey vertical lines mark the beam FWHM. At 1.0$''$ the annuli inside 0.5$''$ are not independent.](figs/fig3_ratio_profile.png)

Several further tests show that the contrast belongs to the nuclei and is not a by-product of the analysis. The ratio does not depend on the continuum brightness (at equal HCO$^+$ it is 0.63 at N and 0.43 at S); it is flat in velocity within each nucleus (0.95–1.09 at N, 0.36–0.43 at S); the contrast is the same in three observing blocks to within 9%; and HCN is not more extended than HCO$^+$ (both have Gaussian cores of 50–65 pc, half-maximum diameters 80–100 pc, marginally resolved). The rise of the ratio in the line wings is explained by beam mixing and the deficit at S exists at all velocities (Appendix A). **What this means:** at the resolution where the nuclei can be separated, the AGN nucleus has the lower HCN/HCO$^+$ ratio, by a factor of two, and this is a deficit of HCN, not an excess of HCO$^+$ or a continuum effect.

### 3.2 No vibrational-pumping signature; the same deficit appears in HNC and CN

**HCN, HNC and CN are all weak relative to HCO$^+$ at S, the effect is similar in the J=1–0 and J=3–2 transitions at flux-complete beams of 1.5$''$ and larger, and no vibrationally excited HCN is detected at either nucleus.**

*Species.* In a single pipeline with matched beams, the N/S contrast of the ratio to HCO$^+$ at 1–0 (2.5$''$ beam, flux-complete) is $1.56\pm0.03$ for HCN, $1.90\pm0.11$ for HNC and $2.10\pm0.04$ for CN (Fig. 4c, Table 1). At J=3–2 the corresponding values from 12-m data alone are $2.14\pm0.05$, $1.84\pm0.14$ and $2.56\pm0.09$ at 1.4$''$, and $1.93\pm0.05$, $2.26\pm0.21$ and $2.99\pm0.14$ at 2.5$''$; as these lack the 7-m data, they are upper values (the HCN/HCO$^+$ contrast at 2.2$''$ drops from 1.65 to 1.44–1.50 with 7-m data). Relative to HCO$^+$ the deficit at S is thus as large or larger for HNC and CN than for HCN, so HCN is not special. The N/S contrast of HNC/HCN is $1.22\pm0.07$ at 1–0 and consistent with unity at 3–2 ($0.86\pm0.06$ at 1.4$''$, $1.17\pm0.10$ at 2.5$''$), and CN/HCN is $1.2$–$1.55$. CN is therefore not enhanced at S relative to HCN, as an X-ray-dominated region would predict. Because these ratios compare lines at similar frequencies, the dust continuum cancels. Normalised instead to the 1.3 mm continuum, HCN, HNC and CN are 2.4–3.1 times weaker at S at J=3–2 (1.1–1.4$''$) and HCO$^+$ is $1.3$–$1.4$ times weaker; this factor scales with the unmeasured dust-temperature ratio $T_{d,N}/T_{d,S}$ (taken as 0.5–2), and we do not rely on it.

![Figure 4. (a) N/S contrast of the HCN/HCO$^+$ $T_b$ ratio at J=1–0 (blue) and J=3–2 (orange) for matched, flux-complete circular beams of 2.2$''$, 1.5$''$ and 2.5$''$; error bars are statistical with estimator spread. (b) Contrast at 1.0$''$ for 3-2 and the 4-3 transition, where the missing 4-3 flux is not constrained and the allowed range is 2.3–3.7. (c) N/S contrast of HNC, CN and HCN relative to HCO$^+$ at 1-0 (2.5$''$, flux-complete) and 3-2 (1.4$''$ and 2.5$''$; hatched, 12-m only, so upper values). Error bars are statistical.](figs/fig4_contrast_J_species.png)

**Table 1.** $T_b$ ratios (RJ, $\pm350$ km s$^{-1}$) and their N/S contrasts, measured at fixed nuclear positions in one pipeline. Errors are statistical. J=3–2 entries are 12-m only.

| ratio | 1–0, 2.5$''$: N | S | N/S | 3–2, 1.4$''$: N | S | N/S |
|---|---|---|---|---|---|---|
| HCN/HCO$^+$ | 0.912 | 0.585 | 1.56 ± 0.03 | 0.833 | 0.389 | 2.14 ± 0.05 |
| HNC/HCO$^+$ | 0.474 | 0.249 | 1.90 ± 0.11 | 0.477 | 0.259 | 1.84 ± 0.14 |
| CN/HCO$^+$ | 1.590 | 0.755 | 2.10 ± 0.04 | 0.893 | 0.348 | 2.56 ± 0.09 |
| HNC/HCN | 0.520 | 0.427 | 1.22 ± 0.07 | 0.572 | 0.666 | 0.86 ± 0.06 |
| CN/HCN | 1.743 | 1.292 | 1.35 ± 0.03 | 1.072 | 0.894 | 1.20 ± 0.04 |

*Rotational excitation.* At matched, flux-complete beams of 1.5–2.5$''$ the contrast is the same in 1–0 and 3–2: its ratio $D=C_{3-2}/C_{1-0}$ is $1.00\pm0.05$ with 12-m+7-m data and $1.07\pm0.05$ with a different estimator, adopted $D=1.04\pm0.07$ (Fig. 4a; 1–0 and 3–2 contrasts are $1.46\pm0.02$ and $1.48\pm0.04$ at 2.2$''$, $1.62\pm0.06$ and $1.60\pm0.07$ at 1.5$''$). Earlier estimates of $D\approx1.1$ and of a lower 3–2/1–0 ratio of HCN at S (0.71) came from 12-m images that lose flux at S, and are superseded (S/N 3–2/1–0 is $0.93\pm0.05$). This is a statement about the envelope-plus-core blend at 300–500 pc, where the contrast is about 1.5; it is not a test of the 2.1 core contrast, because no 1–0 observation at $\le1''$ is flux-complete. For J=4–3, the contrast at 1.0$''$ is 3.7 before and between 2.3 and 3.7 after allowance for missing flux (central value 3.3), against $2.1$ at 3–2; the ratio of the two is bounded between 1.07 and 1.73 (central 1.6), and a single-tuning 4.5$''$ observation gives a blended contrast of 1.6. The data cannot decide whether the contrast rises with J in the core (Appendix A).

*Vibrational pumping.* Vibrationally excited HCN, $v_2=1f$ J=3–2 at 267.199 GHz, is not detected at either nucleus (Fig. 5). The 3σ limits in the 48 pc beam are $\alpha=I(v_2{=}1)/I(v{=}0)<0.05$ at N and $<0.10$ at S, corresponding to vibrational luminosities $L_\mathrm{vib}<29\,L_\odot$ and $<25\,L_\odot$ and surface densities $\Sigma_\mathrm{vib}<0.02\,L_\odot$ pc$^{-2}$, at least 60 times below the CON criterion of $1\,L_\odot$ pc$^{-2}$ (Falstad et al. 2021). The limits on the ratio of the vibrational to ground-state line luminosity are $<0.025$–0.05 at N and $<0.085$–0.10 at S; the CON range in the literature is 0.06–0.18, so this limit excludes pumping as the explanation only for N, and only tentatively. Thus neither nucleus is CON-like, but our data do not exclude moderate pumping ($T_\mathrm{vib}=150$–300 K) or a buried source smaller than 47 pc diluted in the beam. Finally, there is no shock signature in the nuclei: the CS 7–6/HCN 4–3 ratio is uniform and SO/HCO$^+<0.014$.

![Figure 5. Spectra of the HCN $v_2=1f$ J=3–2 line (rest 267.199 GHz, green dotted) at the N (left) and S (right) nuclei in the 0.24$''$ (48 pc) beam. Grey: data (line-free rms 0.42 mJy beam$^{-1}$ per 3.9 MHz channel); blue: model of the HCO$^+$ red wing; orange: wing model plus the HCN $v=0$ profile scaled to the 3σ upper limit (0.033 and 0.090 for this wing model; 0.05 and 0.10 for the widest of four models). No feature is seen above the noise.](figs/fig5_hcn_vib.png)

**What this means:** the nitrogen-bearing species HCN, HNC and CN are all deficient at S relative to HCO$^+$, with no excess of CN or other signature of an X-ray-dominated region or shock; the contrast is shared between J=1–0 and 3–2 on 300–500 pc scales, and there is no CON-level infrared pumping at either nucleus. The effect is thus better described as a nitrogen-species-to-HCO$^+$ contrast than as an effect specific to HCN.

### 3.3 Cause: excitation alone is disfavoured; abundance is favoured but not established

**In a joint non-LTE fit at 2.2$''$, models in which the HCN/HCO$^+$ abundance ratio is equal at N and S fail to reproduce the observed contrast at about the $2\sigma$ level, whereas a higher HCN abundance at N reproduces it; opacity and the $^{12}$C/$^{13}$C ratio are not separated.** The fit uses the 1–0 and 3–2 lines of HCN, HCO$^+$, H$^{13}$CN and H$^{13}$CO$^+$ in a flux-complete 2.2$''$ beam, where the observed contrast is $1.45\pm0.01$ at 1–0 and $1.50\pm0.02$ at 3–2 (statistical). With equal abundances, the best excitation-only and column-density-only models reach contrasts of only 1.1–1.3, and the fit is worse than the full model by $\Delta\chi^2=5$–8 (1.7–2.1σ) at the minimum, rising to 10–18 at the observed value (Fig. 6). An abundance-only change fits with Bayes factors of 25–44 over excitation-only (2.6–2200 across error models, always positive), and implies column density ratios $N(\mathrm{HCN})/N(\mathrm{HCO^+})=3.2\ [2.5, 4.0]$ at N and $2.0\ [1.6, 2.5]$ at S, a factor of about 1.6 between them ($P>0.9$). A radiative-transfer calculation shows that moving the ratio at S from 0.4 to 1 requires a factor of 3–4 in abundance ratio.

![Figure 6. Left: $\Delta\chi^2$ relative to the full model as a function of the model 1–0 HCN/HCO$^+$ N/S contrast when the HCN/HCO$^+$ abundance ratio is forced to be equal at N and S (excitation-only, blue; column-density-only, orange; solid: $^{12}$C/$^{13}$C 40–100, dashed: 20–200). The observed contrast is shown in magenta. Squares/diamonds: abundance-only and full model best fits. Right: natural logarithm of the Bayes factor of abundance-only over excitation-only for different error models and priors.](figs/fig6_hypotheses.png)

The isotopologue data do not decide the issue. H$^{13}$CN/HCN is $0.022\pm0.008$ at N and $<0.03$ (3σ) at S, and H$^{13}$CO$^+$/HCO$^+$ is $0.028\pm0.006$ at N and $0.013\pm0.006$ at S (N/S $=2.2\pm1.2$). An abundance-only difference predicts equal values of the isotopologue ratios, while a different $^{12}$C/$^{13}$C ratio at the two nuclei (about 80 at N and $\gtrsim130$ at S) would produce a ratio of 1.6–2; our data are consistent with both. A single excitation state is rejected ($p<0.04$): the extended gas has $T_k\gtrsim80$ K at $n\sim10^4$ cm$^{-3}$, and a compact HCN-rich core of about 0.5$''$ is required, whose kinetic temperature is undetermined. Continuum and dust background limit the envelope background temperature to $\lesssim3.2$ K, disfavouring radiative pumping. The S core is cool ($T_k=17\pm1$ K), while the N temperature is unconstrained.

**What this means:** the contrast at the scale of the envelope (about 1.5 at 2.2$''$) is not easily explained by a difference in gas excitation between the nuclei, and a modest abundance excess of HCN, HNC and CN at N is the simplest explanation. The result is tentative, applies to the 2.2$''$ beam and not to the 2.1 core contrast, depends on the error model, and is a one-zone calculation.

### 3.4 Implication for the diagnostic: the AGN nucleus is at the bottom of the distribution

**At 50 pc the AGN nucleus S has one of the lowest resolved HCN/HCO$^+$ ratios of any luminous-infrared-galaxy nucleus, while the starburst nucleus N is typical.** Nishimura et al. (2024) list 21 resolved J=3–2 ratios of nuclei in luminous and ultraluminous infrared galaxies, ranging from 0.45 to 2.67 with a median of 1.13 (their ratios are not corrected for missing flux). Our N ratio, 0.99 at 55 pc, lies at the 38th percentile of this sample, and 0.84 and 0.69 at 120 and 250 pc lie at the 38th and 19th (Fig. 7a). S, with 0.43, 0.40 and 0.37 at the same scales, is at the minimum; the lowest literature value, 0.445 for ESO 148-IG002, is equal within the errors, so S is tied with it and not significantly below. The only AGN-dominated nucleus in the sample (IRAS F05189-2524) has 1.63. At 1–0, our ratios (0.84–0.94 at N, 0.57–0.58 at S) agree with earlier values of 0.96 and 0.645, and N is at the mean of the GOALS sample. At J=4–3 the beams are 4.1$''\times2.7''$ (800 pc), the nuclei blend, and N (0.52) and S (0.29) both fall inside the starburst range of Izumi et al. (2016; 0.24–1.03); at 188 pc the values are 0.81 and 0.22. S has the stronger AGN and free-free continuum yet the lower ratio, which disfavours AGN- or XDR-driven enhancement of HCN. At 6.5$''$ the two nuclei blend and give the same ratio (0.44), so that beam hides the N core.

![Figure 7. HCN/HCO$^+$ in context. (a) J=3–2 $T_b$ ratio versus spatial scale: circles are the 21 resolved nuclei of Nishimura et al. (2024) coloured by AGN bolometric fraction (grey <0.10, blue 0.10–0.25, pink diamond ≥0.25; open: not listed); stars are N (orange) and S (blue) at 55, 120 and 250 pc, and open stars are 6.5$''$ blends affected by flux loss. Grey shading marks ratios above 1.1. (b) J=4–3 ratios of N, S and the combined system (800 pc beam) on the AGN and starburst bands of Izumi et al. (2016); grey dots are individual high-resolution sources of that work.](figs/fig7_literature.png)

**What this means:** a single-transition HCN/HCO$^+$ ratio does not flag the AGN in NGC 3256 and instead inverts the order. S is a low-luminosity AGN, so this is a counter-example to a general AGN signature, not a statement about all AGN.

## 4 Discussion

**Significance.** Published ratios for NGC 3256 are blends in a beam of 1.7$''$ in which N and S differ by 1.5 and the lack of enhancement at the AGN was noted (Harada et al. 2018). We add the within-galaxy comparison at 50 pc with the same array and tuning. What a reader should now believe is the following. (i) At a given physical scale the starburst nucleus of NGC 3256 has HCN/HCO$^+$ about twice that of its AGN, with $C=2.1\pm0.2$ at 50–120 pc, although the core value is subject to a one-sided missing-flux uncertainty. (ii) The difference is a deficit of HCN, HNC and CN relative to HCO$^+$ at S, not an HCN-specific X-ray or infrared effect. (iii) The N nucleus has an HCN-rich core of about 100 pc, in which the ratio approaches 1; this core is an ordinary luminous-infrared-galaxy nucleus by the standards of Nishimura et al. (2024). (iv) No CON-level vibrational emission is present. The consequence for the diagnostic is that its sign cannot be assumed: a high ratio at the AGN, as found in sample studies (Izumi et al. 2016; Imanishi et al. 2016), is not produced by AGN activity at S, and in the same galaxy the AGN nucleus has the lowest ratio. This agrees with the lack of trend with AGN luminosity found by Nishimura et al. (2024), and offers a within-galaxy example, at 50 pc, of why such trends may be absent.

**Alternatives tested and excluded.** *Continuum or dust effects*: the ratio does not depend on continuum brightness, and the species ratios are at similar frequencies. *Beam mixing and flux recovery*: the contrast persists in two 12-m configurations to within 0.06 and the changes with beam size are due to core dilution (Section 3.1). *Different HCO$^+$*: HCO$^+$ is equal at the nuclei. *Vibrational pumping*: CON-level pumping is excluded at both nuclei, and for N at the level of the lower CON range (Section 3.2). *XDR chemistry*: CN is not enhanced at S. *Shocks*: CS and SO ratios show no signature in the nuclei. *Excitation alone*: disfavoured at about $2\sigma$ at 2.2$''$; the rotational ladder gives the same contrast in 1–0 and 3–2 at $\ge1.5''$. Disfavoured, but not excluded: about 7% of excitation models with equal abundance reach a 1–0 contrast of 1.35 or more.

**Limits and what would settle them.** (1) The 2.1 core contrast is measured only at J=3–2 with 12-m data at $\le0.6''$, and its missing-flux term is untested (up to $-0.3$). A 12-m+7-m 3–2 image tapered to 1.0$''$ would settle this. (2) The J-independence is established at $\ge1.5''$, where the contrast is about 1.5, and does not constrain excitation in the N core. The 4–3 contrast is 2.3–3.7 and a rise with J is plausible but not established; a one-tuning 4–3 pair at $\le1''$ with short spacings would decide it. (3) Pumping limits exclude CON-level emission only. (4) The cause is not established: abundance, $^{12}$C/$^{13}$C and combinations of temperature and density are not separated; the model preference is about $2\sigma$ and depends on the error model and the $^{12}$C/$^{13}$C prior. An H$^{13}$CN and H$^{13}$CO$^+$ 1–0 measurement about 10 times deeper (uncertainty about 0.003) at about 2$''$, or HN$^{13}$C/HNC, $^{13}$CO/C$^{18}$O in matched beams, would resolve opacity versus abundance. (5) The temperature of the N core is unconstrained; HCN 4–3/3–2 and thermometer lines (CH$_3$CN, CH$_3$OH, H$_2$CO) at $\le0.5''$ would give it. (6) This is one galaxy with an energetically minor AGN, so a calibration requires resolved ratios and $^{13}$C species in a sample of AGN of known luminosity. The distance, the 14% flux scale and image-plane (not uv-plane) matching of beams are further assumptions. Resolved dense gas near S (Appendix B) with ratios of 1–1.4 is consistent with shock-heated outflow gas of the type discussed by Nishimura et al. (2024), but is a single-array, 3–5σ result.

## 5 Conclusions

In NGC 3256 the HCN/HCO$^+$ J=3–2 brightness-temperature ratio is $0.99$ at the starburst nucleus and $0.44$ at the AGN nucleus at 50 pc, a contrast of $2.1\pm0.2$ that falls with increasing beam size because of dilution of an HCN-rich core of about 100 pc at N. The contrast sits in HCN, HNC and CN relative to HCO$^+$, with equal HCO$^+$ brightness at both nuclei; no vibrationally excited HCN is seen at either nucleus, excluding CON-level pumping. Excitation alone is mildly disfavoured and an HCN-family abundance excess at the starburst nucleus is favoured, but opacity and isotope ratios are not separated. A single-line HCN/HCO$^+$ ratio does not identify the AGN here and inverts the expected order, so it should not be used alone as an AGN diagnostic; resolved multi-line and isotopologue data are needed for a calibration.

---

## Appendix A. Tests supporting the main results

**A.1 Flux recovery and array consistency (Section 3.1).** The Band 6 ratio agrees within 9% across three observing blocks of the long-baseline and compact 12-m configurations; flux-recovery and ladder-ratio systematics are $\lesssim10\%$. At 0.5–1.5$''$ the long-baseline configuration loses 15–20% of the flux at S (HCO$^+$ $-16\%$, HCN $-21\%$) but only about 4% of the HCN flux at N; inside 0.5$''$ the arrays agree to 2–8%. The result does not depend on the aperture, window or baseline within the budget of Appendix C.

**A.2 Line wings (Section 3.1).** An apparent rise of HCN/HCO$^+$ in the line wings is beam mixing of the N core with surrounding gas, and the deficit at S is present at all velocities; earlier wing-ratio candidates are not retained.

**A.3 Continuum independence (Section 3.1).** The ratio is not a function of the continuum: at equal HCO$^+$ it is 0.63 at N and 0.43 at S. The free–free H40$\alpha$ recombination line gives a marginal (1.5–3.3σ) comparison, not used in the main text.

**A.4 Saturation and absorption (Section 3.3).** The N core ratio is not tied to saturated lines or absorption (HCN/H$^{13}$CN $>8.3$).

**A.5 Excitation ladders, temperature and density (Sections 3.2, 3.3).** Ladder fits for N and S, including degeneracies between $T_k$ and $n$ and the effect of beam filling, are consistent with the single-zone result of Section 3.3; the claim of a cool S with $T_k\ge80$ K excluded was retracted. Thermometer lines (CO 2–1 and $^{13}$CO cores) give $T_k(\mathrm{S})\ge26$–40 K. A nested-sampling fit gives $T_c(\mathrm{S})=17\pm1$ K and an unconstrained $T_c(\mathrm{N})$; earlier SMC values were retracted. Two-zone models and a cold optically thick branch are disfavoured but depend on priors; H$^{13}$CN 4–3 is uninformative except at $n\ge10^6$ cm$^{-3}$ ($\alpha<0.10$ at N, $<0.14$ at S). The dependence of the N core's HCN-rich, HCO$^+$-poor nature on the choice of CS as reference was revised and does not change the conclusion.

**A.6 4–3 vs 3–2 (Section 3.2).** Matching HCN 4–3 and HCO$^+$ 4–3 (two programmes) to 3–2 at one beam, the 4–3/3–2 $T_b$ ratios at 1.0$''$ (raw) are HCN $0.646\pm0.011$ (N), $0.471\pm0.034$ (S); HCO$^+$ $0.663\pm0.013$ (N), $0.838\pm0.020$ (S). Two bounding corrections for missing 4–3 flux (a uniform pedestal; and a fractional correction as if all lost flux sat in the core) give contrasts of 3.34 and 2.26 at 1.0$''$ and 3.26 and 2.37 at 1.5$''$ (raw: 3.67, 3.85), against 3–2 contrasts of 2.12 and 2.03. The ratio of 4–3 and 3–2 contrasts is 1.73, 1.58 and 1.07 (raw, pedestal, fractional) at 1.0$''$. The recoveries rely on a single-tuning 4.5$''$ observation with an unexplained 1.5$\times$ absolute flux offset; the 1.6 blended contrast and the uv-plane two-component fits (1.4–1.9) are not core measurements. The compact-to-extended behaviour of S HCN (4–3 flux in the high-resolution data is 0.61 of the single tuning) is the crux. The continuum edge ramp biases the HCN 4–3 measurement by 0–10%. An earlier N/S contrast of 1.26 in the 4–3/3–2 line ratio itself is not robust to short-spacing recovery (0.97–1.24).

**A.7 Isotopologues and the joint fit (Section 3.3).** Bayes factors (abundance-only over excitation-only, natural log), which are always positive:

| prior | errors halved | default | errors doubled |
|---|---|---|---|
| $^{12}$C/$^{13}$C 40–100 | 6.70 (BF 809) | 3.71 (41) | 1.92 (6.8) |
| 40–100, physical | 6.68 (794) | 3.56 (35) | 1.78 (6.0) |
| 20–200 | 6.40 (603) | 3.27 (26) | 1.43 (4.2) |
| 20–200, physical | 6.46 (642) | 3.08 (22) | 1.37 (3.9) |

The "physical" priors are $\log n\ge3.8$, size $\le100$ pc and $X(\mathrm{HCO^+})\le10^{-8}$. Best-fit 1–0 contrasts are excitation-only 1.09, column-only 1.13, abundance-only 1.33 and full 1.29 (observed 1.45); the corresponding 3–2 contrasts are 1.23, 1.35, 1.57 and 1.53 (observed 1.50). The abundance ratio of N to S is $1.58\ [1.26,2.5]$ (X 40–100) with $P(>1)=0.95$; for $^{12}$C/$^{13}$C 20–200 it is $4.0\ [3.2,5.0]$ at N and $2.5\ [2.0,3.2]$ at S. With a 384-variant budget, H$^{13}$CN/HCN is $0.0221\pm0.0075$ at N and $0.0014\pm0.0098$ ($<0.031$, 3σ) at S, and H$^{13}$CO$^+$/HCO$^+$ is $0.0284\pm0.0058$ and $0.0128\pm0.0064$; the double ratio $(^{13}/^{12})_N/(^{13}/^{12})_S$ is 5.8 [2.1–13.5] for HCN (95% lower bound 1.2) and 3.1 [1.7–7.9] for HCO$^+$, both above 1. The differences extend to HCO$^+$ and so are not HCN-specific.

**A.8 Pumping and background (Sections 3.2, 3.3).** Continuum and dust background limits ($T_\mathrm{bg}\lesssim3.2$ K) disfavour radiative pumping of the envelope. Limits on $L_\mathrm{vib}$ and $\Sigma_\mathrm{vib}$ at coarser beams are given in Table C3.

**A.9 Superseded estimates.** An older estimate that only about 40% of the log contrast is J-dependent, and an earlier $D=1.10$, used flux-incomplete data and are superseded by Section 3.2.

## Appendix B. Side findings

**B.1 Dense gas in the S outflow.** A blue lobe and a red ridge near S show HCN/HCO$^+\approx1$–1.4 in a single array at 3–5σ, with a dense-gas mass of about $10^7\,M_\odot$ and a rate of 15–60 $M_\odot$ yr$^{-1}$. This is not established in a second array and is absent at J=1–0; we note it as possible shocked, HCN-rich outflow gas, consistent with Nishimura et al. (2024).

**B.2 Methanol.** CH$_3$OH is detected, extended and cold; an enhancement in the bridge was not reproduced.

**B.3 Disk and clump gas.** The HCN/HCO$^+$ ratio in disk and clump gas declines with J (about 11% at 2–1 and 6–10% at 1–0; 3–2 not measurable).

## Appendix C. Measurement conventions and tables

**Conventions.** $T_b$ is the Rayleigh–Jeans brightness temperature with the $\nu^2$ factor applied. The velocity window is $\pm350$ km s$^{-1}$ about 2780 km s$^{-1}$ (N 2785, S 2750 in some tables), except in the fits of Section 3.3 ($\pm150$ km s$^{-1}$). The baseline is a constant (per-pixel median of line-free channels at $<-450$ and $>+600$ km s$^{-1}$); a linear baseline was used for the 4–3, CN and species tables. Apertures have radius equal to the beam FWHM centred on the HCO$^+$ 3–2 peak (N: 156.96346°, −43.90388°; S: 156.96339°, −43.90533°). Noise is from line-free channels. The 0.26$''$ beam is the common circular convolution of a native 0.246$''\times0.232''$ beam. 1$''$ = 194 pc ($D=40$ Mpc).

**Table C1. N/S contrast of HCN/HCO$^+$ J=3–2 by beam.**

| beam | aperture (long-baseline) | stat | single pixel | peak channel | compact 12-m | 12-m+7-m ($\pm350/\pm150$) | missing-flux term |
|---|---|---|---|---|---|---|---|
| 0.26$''$ (50 pc) | 2.26 | 0.06 | 2.09 | 2.13 | – | – | ≤0.06 assumed |
| 0.28$''$ (55 pc) | 2.32 | 0.06 | 2.14 | 2.24 | – | – | ≤0.06 assumed |
| 0.62$''$ (120 pc) | 2.12 | 0.04 | 2.24 | 2.04 | – | – | ≤0.06 assumed |
| 1.0$''$ (194 pc) | 1.99 | 0.04 | 2.16 | 1.84 | 1.94 | – | −0.06 |
| 1.29$''$ (250 pc) | 1.90 | 0.04 | 2.10 | 1.83 | 1.84 | – | −0.06 |
| 1.5$''$ (291 pc) | 1.83 | 0.04 | 2.04 | 1.83 | 1.81 | 1.53 / 1.63 | −0.3 / −0.2 |
| 2.2$''$ (427 pc) | 1.64 | 0.03 | 1.88 | 2.00 | 1.65 | 1.44 / 1.50 | −0.20 / −0.14 |

Adopted: $2.1\pm0.2$ (0.26–0.62$''$; statistical 0.05, estimators $\pm0.1$, missing flux $\le0.06$ downward, quadrature 0.12, widened to 0.2 to include the untested term), $2.0\pm0.2$ (1.0$''$), $1.9\pm0.15$ (1.29$''$), $1.5\pm0.1$ (2.2$''$, 12-m+7-m, $\pm150$ km s$^{-1}$; 1.44 at $\pm350$). The flux scale (14%) cancels in $C$. Changing the window from $\pm350$ to $\pm150$ km s$^{-1}$ changes the flux-complete 2.2$''$ contrast by 0.06 and the 1.5$''$ one by 0.10; a blue-only versus two-sided baseline changes $R$ by $\le0.01$. Individual ratios at 0.26$''$ are $R_N=0.99\pm0.14$, $R_S=0.44\pm0.06$ including the flux scale. The 1–0 contrast also depends on aperture (1.3–1.5).

**Table C2. HCN/HCO$^+$ J=3–2 ratio in annuli (long-baseline data, 0.26$''$ beam).**

| annulus ($''$) | N | S |
|---|---|---|
| 0–0.25 | 0.998 ± 0.019 | 0.443 ± 0.014 |
| 0.25–0.5 | 0.877 ± 0.017 | 0.386 ± 0.017 |
| 0.5–1 | 0.648 ± 0.011 | 0.331 ± 0.021 |
| 1–2 | 0.432 ± 0.017 | 0.340 ± 0.033 |
| 2–3 | 0.377 ± 0.054 | 0.446 ± 0.027 |

At 1.0$''$ beam, the long-baseline / compact configuration give N 0.631/0.656 and S 0.353/0.376 at 0.5–1.5$''$. Annulus values are the mean over pixels, errors the rms noise divided by the square root of the number of independent beams; the outer annuli are affected by missing flux at the 0.03–0.05 level.

**Table C3. Vibrationally excited HCN $v_2=1f$ J=3–2 limits (3σ).**

| beam | N $\alpha$ | N $L_\mathrm{vib}$ ($L_\odot$) | N $\Sigma_\mathrm{vib}$ | S $\alpha$ | S $L_\mathrm{vib}$ ($L_\odot$) | S $\Sigma_\mathrm{vib}$ |
|---|---|---|---|---|---|---|
| 48 pc | <0.050 | <29 | <0.017 | <0.100 | <25 | <0.015 |
| 180 pc | <0.020 | <52 | <0.002 | <0.050 | <43 | <0.002 |
| 1.1 kpc | <0.010 | <95 | <0.0005 | <0.020 | <102 | <0.0003 |
| 800 pc, J=4–3 | <0.025 | <207 | – | <0.090 | <310 | – |

$\Sigma_\mathrm{vib}$ in $L_\odot$ pc$^{-2}$. Limits are profile-likelihood, the widest of four models of the HCO$^+$ wing (exponential-linear, exponential-quadratic, power law, constant), with noise rescaled by $\sqrt{\chi^2_\mathrm{red}}$; the vibrational profile is assumed equal to the HCN $v=0$ profile. The conversion $10^3\,L_\odot=1.64\times10^6$ K km s$^{-1}$ pc$^2$ at 267.199 GHz follows Falstad et al. (2021). CONs have $L_\mathrm{vib}=1.2$–$62\times10^3\,L_\odot$ and vib/HCN $L'$ of 0.06–0.18 (ten beams) or 0.41–0.69 (one beam). Falstad et al.'s own limit for N and S is $<65\,L_\odot$.

**Table C4. 4–3 and 3–2 ratios at matched beams** (statistical errors only; the 4–3 pedestal ratio is in parentheses).

| beam | HCN N | HCN S | HCO$^+$ N | HCO$^+$ S |
|---|---|---|---|---|
| 1.0$''$ 4–3/3–2 raw | 0.646 | 0.471 | 0.663 | 0.838 |
| 1.0$''$ (pedestal) | 0.655 | 0.515 | 0.691 | 0.857 |
| 1.5$''$ 4–3/3–2 raw | 0.608 | 0.413 | 0.644 | 0.830 |
| 1.5$''$ (pedestal) | 0.623 | 0.488 | 0.683 | 0.860 |

HCN/HCO$^+$ at 4–3: N 0.814, S 0.222 at 1.0$''$; N 0.718, S 0.186 at 1.5$''$; the single tuning at 4.5$''$ gives 0.515 and 0.317. The recovery of the high-resolution 4–3 data relative to the single tuning at 4.5$''$ is HCN 0.88 (N) and 0.61 (S), HCO$^+$ 0.77 (N) and 0.87 (S).

**Table C5. Contrast per unit 1.3 mm continuum** (N/S of line integral ÷ N/S of continuum; continuum N/S 1.54 at 2.5$''$, 1.08 at 1.1$''$, 1.26 at 1.4$''$).

| species | 1–0 2.5$''$ | 3–2 1.1$''$ | 3–2 1.4$''$ |
|---|---|---|---|
| HCN | 1.43 | 3.08 | 2.56 |
| HNC | 1.75 | 3.01 | 2.37 |
| CN | 1.93 | 3.18 | 2.98 |
| HCO$^+$ | 0.92 | 1.44 | 1.26 |

Errors about 10%; multiply by $T_{d,N}/T_{d,S}$ (0.5–2) for the contrast per unit dust mass. The 1–0 line values are from a different pipeline from the 3–2 values.

**Line selection for Table 1.** CN 1–0 is the sum of both hyperfine groups (113.144 and 113.491 GHz), each integrated over its own $\pm350$ km s$^{-1}$ window; CN 3–2 is the whole $N=3$–2 complex (340.032 and 340.248 GHz). HNC 1–0 and HCN, HCO$^+$ 1–0 come from the same observation; HNC 3–2 and CN 3–2 come from the programmes listed in Section 2, with 12-m data only; no 7-m data are available for these lines, so their N/S contrasts cannot be corrected for missing flux. $\pm250$ km s$^{-1}$ windows change the ratios by $\le0.03$ and the contrasts by $\le0.06$. Values from a second team at 1.1$''$ (HCN/HCO$^+$ 2.1, HNC/HCO$^+$ 2.1, CN/HCO$^+$ 2.2) agree with these within the uncertainties.

**Table C6. Literature placement** (percentile of the 21 resolved ratios of Nishimura et al. 2024 below the value).

| nucleus | 55 pc | 120 pc | 250 pc |
|---|---|---|---|
| N | 0.99 (38%) | 0.84 (38%) | 0.69 (19%) |
| S | 0.43 (0%) | 0.40 (0%) | 0.37 (0%) |

Restricting to literature nuclei within $\pm50\%$ of our scale gives 25%, 50% and 22% for N and 0% for S. The 6.5$''$ blends (0.43 and 0.39) lie below all 21 resolved values and at the 5th percentile of the 20 global ratios.

## Appendix D. Provenance

Post identifiers refer to the discussion forum of the two teams; snapshot names are the archived script directories.

| result | posts | snapshot scripts and data |
|---|---|---|
| 3.1 contrast, size, profile | p0031, p0005, p0004, p0035, p0026, p0029, p0003, p0008, p0010, p0013 | team-b-s19 (`out/q26/s171_summary.txt`); error budget `r07_budget.py` (`s170_*.json`, `s400_*.json`); annuli `r09_profile.py`; overview `fig1_overview.py`; `fig2_scale.py` |
| 3.2 J-dependence | p0036, p0060, p0034, p0009, p0059, p0053, p0022, p0027, p0025 | team-b-s21 (`out/q34/s192.json`), team-b-s31 (`out/q56/s345.txt`); 4–3 matching `r10_ladder.py`, `r10_corr.py` |
| 3.2 species | p0056, p0058, p0016 | team-b-s30 (`out/q52/s320.json`, `s321`, `s322.txt`); table `r11_table.py` |
| 3.2 pumping | p0057, p0046, p0047, p0048 | team-b-s29 (`out/s42_hcn_vib.csv`, `out/q51/s310.txt`); `fig6.py` |
| 3.3 cause | p0061, p0062, p0063, p0041, p0042, p0044, p0050, p0054, p0051, p0052, p0037, p0039, p0040, p0007, p0011, p0015, p0021, p0028, p0030 | team-b-s32–s34 (`out/q57*`, `out/q59/s350.txt`); `fig4_r3.py` |
| 3.4 literature | p0055, p0053 | `s300_q46_lit_context.py`, `fig5_lit.py`, `percentiles.txt` |
| Appendix B | p0006, p0014, p0017, p0018, p0024; p0001, p0002, p0019, p0020, p0023; p0032, p0033, p0043, p0045, p0049 | – |

Figures 1–3, 5–7 are the teams' figures; Figure 4 was drawn by the paper team from the tabulated values of Table 1, Appendix A.6 and the matched-beam contrasts of Section 3.2 (`paper/scripts/fig4_family.py`).

## References

Aalto, S., et al. 2015, A&A 584, A42. Butterworth, J., et al. 2022, A&A 667, A131. Falstad, N., et al. 2021, A&A (CON vibrational HCN). Harada, N., et al. 2018, ApJ 855, 49. Imanishi, M., et al. 2016, AJ 152, 218. Izumi, T., et al. 2016, ApJ 818, 42. Michiyama, T., et al. 2018, ApJ 868, 95. Nishimura, Y., et al. 2024, A&A 686, A48.
