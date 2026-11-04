# The HCN/HCO$^+$ J=3–2 ratio rises inward at the starburst nucleus of NGC 3256 but stays low at its low-luminosity AGN

**NGC 3256 paper team**

## Abstract

An elevated HCN/HCO$^+$ ratio in the submillimetre lines is used to flag accreting black holes in dusty nuclei, but the claim is contested. We test it within one galaxy, NGC 3256, whose northern nucleus (N) is a starburst and whose southern nucleus (S) hosts an obscured, low-luminosity AGN, using archival ALMA Band 6 data at 50–250 pc. The HCN/HCO$^+$ J=3–2 brightness-temperature ratio at N rises inward from 0.38 at 2–3$''$ to 1.00 inside 0.25$''$ (49 pc), whereas at S it stays at 0.33–0.44. The N/S contrast is 2.26 (+0.11/−0.15) at 0.3$''$ and 2.0 ± 0.1 at 1$''$ in the image plane. It is not a missing-flux artefact there: both lines are recovered equally by two 12-m configurations (to 2%), and adding short spacings changes the contrast by less than 2% at 1.0–1.5$''$. It is not confirmed in the visibility plane, and nothing inside 0.15$''$ is supported. The contrast falls with beam size to 1.7 (apertures ≤0.75$''$) or 1.45 (beam-sized aperture) at 430 pc, where it is the same at J=1–0 and 3–2. HCN, HNC and CN are all weak against HCO$^+$ at S. No vibrationally excited HCN is seen (surface density below 0.02 $L_\odot$ pc$^{-2}$ for a beam-filling source), which excludes CON-level pumping only marginally at N and not at S. Abundance, opacity, $^{12}$C/$^{13}$C and excitation are not separated. In this galaxy the ratio does not flag the AGN: its nucleus has the lowest ratio.

## 1 Introduction

The HCN/HCO$^+$ ratio of rotational transitions is a long-used diagnostic of the state of dense gas in galaxies. Observations of J=3–2 and 4–3 led several groups to report a "submillimetre HCN enhancement": higher HCN/HCO$^+$ in nuclei hosting active galactic nuclei (AGN) than in starbursts, interpreted as the signature of X-ray or mechanical heating, or of AGN-driven chemistry (Izumi et al. 2016; Imanishi et al. 2016). The claim is contested. Resolved ALMA observations of luminous and ultraluminous infrared galaxies find a wide spread of ratios with no trend with AGN luminosity, and associate the highest ratios with shock-heated outflows (Nishimura et al. 2024). High ratios can also arise from abundance differences (Butterworth et al. 2022) or from infrared pumping of HCN through its bending mode in compact obscured nuclei (CONs; Aalto et al. 2015), in which case vibrationally excited HCN lines should be visible.

NGC 3256 ($D\approx40$ Mpc, $1''\approx194$ pc) has two nuclei separated by about 5$''$ (1 kpc): a starburst-dominated northern nucleus (N) and a southern nucleus that hosts a low-luminosity, obscured AGN (S). At 1.7$''$ (330 pc) Harada et al. (2018) give HCN/HCO$^+$ of 0.96 (N) and 0.65 (S) at J=1–0, i.e. N/S $=1.49\pm0.03$ (recomputed from their Table 13 integrals over $\pm150$ km s$^{-1}$; 1.43 from their peak intensities), and 0.78 and 0.41 at J=3–2, N/S 1.9 (as we transcribed them; window and errors for J=3–2 not re-checked), and noted no enhancement at the AGN. J=1–0 values at about 2$''$ were given by Michiyama et al. (2018), and Ledger et al. (2021) give line ratios in a 2.2$''$ beam. The contrast between the nuclei is therefore known at 330–470 pc. What was missing is its radial structure at 50 pc, where each nucleus is resolved from its surroundings; whether the contrast survives missing flux and beam matching; how it depends on beam size and J; and tests of vibrational pumping. The ALMA archive lists no publication of the Band 6 programme used here.

We ask whether an elevated HCN/HCO$^+$ ratio is a signature of AGN heating in dusty nuclei, using NGC 3256 as a within-galaxy test. We find that the ratio is not high at the AGN: at N it rises inward to about 1 while at S it stays near 0.4; that the contrast is about 2.2 at 50–60 pc and is not produced by missing flux in the image plane; that it is not due to CON-level pumping; and that its cause is not separated. This is one galaxy with an energetically minor AGN, so it is a counter-example to a general signature, not a statement about all AGN.

## 2 Data and methods

**Data.** The main data are ALMA Band 6 observations of NGC 3256 (project 2017.1.00379.S) covering HCN and HCO$^+$ J=3–2 (265.9 and 267.6 GHz rest) in three array configurations: a long-baseline 12-m configuration (Xb4b) with a native beam of $0.25''$ (about 50 pc), a more compact 12-m configuration (Xb4d, $0.94''\times0.74''$), and the 7-m array (Xb4f). Complementary data are J=1–0 of HCN, HCO$^+$, H$^{13}$CN, H$^{13}$CO$^+$, HNC and CN (2015.1.00993.S, 2015.1.00412.S), J=3–2 of HNC and CN (2015.1.00412.S, 2018.1.00493.S) and J=4–3 of HCN and HCO$^+$ (2018.1.00493.S, 2018.1.00223.S, and a 7-m tuning, Xb47). We adopt $v_\mathrm{sys}=2780$ km s$^{-1}$, $D=40$ Mpc and a flux-scale uncertainty of 10% per line (14% per ratio). N and S are measured in the same images, so the flux scale of a given line cancels in every N/S comparison.

**Measurements.** Primary-beam-corrected cubes are converted to Rayleigh–Jeans brightness temperature $T_b$, convolved in the image plane with a flux-conserving kernel to a common circular beam, and baseline-subtracted using line-free channels. Velocity-integrated values are taken over $\pm350$ km s$^{-1}$ about the systemic velocity. The ratio $R=T_b(\mathrm{HCN})/T_b(\mathrm{HCO^+})$ is the ratio of integrals; the contrast is $C=R_N/R_S$. Two estimators are used. The *core annulus* sums the moment-0 maps within $r<0.25''$ of each nucleus (and within $r<1''$ for the larger scale); the *beam aperture* has radius equal to the beam FWHM. The nuclei are 1.45$''$ apart, so at beams $\gtrsim1.5''$ a beam-sized aperture on one nucleus contains the other, which lowers $C$ (Sect. 3.2). Errors are from line-free channels and off-line pseudo-moments (null apertures); baseline, window and estimator variations are treated as systematics (Appendix B).

**Missing flux.** Interferometric images lose extended emission; S is more extended and fainter than N and loses more, which would inflate $C$. We use four tests (details in Appendix A.1). (i) At 1.0$''$ and 1.5$''$ the 12-m configuration Xb4d was re-imaged with and without the 7-m data; adding them changes S by $-2.0\%$ and $-0.8\%$, N by $-0.3\%$ and $-1.5\%$, and $C$ by $+1.7\%$ and $-0.7\%$. The two sets use different weightings, and the uv range is not restricted. (ii) The archive cubes of Xb4b (0.25$''$) and Xb4d, both convolved to 1.05$''$, recover both lines equally: the ratio of the ratios is 0.98–0.99 inside 1$''$ and 0.92–0.93 at 1–1.5$''$, so differential HCN/HCO$^+$ loss is $\le2\%$ inside 1$''$ although the recovery of the Xb4d flux differs between the nuclei (N 0.92–1.02, S 0.74–0.93). (iii) The residual-to-restored flux fraction of the clean images is 2–4% at N and 5–8% at S, a low bias of at most 3.5% on $C$. (iv) Between 1.5$''$ and 2.2$''$ the contrast of the flux-complete 12-m+7-m images (uv-tapered) is lower than that of the image-plane convolved 12-m cubes. The 7-m data themselves account for only $-0.01$ to $-0.11$ of this; the remainder (0.1–0.27, estimator dependent) is a difference between image-plane convolution and uv-taper imaging that we cannot separate into weighting and cleaning (A.1).

Three limits remain. No visibilities of Xb4b could be restored (four of four attempts failed for lack of memory), so there is no uv-plane combination with the 7-m data below 1$''$; the compact-scale deconvolution bias is not bounded; and a pessimistic high-pass model (flux below 60–100 k$\lambda$ removed) gives 2.6–3.1, which shows that the ratio rises when extended gas is filtered but is not a bound. Sub-arcsecond values are image-plane results.

**Modelling.** One-zone non-LTE (RADEX) models are fitted to the 1–0 and 3–2 lines of HCN, HCO$^+$ and their $^{13}$C isotopologues in a matched 2.2$''$ beam (Appendix A.7).

## 3 Results

### 3.1 The N ratio rises inward; the S ratio stays low

**HCN/HCO$^+$ J=3–2 at N rises from 0.38 at 2–3$''$ to 1.00 inside 0.25$''$ (49 pc), whereas at S it stays at 0.33–0.44.** In annuli around each nucleus at 0.26$''$ the N ratio is $1.00\pm0.02$ inside 0.25$''$, $0.88\pm0.02$ at 0.25–0.5$''$, $0.65\pm0.01$ at 0.5–1$''$, $0.43\pm0.02$ at 1–2$''$ and $0.38\pm0.05$ at 2–3$''$ (Fig. 2, Table B2). S is $0.44\pm0.01$, $0.39\pm0.02$, $0.33\pm0.02$, $0.34\pm0.03$ and $0.45\pm0.03$. N exceeds S by 2.3 inside 0.5$''$ and 2.0 at 0.5–1$''$, but only by 1.3 at 1–2$''$, and at 2–3$''$ N is below S. The N profile is a gradual rise over about 200 pc (1$''$); we fit no core size and cannot resolve structure inside the beam. S declines by 25% from 0.44 to 0.33 over the inner 1$''$. Its rise to 0.45 at 2–3$''$ is in the outer low-brightness region, where bridge and extended gas dominate and the long-baseline array loses flux; it is not part of the compact S profile. A second analysis at a 0.30$''$ beam, with different annuli (outer annuli restricted to the half-plane away from the other nucleus) and a different baseline treatment, gives the same inner profile (N 0.985, 0.967, 0.868, 0.648 and S 0.449, 0.418, 0.395, 0.321 in annuli 0–0.15, 0.15–0.25, 0.25–0.5, 0.5–1$''$; Table B2); beyond 1$''$ the two conventions differ and N falls to or below S, so only the inner 1$''$ is used. Nothing inside $r\approx0.15''$ (the native beam) is supported. The annulus errors are statistical; the systematic error from array and flux differences is 0.03 inside 1$''$ and 0.05 outside (Appendix A.1). The inner rise does not depend on the continuum baseline: variants of the baseline windows move the N core ratio by at most 0.05, and a blue-only baseline that includes a weak extra line in the HCO$^+$ band is rejected (A.2).

![Figure 1. Peak brightness temperature $T_b$ (RJ, K) of HCO$^+$ J=3–2 (a) and HCN J=3–2 (b), and the HCN/HCO$^+$ $T_b$ ratio from the moment-0 maps (c), NGC 3256, at a common 0.26$''$ (50 pc) circular beam. Pixels below 3σ are black (a, b) or grey (c). Crosses mark N, S and the bridge; offsets are in arcsec (east to the left) and pc. The circle is the beam. Medians within 0.13$''$ of the markers: the ratio map (c) gives 0.99 (N), 0.48 (S) and 0.51 (bridge); the peak maps (a, b) give HCO$^+$/HCN of 1.67/1.54 K (N) and 2.33/1.05 K (S), whose ratios (0.92, 0.45) differ because (c) is built from moment-0 maps.](figs/fig1_overview.png)

![Figure 2. HCN/HCO$^+$ J=3–2 $T_b$ ratio in annuli around N (orange) and S (blue), long-baseline data at 0.26$''$ (Table B2). Thin bars: statistical errors; thick bars: statistical plus array and flux systematic (0.03 inside 1$''$, 0.05 outside).](figs/fig2_profile.png)

**What this means:** at the resolution where the nuclei separate, the ratio is high only at N, and only in the inner few hundred parsecs.

### 3.2 The N/S contrast is 2.2 at 50–60 pc and 2.0 at 200 pc in the image plane, and falls to 1.5–1.7 at 430 pc

**The N/S contrast of HCN/HCO$^+$ J=3–2 is $2.26\,(+0.11/-0.15)$ at 0.3$''$ ($r<0.25''$, 58 pc) and $2.0\pm0.1$ at 1.05$''$ ($r<1''$, 200 pc) in the image plane; it falls with beam size to 1.7 (apertures $\le0.75''$) or $1.44\pm0.02$ (beam-sized aperture) in flux-complete 12-m+7-m data at 2.2$''$ (430 pc).** The budget at 0.3$''$ is statistical 0.09, baseline and window 0.05, array 0.03, residual-to-restored 0 to $-0.08$ (Appendix B, Table B1); an independent analysis with different baseline windows gives $2.2\pm0.2$ at 0.3$''$ and $2.1\pm0.15$ at 1$''$ (window medians 2.13–2.26). Convolving the 0.25$''$ data in the image plane to larger beams (Fig. 3, Table B1b) the core-annulus contrast is 2.27 (native), 2.21 (0.6$''$), 2.09 (1.05$''$), 2.01 (1.5$''$) and 1.89 (2.2$''$), with statistical errors of 0.07–0.11. Dilution alone therefore explains 45% of the drop from 2.26 to the flux-complete 1.44 for this estimator and 73% for the beam-aperture estimator (2.25 to 1.66). The rest is the uv-taper term of Sect. 2 (about 0.2, 3σ).

![Figure 3. N/S contrast of HCN/HCO$^+$ versus beam FWHM. Blue: J=3–2 (open: 12-m only, filled: 12-m+7-m; circles: core annulus or $r=0.5''$ aperture, triangles: aperture $r=$ beam). Green: J=1–0 (flux-complete). Orange: J=4–3 from the 7-m array alone (circles, matched recipe) and a visibility-plane fit (square, no defined beam). Grey band: Harada et al. (2018), 1.7$''$, J=1–0. Errors are statistical except for the 4–3 fit (statistical and systematic in quadrature). Image-plane points below 1$''$ are not confirmed in the visibility plane. Values in Tables B1, B1b, B1c.](figs/fig3_contrast_beam.png)

At 1.5–2.2$''$ the choice of aperture matters as much as the data. For the flux-complete 12-m+7-m 3–2 images the contrast is 1.79 (peak pixel), 1.73 ($r=0.5''$) and 1.54 ($r=1.5''$) at 1.5$''$, and 1.66, 1.67 and 1.45 at 2.2$''$. The drop at $r=$ beam is aperture overlap, not dilution within the beam: at 2.2$''$ the contrast is flat at 1.66–1.70 for $r\le0.75''$. The 1.44 at 2.2$''$ is therefore the beam-aperture value; the overlap-free flux-complete contrast at 430 pc is about 1.7 (Table B1c). The 1.5$''$ number moves by 0.25–0.30 with the estimator alone, and Harada et al.'s 1.49 at 1.7$''$ (1–0) lies in this range; we do not know their aperture.

**The contrast does not depend on J where it was tested.** In flux-complete data at 1.5–2.2$''$ the ratio of the 3–2 and 1–0 contrasts is $D=0.95$–1.06 for every estimator, including the overlap-free ones (peak pixel at 1.5$''$: $1.00\pm0.07$; at 2.2$''$: $1.01\pm0.05$ and $1.04\pm0.03$ for the two 1–0 cubes; Appendix A.6). With the 7-m array alone, at matched beams, the 4–3/3–2 ratio of contrasts is $1.02\pm0.18$ at 6.5$''$ and $1.07\pm0.17$ at 5$''$ (1.14–1.18 at higher uv cuts, within 1σ). At J=4–3 the contrast is 1.8 (+0.3/−0.2 systematic, ±0.12 statistical) in a visibility-plane fit of the 7-m data above 8 k$\lambda$; it cannot be assigned a size. At the 12-m resolution the 4–3 image-plane contrast is higher (2.6, 2.2, 1.7 at 4.5, 5.5, 7$''$) than for the 7-m array (1.6, 1.5, 1.3), because the 12-m cubes resolve out the extended S HCN; the earlier value of 3.7 at 1$''$ is this array and beam effect, not a core measurement (A.6). So the 4–3 contrast rises with resolution in the same sense as 3–2, but a flux-complete 4–3 contrast below 4$''$ does not exist and J-independence is not tested in the core.

**HCN and HCO$^+$ separately.** At the peak, HCN is 2.0 times lower at S ($1.21\pm0.05$ versus $2.45\pm0.05$ K) while HCO$^+$ is comparable (N/S $0.95\pm0.04$). This depends on the estimator (Table 2): in the integrated 0–0.25$''$ annulus HCN N/S is 2.51 and HCO$^+$ 1.11; in the ring medians of Fig. 1 HCO$^+$ N/S is 0.72 and HCN 1.47, so about 46% of $\ln C$ is in HCO$^+$; at 250 pc the HCO$^+$ peak N/S is $1.37\pm0.03$. HCO$^+$ N/S ranges over 0.72–1.37 and HCN N/S over 1.5–2.5. A ratio to HCO$^+$ cannot say which species is the odd one; per unit 1.3 mm continuum HCO$^+$ is also 1.3–1.4 times weaker at S (Table B5).

**Table 2.** N/S of HCN and HCO$^+$ at 0.26$''$ by estimator (statistical errors, Appendix B).

| estimator | HCO$^+$ N/S | HCN N/S | HCN/HCO$^+$ contrast |
|---|---|---|---|
| peak $T_b$ | $0.95\pm0.04$ | $2.03\pm0.09$ | 2.1 |
| ring median (0.13$''$) | 0.72 | 1.47 | 2.06 |
| integrated, 0–0.25$''$ annulus | $1.11\pm0.02$ | $2.51\pm0.08$ | 2.26 |

**What this means:** inside 200 pc the nuclei differ by a factor of about 2, set in the image plane by data that do not lose flux differentially; the contrast grows toward small scales as expected if a compact N component of ratio near 1 is diluted by extended gas of lower ratio. This remains an image-plane result below 1$''$.

### 3.3 HCN, HNC and CN are all weak relative to HCO$^+$ at S

**In flux-complete J=1–0 data at 2.5$''$ the N/S contrast of the ratio to HCO$^+$ is $1.56\pm0.03$ for HCN, $1.90\pm0.11$ for HNC and $2.10\pm0.04$ for CN (Fig. 4, Table 1), so the three nitrogen species are all weaker at S against HCO$^+$.** The N/S contrast of HNC/HCN is $1.22\pm0.07$ and of CN/HCN $1.35\pm0.03$, so HNC and CN are at least as deficient at S as HCN. The errors are statistical only. The J=3–2 values from 12-m data alone agree in sign but are upper values (Appendix A.9).

![Figure 4. N/S contrast of the ratio of HCN, HNC and CN to HCO$^+$, J=1–0 at 2.5$''$, flux-complete 12-m+7-m. Error bars are statistical.](figs/fig3_species.png)

**Table 1.** $T_b$ ratios (RJ, $\pm350$ km s$^{-1}$), flux-complete J=1–0 at 2.5$''$, one pipeline. Errors statistical.

| ratio | N | S | N/S |
|---|---|---|---|
| HCN/HCO$^+$ | 0.912 | 0.585 | 1.56 ± 0.03 |
| HNC/HCO$^+$ | 0.474 | 0.249 | 1.90 ± 0.11 |
| CN/HCO$^+$ | 1.590 | 0.755 | 2.10 ± 0.04 |
| HNC/HCN | 0.520 | 0.427 | 1.22 ± 0.07 |
| CN/HCN | 1.743 | 1.292 | 1.35 ± 0.03 |

CN is not enhanced at S relative to HCN, which is not the signature of an X-ray-dominated region, although CN is also raised by UV fields and shocks, so this does not exclude one. **What this means:** the difference between the nuclei involves the nitrogen species taken against HCO$^+$, not HCN alone.

### 3.4 No CON-level vibrational HCN emission at either nucleus

**Vibrationally excited HCN $v_2=1f$ J=3–2 (267.199 GHz) is not detected; the beam-filling surface density is $<0.02\,L_\odot$ pc$^{-2}$ at both nuclei, and the line-to-ground-state ratio lies below the CON range marginally at N but not at S.** The 3σ limits in the 48 pc beam on $\alpha=I(v_2{=}1)/I(v{=}0)$, the ratio of velocity-integrated brightness assuming the same profile as HCN $v=0$, are $<0.05$ at N and $<0.10$ at S (Fig. 5), i.e. $L_\mathrm{vib}<29\,L_\odot$ and $<25\,L_\odot$ and $\Sigma_\mathrm{vib}<0.017$ and $<0.015\,L_\odot$ pc$^{-2}$, 59–67 times below the CON criterion of $1\,L_\odot$ pc$^{-2}$ (Falstad et al. 2021) for a source filling the beam. The CON range of the vibrational-to-ground-state ratio, 0.06–0.18, is measured in other beams of different size; our limit at N (0.05) is just below it and that at S (0.10) lies inside it. The comparison is therefore beam dependent and the exclusion at N is marginal. The limit also depends on the model of the HCO$^+$ red wing (the widest of four models is used; Table B3). A buried source smaller than 47 pc diluted in the beam, or moderate pumping ($T_\mathrm{vib}=150$–300 K), is not excluded; this matters at S, which JWST shows to be dust-embedded with an obscured 4.7 $\mu$m continuum (Pereira-Santaella et al. 2024). The limits are ours; we compare them only with the CON criterion and range.

![Figure 5. Spectra of the HCN $v_2=1f$ J=3–2 line (rest 267.199 GHz, green dotted) at N (left) and S (right) in the 0.24$''$ (48 pc) beam. Grey: data (line-free rms 0.42 mJy beam$^{-1}$ per 3.9 MHz channel); blue: model of the HCO$^+$ red wing; orange: wing model plus the HCN $v=0$ profile scaled to the 3σ limit (0.033 and 0.090 for this wing model; 0.05 and 0.10 for the widest of four models).](figs/fig4_hcn_vib.png)

**What this means:** neither nucleus has CON-level vibrational surface brightness for a beam-filling source, so such pumping does not explain the contrast; a compact buried source at S is not excluded.

### 3.5 The cause of the contrast is not separated

**The data do not choose between an abundance difference, opacity, a different $^{12}$C/$^{13}$C ratio and excitation, and no model reproduces the observed contrast.** The fits (Appendix A.7) use the 1–0 and 3–2 lines of HCN, HCO$^+$, H$^{13}$CN and H$^{13}$CO$^+$ in a flux-complete 2.2$''$ beam, where the contrast is about 1.45 (beam aperture); they describe the envelope, not the 50–200 pc contrast. The best model reaches 1.29–1.33, so the fit is not acceptable. Forcing equal abundance ratios at N and S gives best-fit contrasts of only 1.1–1.3, and the fit worsens by $\Delta\chi^2\approx3$–11 (depending on hypothesis and $^{12}$C/$^{13}$C prior) when it is forced to the observed 1.45. The preference of abundance-only over excitation-only ($\Delta\chi^2=5$–8; 1.7–2.8σ depending on one or two degrees of freedom; Bayes factors 2.6–2200 over 48 prior and error-model combinations, 16–44 at the default errors) is weak and depends on the error model.

The direct test of opacity is the contrast of the thin lines. The N/S contrast of H$^{13}$CN/H$^{13}$CO$^+$ in the 2.2$''$ 1–0 data is 2.5–2.7 for three estimators, with 16–84% ranges of 0.7–7.7 for two of them because S H$^{13}$CN is consistent with zero ($0.007\pm0.007$), and $P(>1)=0.73$–0.78 (A.7). If both $^{13}$C lines are thin and $^{12}$C/$^{13}$C is equal at both nuclei, this is the contrast of the column-density ratio, about 2.5, larger than the 1.5 observed in the $^{12}$C lines and compatible with 2.0–2.3 at 50–200 pc; that would mean the $^{12}$C ratio at N is suppressed by HCN opacity. The thin assumption is not demonstrated, and the result is a hint of 1–2σ, not a measurement. The double ratios $(^{13}/^{12})_N/(^{13}/^{12})_S$, which are 1 for an abundance-only difference, have a 95% lower bound below 1 for every estimator (Table A3). The opacity is not constrained: a thin-$^{13}$C inversion gives $\tau_{12}\approx0$–3 for HCN at N, RADEX gives $\tau_{10}(\mathrm{HCN})\approx5$, and the two are not reconciled. The kinetic temperature of S is not measured (CO brightness gives $T_k\gtrsim30$–40 K as a lower limit; ladder fits give 17–30 K in one zone, model dependent) and that of the N core is unconstrained (A.5).

**What this means:** the contrast is not easily produced by rotational excitation alone, but abundance, HCN opacity and isotope ratio remain open.

## 4 Discussion

**Significance.** Published ratios for NGC 3256 are blends in beams of 1.7–2.2$''$ (330–470 pc) in which N exceeds S by about 1.5 at J=1–0 and 1.9 at J=3–2 (Harada et al. 2018), and by similar factors in Ledger et al. (2021). This paper adds, relative to that work: (i) the radial profile at 50 pc, with N rising inward to 1.0 and S flat; (ii) a contrast of 2.2 at 50–60 pc and 2.0 at 200 pc, measured with tests against differential flux loss and beam matching (the published 1.5 is a beam-averaged value; our flux-complete value at that beam is 1.5–1.8 depending on aperture); (iii) the demonstration that the contrast is the same at J=1–0 and 3–2 at 430 pc and carried by HCN, HNC and CN together; (iv) a limit on CON-level vibrational HCN at both nuclei. A reader should now believe: the starburst nucleus has the higher ratio at every scale inside 1$''$, by a factor of about 2 at 50–200 pc in the image plane; the contrast is not a flux-loss artefact at $\ge$1$''$ and is probably not one below; it is carried by the nitrogen species against HCO$^+$; and CON-level pumping does not explain it for beam-filling emission. A reader should not believe a cause, a core size, a visibility-plane confirmation below 1$''$, or a J=4–3 core contrast. The consequence for the diagnostic is that its sign cannot be assumed: the AGN nucleus has the lowest ratio in this galaxy. Relative to other nuclei N is ordinary and S is low (below), so the contrast is carried by S being low; whether the AGN lowers the ratio at S, or S is simply a different environment, one galaxy cannot say.

**Placement.** Nishimura et al. (2024) use $L'_\mathrm{HCN}/L'_\mathrm{HCO^+}$ (equal to the $T_b$ ratio within 1.3%) in an aperture the size of the synthesized beam (about 0.3$''$) and a velocity window chosen per galaxy ($\pm150$ to $\pm500$ km s$^{-1}$); ours uses $\pm350$ km s$^{-1}$. NGC 3256 is not in their sample. To match the aperture diameter to one beam we use the ring median within 0.13$''$ (Fig. 1): N 0.99 and S 0.48. Their 21 resolved ratios run from 0.45 to 2.67 (median 1.13). N lies at the 38th percentile, and S is among the lowest, tied with the minimum (0.445, ESO 148-IG002; S is 0.43–0.48 by estimator). The sample has one AGN-dominated nucleus (IRAS F05189-2524, 1.63) and no trend with AGN luminosity, so it cannot test the AGN order; the percentiles carry no uncertainty. The comparison is indicative only (Fig. 6a, Table B6). At J=4–3, Izumi et al. (2016) give high-resolution ranges of 0.24–1.03 for starbursts (including two single-dish points, excluding LMC clumps at 0.11–0.21) and 1.11–3.19 for AGN. N is 0.515 and S 0.317 at 4.5$''$ in the 7-m data, both in the starburst range. In 12-m data at $\le1''$ S is 0.22 (+0.07/−0.04; window and baseline spread 0.14–0.24) and N 0.81 at 1.0$''$, but those cubes recover only 61% of the S HCN flux, so S is biased low (Fig. 6b, Table B4).

![Figure 6. HCN/HCO$^+$ in context. (a) J=3–2 ratio versus spatial scale: circles are the 21 resolved nuclei of Nishimura et al. (2024) coloured by AGN bolometric fraction (open: not listed); stars are N (orange) and S (blue) at 50, 120 and 250 pc (ring median within 0.13$''$ at 50 pc, aperture means beyond). Definitions differ (see text). (b) J=4–3 ratios at 4.5$''$ of N and S, with the S range from the 1.0$''$ high-resolution cubes, on the high-resolution starburst and AGN ranges of Izumi et al. (2016).](figs/fig5_literature.png)

**Comparison with Ledger et al. (2021).** In a 2.2$''$ beam Ledger et al. find HCN 3–2/1–0 of $0.39\pm0.04$ at N and 0.3–0.4 at S, HCN 4–3/3–2 of $0.48\pm0.07$ (N) and $0.46\pm0.07$ (S), and CN/HCN at 1–0 of $1.7\pm0.1$ at N and about 1.1 at S. Our HCN 3–2/1–0 is equal at the two nuclei (S/N $0.93\pm0.05$, equal within about 7%), and our CN/HCN is 1.74 (N) and 1.29 (S) at 2.5$''$, so the CN excess at N is similar. Ledger et al. attribute the variation within the galaxy to excitation; the equal 3–2/1–0 here does not require an N/S excitation difference.

**Alternatives and their status.** *Missing flux*: excluded as the origin of the contrast in the image plane at $\ge$1$''$ and unlikely below; not tested in the visibility plane below 1$''$. *Beam mixing*: lowers the contrast at large beams (45–73% of the drop to 430 pc), cannot raise it. *Dust or continuum*: dust cancels to first order; dependence on continuum brightness within a nucleus was not measured. *Different HCO$^+$*: not at the peak (N/S $0.95\pm0.04$), but 0.72–1.37 in other estimators. *Vibrational pumping*: no CON-level surface brightness; marginal at N, open at S. *XDR chemistry*: CN is not enhanced at S, which argues against but does not exclude it. *Excitation alone*: weakly disfavoured at 2.2$''$. *Abundance, opacity, $^{12}$C/$^{13}$C*: not separated.

**Limits and what would settle them.** (1) The sub-arcsecond contrast is image-plane only; a restore of the long-baseline visibilities on a host with more memory, followed by a 12-m+7-m combination at 0.3–1$''$, would settle it. The unexplained difference between image-plane convolution and uv-tapered images at 1.5$''$ (0.1–0.27) should be understood in the same run. (2) J-independence is established at 1.5–6.5$''$ only; a 4–3 pair at $\le1''$ with short spacings would test the core. (3) Pumping limits exclude CON-level surface brightness for a beam-filling source only. (4) The cause is not established; H$^{13}$CN and H$^{13}$CO$^+$ 1–0 about ten times deeper at about 2$''$, or HN$^{13}$C/HNC and $^{13}$CO/C$^{18}$O in matched beams, would separate opacity from abundance. (5) The temperature of the N core is unconstrained; HCN 4–3/3–2 at $\le0.5''$ would give it. (6) One galaxy with an energetically minor AGN: calibration needs resolved ratios and $^{13}$C species in a sample of AGN of known luminosity.

## 5 Conclusions

In NGC 3256 the HCN/HCO$^+$ J=3–2 brightness-temperature ratio at the starburst nucleus rises inward from 0.38 at 2–3$''$ to 1.00 inside 0.25$''$, while at the low-luminosity AGN it stays at 0.33–0.44. The N/S contrast is 2.26 (+0.11/−0.15) at 0.3$''$ and 2.0 ± 0.1 at 1$''$ in the image plane, where tests show no differential flux loss above 2%; it is not confirmed in the visibility plane and not measured inside 0.15$''$. It falls to 1.5–1.7 at 430 pc, partly by beam dilution, and is the same at J=1–0 and 3–2 there. HCN, HNC and CN are all weak against HCO$^+$ at S. No CON-level vibrational HCN surface brightness is seen at either nucleus. Abundance, opacity, $^{12}$C/$^{13}$C and excitation are not separated. In this one galaxy the HCN/HCO$^+$ ratio does not flag the AGN; a calibration needs resolved multi-line and isotopologue data in a sample of AGN.

---

## Appendix A. Tests supporting the main results

**A.1 Flux recovery and array consistency (Sections 2, 3.1, 3.2).** *Xb4d with and without the 7-m data* (peak pixel, 12-m: natural weighting with 0.5$''$ or 1.2$''$ taper; 12-m+7-m: Briggs robust 0.5 with 0.3$''$ or 1.2$''$ taper; same window, baseline and circular beam; errors rms/peak, conservative):

| beam | data | N ratio | S ratio | N/S |
|---|---|---|---|---|
| 1.0$''$ | 12-m | $0.867\pm0.050$ | $0.427\pm0.041$ | 2.03 |
| 1.0$''$ | 12-m+7-m | $0.864\pm0.038$ | $0.419\pm0.030$ | 2.06 |
| 1.5$''$ | 12-m | $0.792\pm0.048$ | $0.411\pm0.045$ | 1.93 |
| 1.5$''$ | 12-m+7-m | $0.780\pm0.038$ | $0.408\pm0.034$ | 1.91 |

The S ratio changes by $-2.0\%$ (1.0$''$) and $-0.8\%$ (1.5$''$) and N by $-0.3\%$ and $-1.5\%$; the 9–11% shifts in the same outputs come from the baseline and continuum choice (N 0.867 to 0.966 at 1.0$''$), and the two 12-m execution blocks differ by S 0.387 versus 0.428 (N/S 1.98 versus 2.24), so baseline and execution-block scatter (S $\pm5\%$, N $\pm5$–10%) exceeds the array difference. A second 12-m+7-m imaging of the HCN cube diverged (N/S 0.53/0.35) and is excluded. *Xb4b versus Xb4d* at a common 1.05$''$ beam: recovery N 0.92–1.02, S 0.74–0.93 inside 1$''$; ratio of ratios 0.98–0.99 inside 1$''$ and 0.92–0.93 at 1–1.5$''$; at 0–0.15$''$ the Xb4d recovery of 1.5 is a pixel-size artefact. In the image-plane long-baseline data, at 0.5–1.5$''$ the loss at S is 15–20% (HCO$^+$ $-16\%$, HCN $-21\%$) against about 4% of the HCN flux at N; inside 0.5$''$ the arrays agree to 2–8%. *Residual-to-restored fraction* of the clean images: N 0.02–0.04, S 0.05–0.08, giving a low bias of $\le3.5\%$ on $C$ (measured for the cores at about 1.3$''$ only). *Not tested*: Xb4b visibilities (not restorable), the Xb4b clean depth, any JvM-type correction, and deconvolution bias at 0.3$''$. A Fourier high-pass model of the data (uv cut 60–100 k$\lambda$) gives 2.6–3.1 and is pessimistic. In the visibility plane, a six-component model of Xb4d+7-m gives N/S 2.20 at 0.6$''$, 2.04 at 1.0$''$ and 1.90 at 1.5$''$ (image: 2.26, 2.14, 2.03), but extrapolates to 2.44 at 0.4$''$ and 3.16 at 0.25$''$ with a model-form spread of 0.98–1.5 on N(0.25$''$); it is not used below 0.6$''$. The N HCN/HCO$^+$ amplitude ratio in annuli of uv distance is 0.59–0.92 at 3–90 k$\lambda$ but $1.29\pm0.08$ at 90–130 and $1.49\pm0.28$ at 130–200 k$\lambda$; this is unexplained (non-Gaussian N shape or annulus blending are candidates) and enters only as a $+0.1$ allowance on the 0.3$''$ contrast. *The 1.5$''$ estimator spread* (Table B1c): at 1.5$''$ the 12-m+7-m, 12-m-only and Xb4b-convolved images give 1.79, 1.90 and 2.03 at the peak pixel and 1.54, 1.59 and 1.83 at $r=$ beam. The 7-m data account for $-0.11$ (peak), $-0.10$ ($r=0.5''$) and $-0.05$ ($r=$ beam) at 1.5$''$, and $-0.03$, $-0.03$, $-0.01$ at 2.2$''$. The Xb4d image convolved to 1.5$''$ gives 1.81 at $r=$ beam against 1.59 for the 12-m uv taper of the same configuration (0.23); at the peak pixel the re-imaged Xb4d 12-m (1.93) and the uv taper (1.90) agree, so the term depends on the estimator and we do not know whether weighting or cleaning causes it. Earlier statements that adding 7-m data lowers the contrast by 0.14–0.30 at 1.5–2.2$''$ compared the 7-m-including taper with image-plane convolved cubes and so mixed this term with the short-spacing effect; they are replaced by the numbers above.

**A.2 Line wings, velocity dependence and baseline (Sections 3.1, 3.2).** An apparent rise of HCN/HCO$^+$ in the line wings is beam mixing of the N core with surrounding gas. In 50 km s$^{-1}$ bins at a 1.05$''$ beam (Xb4d, statistical errors) the ratio is, for $v-v_\mathrm{sys}$ from $-150$ to $+100$ km s$^{-1}$, N $0.79$–$0.84$ ($\pm0.01$, $\pm0.06$ in the outermost bin) and S $0.36$–$0.41$ ($\pm0.01$), i.e. flat within a few per cent; the red bin at 100–150 km s$^{-1}$ gives N $1.27\pm0.08$ and S $0.45\pm0.06$, where HCO$^+$ is weak and the HCN $v_2=1f$ line lies, so it is not interpreted. Bins beyond $\pm150$ km s$^{-1}$ are limits. *Baseline*: the HCO$^+$ band holds a weak extra line at $-1100$ to $-900$ km s$^{-1}$. A blue-only baseline in the visibility fit that contains it is rejected ($\chi^2_\nu=47.7$ versus 5.5; N(0.4$''$) 1.12 versus 0.96). In the image plane the N core ratio ($r<0.15''$, 0.30$''$) is 0.97 (adopted), 0.98 (order 0), 0.97 (extra window), 0.99 (blue-only, which excludes the extra line), 0.94 (red-only) and 0.99 (pipeline subtraction only); S is 0.446–0.461 for all. The residual baseline at line centre is $\le0.06$ mJy beam$^{-1}$ against line peaks of about 10 (N) and 5–11 (S) mJy beam$^{-1}$; core continuum peaks are 0.58–0.61 (N) and 1.87–1.89 (S) mJy beam$^{-1}$ at 0.24$''$.

**A.3 Continuum (Section 3.2).** At equal HCO$^+$ the ratio is 0.63 at N and 0.43 at S; this shows that the nuclei differ at equal HCO$^+$, but not that the ratio is independent of continuum brightness within a nucleus, which was not measured. The H40$\alpha$ comparison is marginal (1.5–3.3σ) and not used.

**A.4 Saturation (Section 3.5).** The limit HCN/H$^{13}$CN $>9$ at N and $>6.7$ at S (3σ, from H$^{13}$CN 4–3, depending on the CO 3–2 wing model) is far below any plausible $^{12}$C/$^{13}$C and says nothing about saturation; opacity is discussed in Section 3.5.

**A.5 Excitation, temperature and density (Section 3.5).** Ladder fits for N and S, including degeneracies between $T_k$ and $n$ and beam filling, are consistent with the one-zone result. The kinetic temperature of S is not measured. CO(2–1) brightness at 0.2$''$ (peak 31 K, 20 K velocity-matched) gives $T_k\gtrsim30$–40 K as a lower limit if the CO is thick and co-spatial (filling factor 0.6–0.9; low confidence). A one-zone nested-sampling fit gives $T_c(\mathrm{S})=17\pm1$ K, and an earlier fit at 1$''$ allowed 25–80 K once the deblended 6$''$ 4–3/3–2 values were varied; these depend on the model, envelope and band-scale priors and are not reconciled with the CO limit unless the dense gas is cooler or sub-thermal relative to the CO. A shared extended (10$''$) slab prefers 80–130 K, which does not apply to the S core. $T_c(\mathrm{N})$ is unconstrained (12–95 K). Two-zone models and a cold optically thick branch are disfavoured but depend on priors; H$^{13}$CN 4–3 is uninformative except at $n\ge10^6$ cm$^{-3}$. HCN 4–3/3–2 at $\le0.5''$ (0.34 cold versus 0.60 warm) or $^{13}$CO/C$^{18}$O would decide.

**A.6 J-dependence: 1–0, 3–2 and 4–3 (Section 3.2).** *1–0 versus 3–2, flux-complete* ($\pm350$ km s$^{-1}$; Table B1c). For the 2.2$''$ 1–0 cubes of both programmes and the 3–2 12-m+7-m images, $D=C_{3-2}/C_{1-0}$ is 1.01–1.06 at the peak pixel and $r\le0.75''$ and 0.98–1.00 at $r=$ beam (errors 0.02–0.05); at 1.5$''$ it is 1.00 ± 0.07 (peak), 0.98 ± 0.06 ($r=0.5''$) and 0.95 ± 0.04 ($r=$ beam), and with 12-m-only 3–2, 1.06 ± 0.07, 1.03 ± 0.06 and 0.98 ± 0.04. All are consistent with 1 (adopted: equal within about 7%). The 3–2/1–0 HCN ratio at S/N is $0.93\pm0.05$. At 1.7$''$ the 12-m-only 3–2 contrast is 1.98 (N 0.748, S 0.378) against 1.49 at 1–0 from Harada et al.; the apertures differ and we do not use this pair. None of this tests the core. *ACA only, matched recipe* (same circular beam, positions, window, baseline; uv cuts scaled by frequency; errors from null apertures): N/S at 3–2 and 4–3 are $1.23\pm0.16$ and $1.25\pm0.14$ at 6.5$''$ (ratio $1.02\pm0.18$), $1.36\pm0.18$ and $1.45\pm0.13$ at 5$''$ (1.07 ± 0.17), $1.40\pm0.19$ and $1.60\pm0.17$ with cuts of 11.3 and 15 k$\lambda$ (1.14 ± 0.19), $1.56\pm0.33$ and $1.84\pm0.26$ with 15 and 20 k$\lambda$ (1.18 ± 0.30). The HCN contrast is $1.23$/$1.48$ (3–2) and $1.38$/$1.76$ (4–3) at 6.5$''$/5$''$ and HCO$^+$ $1.00$/$1.09$ and $1.10$/$1.21$. Deblending N and S with three two-component models changes the 3–2 contrast to 1.34–1.72 and the 4–3 to 1.24–1.34 at 6.5$''$ and the ratio to 0.78–0.93, so blending does not create a 4–3 excess. Earlier 3–2 values of 1.07 at 6.5$''$ and 4–3 values of 1.62, 1.49, 1.32 at 4.5, 5.5, 7$''$ came from other apertures and baselines and are not paired with these. *Visibility plane.* A two-Gaussian and a three-component fit of the eight usable 7-m 4–3 execution blocks above 8 k$\lambda$ gives a contrast of 1.8 (bootstrap $\pm0.12$, systematic $+0.3/-0.2$ from model and baseline; 16 cells 1.60–2.25, median 1.71; half-widths 250/500 km s$^{-1}$ give 1.72/1.54; including a ninth block, whose HCN N/S inversion is unexplained, gives 1.59). The fitted size changes nothing (FWHM 0.3–4.8$''$: 1.67, 1.66, 1.64, 1.60, 1.73), so no scale can be attached, and the HCN flux of S differs by 70% between near-identical blocks (an uncertainty of 0.3–0.4 in the contrast). *12-m 4–3.* The 12-m cubes (HCN and HCO$^+$ from two programmes) give 2.64, 2.18, 1.72 at 4.5, 5.5, 7$''$ against 1.62, 1.49, 1.32 for the 7-m data, driven by HCN N/S (2.98 versus 2.06 at 4.5$''$); the 12-m S HCN recovery is 0.61 of the 7-m value (HCO$^+$ 0.87). The earlier 3.7 at 1$''$ and the 1.5$\times$ offset in absolute flux between the two programmes and the 7-m tuning were traced to this loss and to zero-filled edge channels entering the baseline fit (18 and 72 channels); with those removed the offset is 1.06–1.15, and the offset cancels in N/S. Matching HCN 4–3 and HCO$^+$ 4–3 to 3–2 at one beam, the 4–3/3–2 $T_b$ ratios at 1.0$''$ (raw) are HCN $0.646\pm0.011$ (N), $0.471\pm0.034$ (S); HCO$^+$ $0.663\pm0.013$ (N), $0.838\pm0.020$ (S). Bounding corrections for missing 4–3 flux give 4–3 contrasts of 2.3–3.3 at 1.0$''$, against 3–2 contrasts of 2.1; these are not core measurements. The continuum edge ramp biases the HCN 4–3 measurement by 0–10%.

**A.7 Isotopologues and the joint fit (Section 3.5).** The $^{12}$C/$^{13}$C ratio is a free parameter in the fits (grids of 40–100 and 20–200, the first spanning starburst-nucleus values measured in other galaxies by Tang et al. 2019); it is not measured here and no value is adopted for N. A thin-line inversion of the 1–0 ratios gives $^{12}$C/$^{13}$C of 45 (+27/−14) from HCN and 36 (+8/−6) from HCO$^+$ at N; if the $^{12}$C lines are thick these are lower limits and imply $\tau_{12}\approx1$–3 (HCN) and 1–4 (HCO$^+$) for $^{12}$C/$^{13}$C of 60–135. The RADEX posterior at N gives $\tau_{10}(\mathrm{HCN})=4.6\,[3.0,8.7]$ and $\tau_{10}(\mathrm{HCO^+})=0.5\,[0.2,4.7]$, and $\tau_{10}(\mathrm{HCN})\approx5$–12 at both nuclei. A posterior for $X_S/X_N=(^{12}\mathrm{C}/^{13}\mathrm{C})_S/(^{12}\mathrm{C}/^{13}\mathrm{C})_N$ has median 0.98 and 95% interval 0.45–2.55 (medians 0.72–1.39 and 97.5% upper limits 1.9–3.4 depending on estimator and prior). The $\Delta\chi^2$ between abundance-only and excitation-only fits is 5–8 (the best fits of the two hypotheses, 7.7 and 5.0 above the unconstrained fit for the $^{12}$C/$^{13}$C 40–100 prior); as a Gaussian-equivalent significance this is 2.2–2.8σ for one degree of freedom and 1.7–2.4σ for two, and no single number is quoted. When the model contrast is forced to the observed 1.45 the profile is 11.0 (excitation-only) and 6.0 (column-only) above its own minimum for $^{12}$C/$^{13}$C 40–100 and 3.0 and 4.5 for 20–200. Fig. A1 shows the fits. Bayes factors (abundance-only over excitation-only; natural log in the table) are always positive but strongly error-model dependent; over all 48 cells (3 error models × 2 $^{12}$C/$^{13}$C priors × 8 density and path-length priors) they run from 2.6 to 2200 (halved errors 371–2200, default 16–44, doubled 2.6–6.9). Four of the cells:

| prior | errors halved | default | errors doubled |
|---|---|---|---|
| $^{12}$C/$^{13}$C 40–100 | 6.70 (BF 809) | 3.71 (41) | 1.92 (6.8) |
| 40–100, physical | 6.68 (794) | 3.56 (35) | 1.78 (6.0) |
| 20–200 | 6.40 (603) | 3.27 (26) | 1.43 (4.2) |
| 20–200, physical | 6.46 (642) | 3.08 (22) | 1.37 (3.9) |

The "physical" priors are $\log n\ge3.8$, size $\le100$ pc and $X(\mathrm{HCO^+})\le10^{-8}$. Best-fit 1–0 contrasts are excitation-only 1.09, column-only 1.13, abundance-only 1.33 and full 1.29 (observed 1.45, $\pm150$ km s$^{-1}$, statistical error 0.015); the 3–2 contrasts are 1.23, 1.35, 1.57, 1.53 (observed 1.50; no 3–2 profile was computed). $N(\mathrm{HCN})/N(\mathrm{HCO^+})=3.2\,[2.5,4.0]$ at N and $2.0\,[1.6,2.5]$ at S; moving the ratio at S from 0.4 to 1 needs a factor of 3–4 in abundance ratio.

![Figure A1. $\Delta\chi^2$ relative to the full model as a function of the model 1–0 HCN/HCO$^+$ N/S contrast when the abundance ratio is forced equal at N and S (excitation-only blue, column-density-only orange; solid: $^{12}$C/$^{13}$C 40–100, dashed: 20–200). The observed contrast is shown in magenta. Squares/diamonds: abundance-only and full-model best fits. Right: natural logarithm of the Bayes factor of abundance-only over excitation-only for error models and priors.](figs/figA1_hypotheses.png)

**Table A3.** $^{13}$C/$^{12}$C ratios at N and S (1–0, 2.2$''$ aperture) and double ratio (13/12)$_N$/(13/12)$_S$ by estimator (T template, A aperture, P peak pixel; 1σ errors; MC with Gaussian errors).

| line | est. | N | S | double ratio, nominal | MC median [16–84%] | P(≤1) |
|---|---|---|---|---|---|---|
| H$^{13}$CN/HCN | T | 0.0217±0.0079 | 0.0070±0.0068 | 3.1 | 2.6 [1.3–7.0] | 0.08 |
| H$^{13}$CN/HCN | A | 0.0204±0.0058 | 0.0060±0.0091 | 3.4 | 2.3 [1.1–6.6] | 0.09 |
| H$^{13}$CN/HCN | P | 0.0319±0.0044 | 0.0216±0.0061 | 1.5 | 1.5 [1.1–2.1] | 0.08 |
| H$^{13}$CO$^+$/HCO$^+$ | T | 0.0192±0.0039 | 0.0114±0.0092 | 1.7 | 1.5 [0.9–3.7] | 0.22 |
| H$^{13}$CO$^+$/HCO$^+$ | A | 0.0276±0.0049 | 0.0144±0.0054 | 1.9 | 1.9 [1.3–3.1] | 0.04 |
| H$^{13}$CO$^+$/HCO$^+$ | P | 0.0306±0.0031 | 0.0327±0.0043 | 0.9 | 0.9 [0.8–1.1] | 0.65 |

The MC ratio is ill defined where S is consistent with 0 (S≤0 in 11–26% of draws). Opacity at N for thin $^{13}$C lines (single zone, same excitation), for assumed $^{12}$C/$^{13}$C = 40, 60, 90, 135: $\tau_{12}(\mathrm{HCN})=0$, 0.6, 1.6, 2.8 and $\tau_{12}(\mathrm{HCO^+})=0.23$, 1.2, 2.3, 3.7.

**Thin-line contrast.** $(\mathrm{H^{13}CN/H^{13}CO^+})_N/(\ldots)_S$ in the 2.2$''$ 1–0 beam, $T_b$ scale: template 2.7 (MC median 2.4 [0.75, 7.7]; 24% of draws dropped for non-positive S; $P(>1)=0.78$), aperture 2.6 (1.7 [0.72, 5.4]; $P(>1)=0.73$), peak pixel 2.5 (2.5 [1.8, 3.7]; biased high at S). The $^{12}$C HCN/HCO$^+$ contrast in the same apertures is 1.48 (T, A) and 1.57 (P).

**A.8 Pumping and background (Sections 3.4, 3.5).** Continuum and dust background limits ($T_\mathrm{bg}\lesssim3.2$ K) disfavour radiative pumping of the envelope. Limits at coarser beams are in Table B3.

**A.9 12-m-only J=3–2 species contrasts (Section 3.3).** HCN/HCO$^+$, HNC/HCO$^+$, CN/HCO$^+$ at 3–2: 1.4$''$ $2.14\pm0.05$, $1.84\pm0.14$, $2.56\pm0.09$ (N 0.833, 0.477, 0.893; S 0.389, 0.259, 0.348); 2.5$''$ $1.93\pm0.05$, $2.26\pm0.21$, $2.99\pm0.14$. HNC/HCN is $0.86\pm0.06$ (1.4$''$) and $1.17\pm0.10$ (2.5$''$), CN/HCN 1.2–1.55. These lack 7-m data and are upper values of the contrast; the HNC value at 1.4$''$ is below HCN and is a recovery effect. The 1.4$''$ HCN/HCO$^+$ value is 2.14 in one analysis and 2.01 ± 0.04 in another (N 0.833 versus 0.791, a baseline and window difference); the 1.58 sometimes compared with it is the 1–0 contrast at 2.5$''$ (N 0.906, S 0.575) and is a different line, beam and estimator. Normalised to the 1.3 mm continuum, HCN, HNC and CN are 2.4–3.1 times and HCO$^+$ 1.3–1.4 times weaker at S (Table B5). That factor scales with the unmeasured dust-temperature ratio $T_{d,N}/T_{d,S}$ (0.5–2; line per dust column scales as $T_{d,N}/T_{d,S}$) and is not used. Non-dust continuum at S would make the line deficit smaller, not larger. None of the entries is uv-matched: the 1–0 lines are flux-complete and the continuum and 3–2 lines are 12-m only.

## Appendix B. Measurement conventions and tables

**Conventions.** $T_b$ is the Rayleigh–Jeans brightness temperature with the $\nu^2$ factor applied (moment-0 ratios in the 0.30$''$ analysis lack the factor 1.013, which cancels in N/S). The velocity window is $\pm350$ km s$^{-1}$ about 2780 km s$^{-1}$ (N 2785, S 2750 in some tables), except in the fits of Section 3.5 ($\pm150$ km s$^{-1}$). The baseline is a constant (per-pixel median of line-free channels at $<-450$ and $>+600$ km s$^{-1}$, or 600–1500 km s$^{-1}$ from line centre); a linear baseline was used for the 4–3, CN and species tables. Beam-aperture values are centred on the HCO$^+$ 3–2 peak (N: 156.96346°, −43.90388°; S: 156.96339°, −43.90533°). Noise is from line-free channels or off-line pseudo-moments. The 0.26$''$ beam is the common circular convolution of a native 0.246$''\times0.232''$ beam. 1$''$ = 194 pc ($D=40$ Mpc).

**Table B1. Error budget of the N/S contrast of HCN/HCO$^+$ J=3–2 (core annulus).**

| beam | $C$ | statistical | baseline/window | array (Xb4b/Xb4d) | residual/restored |
|---|---|---|---|---|---|
| 0.3$''$, $r<0.25''$ | 2.26 | 0.09 | 0.05 | 0.03 | 0 / −0.08 |
| 1.05$''$, $r<1''$ | 1.98 | 0.06 | 0.01–0.05 | 0.04 | 0 / −0.08 |

The totals are 2.26 (+0.11/−0.15) and 2.0 ± 0.1. Estimators at 0.26$''$: aperture 2.26, single pixel 2.09, peak channel 2.13. Apertures $r=$ beam at 0.28, 0.62 and 1.29$''$ (long-baseline data): 2.32, 2.11, 1.90 ($\pm0.2$ estimator spread). The flux scale (14%) cancels in $C$. Changing the window from $\pm350$ to $\pm150$ km s$^{-1}$ changes the flux-complete 2.2$''$ contrast by 0.06 and the 1.5$''$ one by 0.10; a blue-only versus two-sided baseline changes $R$ by $\le0.01$. Individual ratios at 0.26$''$ are $R_N=0.99\pm0.14$, $R_S=0.44\pm0.06$ including the flux scale. A second analysis (window medians) gives 2.13 [1.92–2.25] at 0.3$''$, 2.21 at 0.4$''$, 2.26 at 0.6$''$, 2.14 at 1.0$''$ and 2.03 at 1.5$''$.

**Table B1b. N/S contrast versus beam, Xb4b convolved in the image plane** (statistical errors from off-line pseudo-moments).

| beam | core $r<0.25''$ | $r<1''$ | aperture $r=$ beam |
|---|---|---|---|
| native 0.25$''$ | 2.27 ± 0.11 | 2.03 ± 0.10 | 2.23 ± 0.11 |
| 0.3$''$ | 2.26 ± 0.11 | 2.03 ± 0.10 | 2.25 ± 0.13 |
| 0.6$''$ | 2.21 ± 0.11 | 2.02 ± 0.09 | 2.11 ± 0.11 |
| 1.05$''$ | 2.09 ± 0.10 | 1.98 ± 0.08 | 1.98 ± 0.08 |
| 1.5$''$ | 2.01 ± 0.08 | 1.94 ± 0.07 | 1.85 ± 0.07 |
| 2.2$''$ | 1.89 ± 0.07 | 1.85 ± 0.07 | 1.66 ± 0.07 |

**Table B1c. N/S contrast at 1.5$''$ and 2.2$''$ by estimator** (uv-tapered images, $\pm350$ km s$^{-1}$; 1σ from line-free channels; the estimators share pixels, so differences are far more significant than the errors).

| beam | data | peak pixel | $r=0.25''$ | $r=0.5''$ | $r=0.75''$ | $r=$ beam |
|---|---|---|---|---|---|---|
| 1.5$''$ | 3–2 12-m+7-m | 1.79 ± 0.04 | 1.76 | 1.73 | 1.70 | 1.54 ± 0.03 |
| 1.5$''$ | 3–2 12-m only | 1.90 ± 0.04 | 1.87 | 1.83 | 1.78 | 1.59 ± 0.03 |
| 1.5$''$ | 3–2 Xb4b convolved | 2.03 ± 0.05 | 2.00 | 1.98 | 1.95 | 1.83 ± 0.04 |
| 1.5$''$ | 1–0 flux-complete | 1.79 ± 0.12 | 1.79 | 1.77 | 1.74 | 1.62 ± 0.06 |
| 2.2$''$ | 3–2 12-m+7-m | 1.66 ± 0.03 | 1.67 | 1.67 | 1.66 | 1.45 ± 0.02 |
| 2.2$''$ | 3–2 12-m only | 1.69 ± 0.03 | 1.70 | 1.69 | 1.67 | 1.46 ± 0.03 |
| 2.2$''$ | 1–0 flux-complete | 1.60–1.64 ± 0.04–0.07 | 1.60–1.63 | 1.59–1.62 | 1.59–1.60 | 1.47–1.48 ± 0.02–0.04 |

For the 1.5$''$ 12-m+7-m images N is 0.742, 0.740, 0.733, 0.719, 0.661 and S 0.414, 0.420, 0.423, 0.424, 0.429 for the estimators in order. The 1–0 rows are the two flux-complete programmes (the range is over them); $\pm150$ km s$^{-1}$ windows give 1.74 (1.5$''$ peak) and 1.58–1.63 (2.2$''$). Also, flux-complete 1–0 at 2.5$''$: 1.57 (peak), 1.43 ($r=$ beam).

**Table B2. HCN/HCO$^+$ J=3–2 ratio in annuli** (long-baseline data, 0.26$''$ beam; and second analysis at 0.30$''$ with half-plane outer annuli).

| annulus ($''$) | N | S | N/S |
|---|---|---|---|
| 0–0.25 | 0.998 ± 0.019 | 0.443 ± 0.014 | 2.25 |
| 0.25–0.5 | 0.877 ± 0.017 | 0.386 ± 0.017 | 2.27 |
| 0.5–1 | 0.648 ± 0.011 | 0.331 ± 0.021 | 1.96 |
| 1–2 | 0.432 ± 0.017 | 0.340 ± 0.033 | 1.27 |
| 2–3 | 0.377 ± 0.054 | 0.446 ± 0.027 | 0.85 |

Annulus values are the mean over pixels, errors the rms noise divided by the square root of the number of independent beams (statistical only). Systematic terms: 0.03 inside 1$''$ (long-baseline versus compact configuration at a 1.0$''$ beam: N 0.631/0.656 and S 0.353/0.376 at 0.5–1.5$''$) and 0.05 at 1–3$''$, where missing flux affects the outer annuli. Second analysis (0.30$''$; annuli 0–0.15, 0.15–0.25, 0.25–0.5, 0.5–1, 1–1.5, 1.5–2$''$; statistical plus baseline error; inside 1$''$ the half-plane restriction is not applied): N 0.985, 0.967, 0.868, 0.648, 0.342, 0.245 and S 0.449, 0.418, 0.395, 0.321, 0.371, 0.502; N/S 2.19, 2.31, 2.20, 2.02, 0.92, 0.49. Errors on the individual inner ratios are 0.03–0.08; beyond 1$''$ they are 0.06–0.28. The 0–0.15$''$ annulus lies inside the native beam and is unsupported.

**Table B3. Vibrationally excited HCN $v_2=1f$ J=3–2 limits (3σ).**

| beam | N $\alpha$ | N $L_\mathrm{vib}$ ($L_\odot$) | N $\Sigma_\mathrm{vib}$ | S $\alpha$ | S $L_\mathrm{vib}$ ($L_\odot$) | S $\Sigma_\mathrm{vib}$ |
|---|---|---|---|---|---|---|
| 48 pc | <0.050 | <29 | <0.017 | <0.100 | <25 | <0.015 |
| 180 pc | <0.020 | <52 | <0.002 | <0.050 | <43 | <0.002 |
| 1.1 kpc | <0.010 | <95 | <0.0005 | <0.020 | <102 | <0.0003 |
| 800 pc, J=4–3 | <0.025 | <207 | – | <0.090 | <310 | – |

$\Sigma_\mathrm{vib}$ in $L_\odot$ pc$^{-2}$. Limits are profile-likelihood, the widest of four models of the HCO$^+$ wing (exponential-linear, exponential-quadratic, power law, constant), with noise rescaled by $\sqrt{\chi^2_\mathrm{red}}$; the vibrational profile is assumed equal to the HCN $v=0$ profile. The conversion $10^3\,L_\odot=1.64\times10^6$ K km s$^{-1}$ pc$^2$ at 267.199 GHz follows Falstad et al. (2021). CONs have $L_\mathrm{vib}=1.2$–$62\times10^3\,L_\odot$ and vib/HCN $L'$ of 0.06–0.18 (ten sources, beams differ from ours) or 0.41–0.69 (one source).

**Table B4. 4–3 and 3–2 ratios at matched beams** (statistical errors only; the 4–3 pedestal ratio is in parentheses).

| beam | HCN N | HCN S | HCO$^+$ N | HCO$^+$ S |
|---|---|---|---|---|
| 1.0$''$ 4–3/3–2 raw | 0.646 | 0.471 | 0.663 | 0.838 |
| 1.0$''$ (pedestal) | 0.655 | 0.515 | 0.691 | 0.857 |
| 1.5$''$ 4–3/3–2 raw | 0.608 | 0.413 | 0.644 | 0.830 |
| 1.5$''$ (pedestal) | 0.623 | 0.488 | 0.683 | 0.860 |

HCN/HCO$^+$ at 4–3 in 12-m data (fixed window, linear baseline): N 0.81 ± 0.02, S 0.22 ± 0.02 at 1.0$''$ (Monte Carlo 0.80 and 0.23; S range 0.14–0.24 with window and baseline, because S HCO$^+$ has blue and red horns and HCN only a blue one); N 0.67, S 0.18 at 1.9$''$; both lower limits for S because of the HCN loss. In the 7-m data the single tuning at 4.5$''$ gives N 0.515, S 0.317 ($\pm0.011$), and the native-beam peak pixel (4.14$''\times2.70''$) gives 0.52 and 0.29. The 0.29 sometimes quoted for S at 56 pc is an extrapolation, not a measurement. The recovery of the 12-m 4–3 data relative to the single tuning at 4.5$''$ is HCN 0.88 (N) and 0.61 (S), HCO$^+$ 0.77 (N) and 0.87 (S).

**Table B5. Contrast per unit 1.3 mm continuum** (N/S of line integral ÷ N/S of continuum; continuum N/S 1.54 at 2.5$''$, 1.08 at 1.1$''$, 1.26 at 1.4$''$).

| species | 1–0 2.5$''$ | 3–2 1.1$''$ | 3–2 1.4$''$ |
|---|---|---|---|
| HCN | 1.43 | 3.08 | 2.56 |
| HNC | 1.75 | 3.01 | 2.37 |
| CN | 1.93 | 3.18 | 2.98 |
| HCO$^+$ | 0.92 | 1.44 | 1.26 |

Errors about 10%; multiply by $T_{d,N}/T_{d,S}$ (0.5–2) for the contrast per unit dust mass. The 1–0 line values are from a different pipeline from the 3–2 values; a second pipeline gives 2.49, 2.14, 2.98 and 1.16 at 3–2 1.4$''$ (differences up to 10%, same sign).

**Line selection for Table 1.** CN 1–0 is the sum of both hyperfine groups (113.144 and 113.491 GHz), each integrated over its own $\pm350$ km s$^{-1}$ window; CN 3–2 is the whole $N=3$–2 complex (340.032 and 340.248 GHz). HNC 1–0 and HCN, HCO$^+$ 1–0 come from the same observation; HNC 3–2 and CN 3–2 come from the programmes listed in Section 2, with 12-m data only; no 7-m data are available for these lines, so their N/S contrasts cannot be corrected for missing flux. $\pm250$ km s$^{-1}$ windows change the ratios by $\le0.03$ and the contrasts by $\le0.06$.

**Table B6. Literature placement** (J=3–2; percentile of the 21 resolved ratios of Nishimura et al. 2024 below the value; indicative only, definitions differ).

| nucleus | 50 pc (ring median, 0.13$''$ radius) | 120 pc (aperture) | 250 pc (aperture) |
|---|---|---|---|
| N | 0.99 (38%) | 0.84 (38%) | 0.69 (19%) |
| S | 0.48 (5%) | 0.40 (0%) | 0.37 (0%) |

The 55 pc aperture-mean values are N 0.99 (38%) and S 0.43 (0%). Restricting to literature nuclei within $\pm50\%$ of our scale gives 25%, 50% and 22% for N and 0% for S. The 6.5$''$ blends (0.43 and 0.39) lie below all 21 resolved values and at the 5th percentile of the 20 global ratios. Peak $T_b$ N/S at 50, 55, 120 and 250 pc: HCN $2.03\pm0.09$, $2.02\pm0.09$, $2.02\pm0.05$, $2.50\pm0.06$; HCO$^+$ $0.95\pm0.04$, $0.90\pm0.04$, $0.99\pm0.03$, $1.37\pm0.03$ (the Monte Carlo bias of the peak, 0.01–0.05 K, is not subtracted). The comparison is made at matched physical scale, not matched angular beam at 250 pc (1.29$''$ against 0.3$''$); 8 of the 21 have beams below 110 pc (median ratio 1.08). Nishimura's AGN fractions and CON flags are second-hand to us.

## Appendix C. Provenance

Post ids refer to the project forum; snapshot scripts are in `/snapshots/<id>/scripts/` unless a path is given; `inbox/rNN` are the figure and table scripts supplied to the writing team in earlier rounds (`/workspace/inbox/`), and `new/rNN` those of the latest round (`/workspace/new_evidence/inbox/`).

| result | posts | scripts and outputs |
|---|---|---|
| Fig. 1 overview | p0004, p0005, p0026 | `inbox/r01/fig1_overview.py`, `fig1_numbers.txt` |
| Annulus profile (Fig. 2, Table B2, A.1 array differences) | p0003, p0005, p0008, p0010 | `inbox/r09/r09_profile.py`, `r09_profile.json`; `paper/scripts/fig2_profile.py` |
| 0.30$''$ profile, baseline variants, velocity bins (A.2, Table B2) | p0064, p0065, p0005-r08/r09, p0064-r01 | `team-b-s35`: `s410_q64_imageplane_calib.py`, `s411`, `s415`, `out/q64r2/calib.json`, `final_budget.txt`; `team-a-s33`: `q64r2a_spec.py`, `q64r2c_images.py`, `q64r2d_summary.py`, `q64r2e_null.py`; `new/r01` (`f1_compute.py`), `new/r04` (`f5_baseline_check.py`); velocity bins `team-a-s02/scripts/q1_velratio.py` via `new/r10` |
| Core contrast and budget (Tables B1, B1b) | p0031, p0064, p0065 | `team-b-s35`: `s400_q64_casa_image.py`, `s412_q64_null_and_fft.py`, `s413_q64_contrast_table.py`, `s414_q64_beamcheck_resid.py`, `s415_q64_final_budget.py`, `team/checks/q64r2_check.py`; `team-a-s33`: `q64c_image.py`, `q64i_uvfit3.py`, `q64g_uvannuli.py`; `inbox/r07/r07_budget.py`, `budget_table.md`; `new/q01` (`q01_compute.py`), `new/r03` |
| 12-m versus 12-m+7-m (A.1) | p0064 | `team-a-s33/scripts/q64c_image.py` (sets `main`, `var`); `new/q02` (`q02_table.py`) |
| Estimator spread, flux-complete 1.5$''$/2.2$''$ (Table B1c, Fig. 3) | p0034, p0036 | `team-b-s21`: `s191_q34_casa_image.py`, `s192_q34_measure.py`; `new/r11` (`r11_apertures.py`, `r11_xb4b_imgtaper.py`), `new/r12` (`r12_apertures_10.py`, `r12_D.py`, `r12_tb_err.py`) |
| Contrast versus beam (Fig. 3) | p0031, p0036, p0064, p0066 | `paper/scripts/fig3_contrast_beam.py`; `new/r02` (`f2_plot.py`, `f2_values.csv`) |
| J-dependence 1–0/3–2 (A.6) | p0034, p0036, p0060 | as Table B1c; `new/r12` |
| ACA 3–2 versus 4–3, 4–3 visibility fit, 12-m 4–3 (A.6) | p0053, p0055, p0060, p0066 | `team-a-s34`: `q65a_dump.py`, `q65b_tclean.py`, `q65c_uvfit.py`, `q65d_summary.py`, `q65k_offset.py`, `q65l_fwhmscan.py`, `q65m_tables.py`, `out/q60/imgplane.csv`; `team-a-s30`: `q53g_fixedwin_all.py`, `q53h_mc.py`; `new/r09` (`r09_table.py`; source `s438_q66_32_vs_43.py` not in a snapshot); `new/r08` |
| Earlier 4–3 and J ladder (Table B4) | p0022, p0025, p0027, p0053 | `inbox/r10/r10_ladder.py`, `r10_corr.py`, `ladder.csv`, `ladder_corr.csv` |
| Species (Fig. 4, Tables 1, B5) | p0056, p0058 | `team-b-s30`: `s320_q52_chem_ratios.py`, `s321_q52_ratios.py`, `s322_q52_cont_ref.py`; `team-a-s31`: `q54d_vs_hco.py`; `inbox/r11/r11_table.py`; `paper/scripts/fig3_species.py`; `new/r05` (`f3b_plot.py`) |
| Vibrational HCN limits (Fig. 5, Table B3, A.8) | p0057, p0046–p0048 | `inbox/r06/fig6.py`, `hcn_vib_limits.png` |
| Fits and Bayes factors (Fig. A1, A.7) | p0041, p0042, p0044, p0050, p0054, p0059, p0062, p0063 | `team-b-s33`: `s353_q57_joint_bracket.py`, `s356_q57_physical.py`; `inbox/r04/fig4_r3.py`, `lnBF_table.txt`; `new/r06` (`f4_profile.py`) |
| Isotopologue ratios and thin-line contrast (Table A3) | p0007, p0061 | `inbox/q02/q02_calc.py`, `inbox/q08/q08_calc.py` |
| Temperatures (A.5) | p0012, p0015, p0037, p0040, p0051, p0052, p0063 | `new/r08`; `team-a` snapshots s17–s19, s26, s28 |
| Literature (Fig. 6, Table B6) | p0053, p0055 | `team-b-s28/scripts/s300_q46_lit_context.py`; `inbox/r05/fig5_lit.py`, `percentiles.txt`; `paper/scripts/fig5_literature.py`; `new/r07` (`f3a_plot.py`) |
| Harada et al. N/S at 1.7$''$ | p0058 | `new/r10` (`r10a_harada_ratio.py`, `.txt`); `team-a-s32/scripts/q55f_r10_core.py` |

Two literature checks disagree on whether Falstad et al. (2021) list NGC 3256; no published limit for NGC 3256 is used here.

## References

Aalto, S., et al. 2015, A&A 584, A42 (2015A&A...584A..42A). Butterworth, J., et al. 2022, A&A 667, A131 (2022A&A...667A.131B). Falstad, N., et al. 2021, A&A 649, A105 (2021A&A...649A.105F). Harada, N., et al. 2018, ApJ 855, 49 (2018ApJ...855...49H). Imanishi, M., et al. 2016, AJ 152, 218 (2016AJ....152..218I). Izumi, T., et al. 2016, ApJ 818, 42 (2016ApJ...818...42I). Ledger, B., et al. 2021, MNRAS 504, 5863 (2021MNRAS.504.5863L). Michiyama, T., et al. 2018, ApJ 868, 95 (bibcode not verified). Nishimura, Y., et al. 2024, A&A 686, A48 (2024A&A...686A..48N). Pereira-Santaella, M., et al. 2024, A&A 681, A117 (2024A&A...681A.117P). Tang, X. D., et al. 2019, A&A 629, A6 (2019A&A...629A...6T; bibcode by ADS pattern, not opened).
