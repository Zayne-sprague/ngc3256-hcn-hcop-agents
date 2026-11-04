# HCN/HCO$^+$ J=3–2 is about twice as high at the starburst nucleus of NGC 3256 as at its low-luminosity AGN at 50–120 pc

**Agent teams A and B (forum), written by the paper team**

## Abstract

An elevated HCN/HCO$^+$ ratio in the submillimetre lines is used to flag accreting black holes in dusty nuclei, but the claim is contested and has rarely been tested where two nuclei of one galaxy can be separated. We use archival ALMA Band 6 data of the merger NGC 3256, which has a starburst-dominated northern nucleus (N) and an obscured, low-luminosity AGN in the southern nucleus (S), at 50–250 pc resolution. The HCN/HCO$^+$ J=3–2 brightness-temperature ratio is 0.99 ± 0.14 at N and 0.44 ± 0.06 at S (aperture means at 50 pc; flux-scale error included), a contrast of 2.1 (+0.1/−0.3) at 50–120 pc; the lower error is a one-sided, still untested allowance for missing short spacings, and the best-controlled, flux-complete contrast is 1.44–1.50 at 2.2$''$ (430 pc), compared with about 1.5 in the 1.7$''$ beam of earlier work. At the peak the contrast is mostly in HCN, with HNC and CN also weak relative to HCO$^+$ at S. No vibrationally excited HCN is seen at either nucleus (beam-filling surface density below 0.02 $L_\odot$ pc$^{-2}$), but the limit on the vibrational-to-ground-state luminosity ratio excludes the compact-obscured-nucleus range only at N. Non-linear radiative-transfer fits to the 2.2$''$ envelope disfavour excitation alone at about 2σ, but abundance, opacity and $^{12}$C/$^{13}$C are not separated, so the cause is not established. In this one galaxy, whose AGN is energetically minor, the HCN/HCO$^+$ ratio does not flag the AGN.

## 1 Introduction

The HCN/HCO$^+$ ratio of rotational transitions is a long-used diagnostic of the state of dense gas in galaxies. Observations of J=3–2 and 4–3 led several groups to report a "submillimetre HCN enhancement": higher HCN/HCO$^+$ in nuclei hosting active galactic nuclei (AGN) than in starbursts, interpreted as the signature of X-ray or mechanical heating, or of AGN-driven chemistry (Izumi et al. 2016; Imanishi et al. 2016). The claim is contested. Resolved ALMA observations of luminous and ultraluminous infrared galaxies find a wide spread of ratios with no trend with AGN luminosity or star-formation rate, and associate the highest ratios with shock-heated outflows (Nishimura et al. 2024). High ratios can also arise from abundance differences (Butterworth et al. 2022) or from infrared pumping of HCN through its bending mode in compact obscured nuclei (CONs; Aalto et al. 2015), in which case vibrationally excited HCN lines should be visible.

NGC 3256 ($D\approx40$ Mpc, $1''\approx194$ pc) has two nuclei separated by about 5$''$ (1 kpc): a starburst-dominated northern nucleus (N) and a southern nucleus that hosts a low-luminosity, obscured AGN (S). Published ratios come from beams of about 1.7$''$ (330 pc), in which N exceeds S by only a factor of about 1.5 and no enhancement at the AGN was seen (Harada et al. 2018); J=1–0 values at about 2$''$ were given by Michiyama et al. (2018). Missing were the comparison at 50 pc, where each nucleus is resolved from its surroundings, the line (HCN or HCO$^+$) in which the nuclei differ, and tests of vibrational pumping. The ALMA archive lists no publication of the Band 6 programme used here.

We ask whether an elevated HCN/HCO$^+$ ratio is a signature of AGN heating in dusty nuclei, using NGC 3256 as a within-galaxy test. Our answer is that the ratio is higher at the starburst nucleus than at the AGN, by about 2 at 50–120 pc and at least 1.5 when flux-complete, that the difference sits in the nitrogen species and is not due to CON-level pumping, and that its cause (abundance, opacity or isotope ratio) is not separated. The ratio therefore does not flag the AGN in this galaxy. This is one galaxy with an energetically minor AGN, so it is a counter-example to a general signature, not a statement about all AGN.

## 2 Data and methods

**Data.** The main data are ALMA Band 6 observations of NGC 3256 (project 2017.1.00379.S) covering HCN and HCO$^+$ J=3–2 (265.9 and 267.6 GHz rest) in three array configurations: a long-baseline 12-m configuration with a native beam of $0.25''$ (about 50 pc), a more compact 12-m configuration, and the 7-m array. Complementary data are HCN, HCO$^+$, H$^{13}$CN, H$^{13}$CO$^+$, HNC and CN J=1–0 (projects 2015.1.00993.S and 2015.1.00412.S), HNC and CN J=3–2 (2015.1.00412.S, 2018.1.00493.S) and HCN and HCO$^+$ J=4–3 (2018.1.00493.S, 2018.1.00223.S). We adopt $v_\mathrm{sys}=2780$ km s$^{-1}$, $D=40$ Mpc and a flux-scale uncertainty of 10% per line (14% per ratio). N and S are measured in the same images, so the flux scale of a given line multiplies both nuclei and cancels in every N/S comparison.

**Measurements.** Primary-beam-corrected cubes are converted to Rayleigh–Jeans brightness temperature $T_b$, convolved with a flux-conserving kernel to a common circular beam, and baseline-subtracted using line-free channels. Velocity-integrated values are taken over $\pm350$ km s$^{-1}$ about the systemic velocity. The ratio $R=T_b(\mathrm{HCN})/T_b(\mathrm{HCO^+})$ is the ratio of integrals in an aperture of radius equal to the beam FWHM centred on the HCO$^+$ peak of each nucleus; the contrast is $C=R_N/R_S$. Unless stated, "0.99", "0.44" and "2.26" are these aperture values at 0.26$''$; the adopted 2.1 is the central value over the estimators (aperture 2.26, single pixel 2.09, peak channel 2.13). Uncertainties are from line-free channels; estimator, window and baseline variations are treated as systematics (Appendix C).

**Missing flux.** Interferometric images lose extended emission. S is more extended and fainter than N and loses more flux, which inflates $C$, so the term is one-sided. Comparing the two 12-m configurations gives $-0.06$ at 1.0$''$. For beams of 1.5–2.2$''$ flux-complete 12-m+7-m images lower $C$ by 0.14–0.30. No flux-complete test exists at $\le0.6''$; it is the main untested systematic of the core contrast.

**Modelling.** One-zone non-LTE (RADEX) models are fitted to the 1–0 and 3–2 lines of HCN, HCO$^+$ and their $^{13}$C isotopologues in a matched 2.2$''$ beam (Appendix A.7).

## 3 Results

### 3.1 The N/S contrast is about 2 at 50–120 pc and at least 1.5 when flux-complete

**HCN/HCO$^+$ J=3–2 is 0.99 at N and 0.44 at S at 50 pc; the contrast is 2.1 (+0.1/−0.3), falls with beam size, and is 1.44–1.50 flux-complete at 2.2$''$.** The ratio map at 0.26$''$ shows a compact high-ratio region centred on N and uniformly low values at S and along the bridge (Fig. 1).

![Figure 1. Peak brightness temperature $T_b$ (RJ, K) of HCO$^+$ J=3–2 (left) and HCN J=3–2 (middle), and the HCN/HCO$^+$ $T_b$ ratio (right), NGC 3256, at a common 0.26$''$ (50 pc) circular beam. Pixels below 3σ are black (left, middle) or grey (right). Crosses mark N, S and the bridge; offsets are in arcsec (east to the left) and pc. The circle is the beam. Ring medians within 0.13$''$ of the markers give ratios of 0.99 (N), 0.48 (S) and 0.51 (bridge); the peak-channel values are maxima and are higher.](figs/fig1_overview.png)

The contrast depends on beam size (Fig. 2): in equal circular beams it is $2.26\pm0.06$ at 0.26$''$ (50 pc), $2.12\pm0.04$ at 120 pc, $1.99\pm0.04$ at 1.0$''$ (194 pc) and $1.64$ at 2.2$''$ (430 pc) with 12-m data only, and $1.44$–$1.50$ at 2.2$''$ with flux-complete 12-m+7-m data (Table C1). The 2.1 adopted at 50–120 pc has a statistical error of about 0.05, an estimator spread of $\pm0.1$ and a one-sided missing-flux term of up to $-0.3$, the value measured at 1.5–2.2$''$ and not yet tested below 1$''$. We therefore quote $C=2.1\,(+0.1/-0.3)$ and never as scale-independent. At 6.5$''$ the nuclei and bridge blend and the ratio is about 0.4 everywhere.

![Figure 2. (a) HCN/HCO$^+$ J=3–2 $T_b$ ratio at N (orange), S (blue) and the bridge (green) versus beam FWHM, in an aperture of radius equal to the beam. Filled circles: long-baseline 12-m data; open squares: compact 12-m configuration at 1.29$''$; triangles at 6.5$''$: compact-array blends. Grey bands: ±14% flux scale, which cancels in N/S. (b) Peak $T_b$ of HCO$^+$ (solid) and HCN (dashed, open) at N and S.](figs/fig2_ratio_vs_beam.png)

**HCN and HCO$^+$ separately.** At the peak the HCO$^+$ brightness is the same at the two nuclei (2.52 versus 2.65 K, N/S 0.95) while the HCN peak is 2.0 times lower at S (1.2 versus 2.5 K). This is estimator dependent. In the ring medians of Fig. 1 HCO$^+$ is 1.67 K at N and 2.33 K at S (N/S 0.72) and HCN 1.54 and 1.05 K (1.47); in the integrated 0–0.25$''$ annulus HCO$^+$ is $236$ and $212$ K km s$^{-1}$ (N/S 1.11); at 250 pc the HCO$^+$ peak N/S is 1.37. The contrast is therefore mostly in HCN at the peak, with a smaller HCO$^+$ share in other estimators (Table 2).

**Table 2.** N/S of HCN and HCO$^+$ at 0.26$''$ by estimator (errors: see Appendix C; the same table at uv-matched 0.6$''$ and 1.0$''$ awaits the q01 images).

| estimator | HCO$^+$ N/S | HCN N/S | HCN/HCO$^+$ contrast |
|---|---|---|---|
| peak $T_b$ | 0.95 | 2.0 | 2.1 |
| ring median (0.13$''$) | 0.72 | 1.47 | 2.06 |
| integrated, 0–0.25$''$ annulus | 1.11 | – | – |

**The N ratio rises inward.** In annuli around each nucleus at 0.26$''$, the N ratio is $1.00\pm0.02$ inside 0.25$''$ (49 pc), $0.88\pm0.02$ at 0.25–0.5$''$, $0.65\pm0.01$ at 0.5–1$''$, $0.43\pm0.02$ at 1–2$''$ and $0.38\pm0.05$ at 2–3$''$, whereas S is $0.44\pm0.01$, $0.39\pm0.02$, $0.33\pm0.02$ and $0.34\pm0.03$ (Fig. 3, Table C2). The N ratio falls steadily from the centre over the inner 200 pc and reaches the S value only beyond 1$''$ (190 pc); S is nearly flat. Between 55 and 250 pc the N ratio falls from 0.99 to 0.69 (index $-0.23\pm0.03$) and the S ratio from 0.43 to 0.37 (index $-0.10\pm0.03$). We do not fit a core size; structure inside 50 pc is unresolved.

![Figure 3. HCN/HCO$^+$ J=3–2 $T_b$ ratio (left) and annulus-mean integrated $T_b$ of HCN (solid) and HCO$^+$ (dotted) at N (middle) and S (right) as a function of radius from the HCO$^+$ peak. Blue: long-baseline data at 0.26$''$; orange: convolved to 1.0$''$; green: compact 12-m configuration at 1.0$''$ (flux cross-check). Error bars are statistical. Grey vertical lines mark the beam FWHM.](figs/fig3_ratio_profile.png)

Further checks are in Appendix A.1–A.3: the two 12-m configurations agree within 0.03 inside 1$''$, and the contrast is the same in three observing blocks to within 9%. **What this means:** at the resolution where the nuclei separate the AGN nucleus has the lower ratio by about a factor of two, the difference is mostly in HCN at the peak, and it shrinks to at least 1.5 in the flux-complete 2.2$''$ beam.

### 3.2 HCN, HNC and CN are all weak relative to HCO$^+$ at S

**In flux-complete J=1–0 data at 2.5$''$ the N/S contrast of the ratio to HCO$^+$ is $1.56\pm0.03$ for HCN, $1.90\pm0.11$ for HNC and $2.10\pm0.04$ for CN (Fig. 4b, Table 1), so the contrast belongs to the nitrogen species, not to HCN alone.** The J=3–2 values from 12-m data alone (1.4$''$: $2.14\pm0.05$, $1.84\pm0.14$, $2.56\pm0.09$; 2.5$''$: $1.93$, $2.26$, $2.99$) are upper values of the contrast, since the HCN/HCO$^+$ contrast at 2.2$''$ drops from 1.65 to 1.44–1.50 with 7-m data and no 7-m data exist for HNC and CN. They agree in sign and are in Appendix A.10; they are not compared with the 1–0 column. The N/S contrast of HNC/HCN is $1.22\pm0.07$ at 1–0 and of CN/HCN $1.35\pm0.03$, so HNC and CN are not identical to HCN, but are at least as deficient at S as HCN against HCO$^+$ in the flux-complete data. The errors are statistical only.

![Figure 4. Flux-complete N/S contrasts only. (a) HCN/HCO$^+$ contrast at J=1–0 (blue) and J=3–2 (orange) for matched 12-m+7-m circular beams of 2.2$''$, 1.5$''$ and 2.5$''$; this probes the 300–500 pc blend, not the core. (b) Contrast of HCN, HNC and CN relative to HCO$^+$ at J=1–0, 2.5$''$. Error bars are statistical.](figs/fig4_contrast_J_species.png)

**Table 1.** $T_b$ ratios (RJ, $\pm350$ km s$^{-1}$), flux-complete J=1–0 at 2.5$''$, one pipeline. Errors statistical.

| ratio | N | S | N/S |
|---|---|---|---|
| HCN/HCO$^+$ | 0.912 | 0.585 | 1.56 ± 0.03 |
| HNC/HCO$^+$ | 0.474 | 0.249 | 1.90 ± 0.11 |
| CN/HCO$^+$ | 1.590 | 0.755 | 2.10 ± 0.04 |
| HNC/HCN | 0.520 | 0.427 | 1.22 ± 0.07 |
| CN/HCN | 1.743 | 1.292 | 1.35 ± 0.03 |

The same contrast at 1–0 and 3–2 in flux-complete 1.5–2.5$''$ beams (Fig. 4a; $D=C_{3-2}/C_{1-0}=1.04\pm0.07$, where the $\pm0.07$ spans two estimators) is a statement about the 300–500 pc blend, where the contrast is about 1.5. The core is untested: no 1–0 observation at $\le1''$ is flux-complete and the 4–3 contrast at 1.0$''$ is undecided between 2.3 and 3.7 (Appendix A.6). CN is not enhanced at S relative to HCN, which is not the signature of an X-ray-dominated region, although CN is also raised by UV fields and shocks and so this does not exclude one. **What this means:** the difference between the nuclei is a nitrogen-species-to-HCO$^+$ contrast, not an effect specific to HCN.

### 3.3 No CON-level vibrational HCN emission at either nucleus

**Vibrationally excited HCN $v_2=1f$ J=3–2 (267.199 GHz) is not detected; the beam-filling surface density is $<0.02\,L_\odot$ pc$^{-2}$ at both nuclei, but the luminosity-ratio limit excludes the CON range only at N.** The 3σ limits in the 48 pc beam are $\alpha=I(v_2{=}1)/I(v{=}0)<0.05$ at N and $<0.10$ at S (Fig. 5), i.e. $L_\mathrm{vib}<29\,L_\odot$ and $<25\,L_\odot$ and $\Sigma_\mathrm{vib}<0.017$ and $<0.015\,L_\odot$ pc$^{-2}$, 59–67 times below the CON criterion of $1\,L_\odot$ pc$^{-2}$ (Falstad et al. 2021) for a source filling the beam. The vibrational-to-ground-state luminosity ratio is $<0.025$–0.05 at N, below the CON range 0.06–0.18, and $<0.085$–0.10 at S, which overlaps it. The limit therefore excludes CON-level pumping tentatively at N only. It assumes that the vibrational profile equals the HCN $v=0$ profile and depends on the model of the HCO$^+$ red wing (the widest of four models is used; Table C3). A buried source smaller than 47 pc diluted in the beam, or moderate pumping ($T_\mathrm{vib}=150$–300 K), is not excluded. This matters at S, which JWST observations show to be dust-embedded with an obscured 4.7 $\mu$m continuum (Pereira-Santaella et al. 2024), where a compact pumping source diluted in our beam cannot be ruled out. There is no shock signature in the nuclei (CS 7–6/HCN 4–3 uniform, SO/HCO$^+<0.014$; beams and errors pending).

![Figure 5. Spectra of the HCN $v_2=1f$ J=3–2 line (rest 267.199 GHz, green dotted) at N (left) and S (right) in the 0.24$''$ (48 pc) beam. Grey: data (line-free rms 0.42 mJy beam$^{-1}$ per 3.9 MHz channel); blue: model of the HCO$^+$ red wing; orange: wing model plus the HCN $v=0$ profile scaled to the 3σ limit (0.033 and 0.090 for this wing model; 0.05 and 0.10 for the widest of four models).](figs/fig5_hcn_vib.png)

**What this means:** neither nucleus has CON-level vibrational surface brightness, so pumping of a beam-filling region does not explain the contrast; a compact buried source at S is not excluded.

### 3.4 The cause of the contrast is not separated

**Excitation alone is disfavoured at about 2σ in the 2.2$''$ envelope, but the data do not choose between an abundance difference, opacity and a different $^{12}$C/$^{13}$C ratio, and no model reproduces the observed contrast.** The fit uses the 1–0 and 3–2 lines of HCN, HCO$^+$, H$^{13}$CN and H$^{13}$CO$^+$ in a flux-complete 2.2$''$ beam, where the contrast is $1.45\pm0.01$ (1–0) and $1.50\pm0.02$ (3–2); it applies to the envelope, not to the 2.1 core contrast. With equal abundance ratios, excitation-only and column-density-only models reach contrasts of 1.1–1.3 and are worse than the full model by $\Delta\chi^2=5$–8 (1.7–2.1σ). Even the best abundance-only model reaches only 1.33 (full model 1.29) against 1.45 observed (Fig. 6, Appendix A.7); the Bayes-factor preference for it over excitation-only runs from 4 to 800 depending on the assumed errors and is not used as a measurement.

The isotopologue data do not exclude an abundance-only difference. In the 2.2$''$ 1–0 data H$^{13}$CN/HCN is $0.022\pm0.008$ at N (and consistent with zero, $<0.03$ at 3σ, at S) and H$^{13}$CO$^+$/HCO$^+$ is $0.028\pm0.005$ at N and $0.014\pm0.005$ at S (aperture estimator). The double ratio $(^{13}/^{12})_N/(^{13}/^{12})_S$, which is 1 for an abundance-only difference, is 1.5–3.4 for HCN and 0.9–1.9 for HCO$^+$ depending on estimator, with $P(\le1)=0.04$–0.65 and the 95% lower bound below 1 for every estimator (Table A3). A posterior for $X_S/X_N=(^{12}\mathrm{C}/^{13}\mathrm{C})_S/(^{12}\mathrm{C}/^{13}\mathrm{C})_N$ has median 0.98 and 95% interval 0.45–2.55. There is thus no evidence that the isotope ratios differ, and also no evidence against it; the data lean weakly toward higher $^{13}$C/$^{12}$C at N for HCN. If the $^{13}$C lines are thin and the $^{12}$C lines thin, $^{12}$C/$^{13}$C at N is 45 (+27/−14) from HCN and 36 (+8/−6) from HCO$^+$, below the value of about 80 adopted for N; if the $^{12}$C lines are thick these are lower limits and imply $\tau_{12}\approx1$–3 (HCN) and 1–4 (HCO$^+$) for $^{12}$C/$^{13}$C of 60–135. The RADEX posterior at N gives $\tau_{10}(\mathrm{HCN})=4.6\,[3.0,8.7]$ and $\tau_{10}(\mathrm{HCO^+})=0.5\,[0.2,4.7]$, and $\tau_{10}(\mathrm{HCN})\approx5$–12 at both nuclei. Opacity of the HCN lines is therefore likely to matter, and the "abundance excess" reading is one of three that the data allow. A single excitation state is rejected ($p<0.04$): the extended gas has $T_k\gtrsim80$ K at $n\sim10^4$ cm$^{-3}$ and a compact HCN-rich component of about 0.5$''$ is required, whose temperature is undetermined; S is cool ($T_k=17\pm1$ K). Continuum and dust background limit the envelope background temperature to $\lesssim3.2$ K, disfavouring radiative pumping.

![Figure 6. $\Delta\chi^2$ relative to the full model as a function of the model 1–0 HCN/HCO$^+$ N/S contrast when the abundance ratio is forced equal at N and S (excitation-only blue, column-density-only orange; solid: $^{12}$C/$^{13}$C 40–100, dashed: 20–200). The observed contrast is shown in magenta. Squares/diamonds: abundance-only and full-model best fits. Right: natural logarithm of the Bayes factor of abundance-only over excitation-only for error models and priors.](figs/fig6_hypotheses.png)

**What this means:** the envelope contrast is not easily produced by a difference in rotational excitation alone, but abundance, HCN opacity and isotope ratio remain open, and the result concerns the 2.2$''$ beam.

### 3.5 Placement in the distribution of resolved ratios

**At 50 pc S is tied with the lowest resolved HCN/HCO$^+$ J=3–2 ratio of any luminous-infrared-galaxy nucleus, while N is typical.** Nishimura et al. (2024) list 21 resolved J=3–2 ratios of U/LIRG nuclei from 0.45 to 2.67 (median 1.13; not corrected for missing flux; their definition of the ratio, brightness or integrated flux and its velocity window, is being checked against ours). Our N ratio of 0.99 at 55 pc lies at the 38th percentile, and 0.84 and 0.69 at 120 and 250 pc at the 38th and 19th (Fig. 7a; percentiles carry no uncertainty and rest on 4–10 nuclei of comparable scale). S, with 0.43, 0.40 and 0.37, is at the minimum, but the lowest literature value (0.445, ESO 148-IG002) is equal within errors, so S is tied with it. The sample contains one AGN-dominated nucleus (IRAS F05189-2524, 1.63), so it cannot test whether the AGN order is reversed. At 1–0 our ratios (0.84–0.94 at N, 0.57–0.58 at S) agree with earlier values of 0.96 and 0.645. At J=4–3 the beams are about 4$''$ (800 pc), the nuclei blend, and N (0.52) and S (0.29) both fall within the starburst range of Izumi et al. (2016; 0.24–1.03, quoted from the paper and being re-checked); at 188 pc the values are 0.81 and 0.22. S has the stronger AGN and free-free continuum yet the lower ratio, which disfavours AGN- or XDR-driven enhancement of HCN there.

![Figure 7. HCN/HCO$^+$ in context. (a) J=3–2 $T_b$ ratio versus spatial scale: circles are the 21 resolved nuclei of Nishimura et al. (2024) coloured by AGN bolometric fraction (grey <0.10, blue 0.10–0.25, pink diamond ≥0.25; open: not listed); stars are N (orange) and S (blue) at 55, 120 and 250 pc; open stars are 6.5$''$ blends affected by flux loss. (b) J=4–3 ratios of N, S and the combined system (800 pc beam) on the AGN and starburst bands of Izumi et al. (2016).](figs/fig7_literature.png)

**What this means:** the ratio at the AGN is not high; this is a counter-example to a general AGN signature from one galaxy with a low-luminosity AGN, not a statement about all AGN.

## 4 Discussion

**Significance.** Published ratios for NGC 3256 are blends in a 1.7$''$ beam in which N and S differ by 1.5 and no enhancement at the AGN was noted (Harada et al. 2018, a different programme from ours). We add the within-galaxy comparison at 50 pc. A reader should now believe: (i) the starburst nucleus has HCN/HCO$^+$ about twice that of the AGN at 50–120 pc (2.1, +0.1/−0.3) and at least 1.44 flux-complete at 430 pc; (ii) the difference involves HCN, HNC and CN against HCO$^+$ and is mostly in HCN at the peak; (iii) the N ratio rises inward over the inner 200 pc; (iv) no CON-level vibrational surface brightness is present. The consequence for the diagnostic is that its sign cannot be assumed: the AGN nucleus has the lowest ratio in this galaxy. This agrees with the lack of trend found by Nishimura et al. (2024) and offers a within-galaxy example of why such trends may be absent. Ledger et al. (2021) discussed CN/HCN and HCN 3–2/1–0 differences between the two nuclei; their numbers are not yet compared with ours.

**Alternatives and their status.** *Continuum or dust*: the species ratios compare lines at similar frequencies, so dust cancels; no ratio-versus-continuum test within a nucleus is yet available (A.3). *Beam mixing and missing flux*: not excluded; the one-sided term is up to −0.3 below 1$''$. *Different HCO$^+$*: not at the peak (N/S 0.95), but 0.72–1.11 in other estimators. *Vibrational pumping*: no CON-level surface brightness (<0.02 $L_\odot$ pc$^{-2}$); luminosity-ratio limit excludes the CON range only at N, and a buried compact source is open. *XDR chemistry*: CN is not enhanced at S, which argues against but does not exclude it. *Shocks*: no CS or SO signature at the nuclei (numbers pending). *Excitation alone*: disfavoured at about 2σ at 2.2$''$ and about 7% of equal-abundance excitation models reach a 1–0 contrast of 1.35. *Abundance, opacity, $^{12}$C/$^{13}$C*: not separated.

**Limits and what would settle them.** (1) The 2.1 core contrast is measured at J=3–2 with 12-m data at $\le0.6''$ and its missing-flux term is untested (up to −0.3); uv-matched 12-m+7-m images tapered to 0.6$''$ and 1.0$''$ (request q01) would settle it. (2) J-independence is established at $\ge1.5''$ only; a 4–3 pair at $\le1''$ with short spacings would test the core. (3) Pumping limits exclude CON-level surface brightness for a beam-filling source only. (4) The cause is not established; an H$^{13}$CN and H$^{13}$CO$^+$ 1–0 measurement about 10 times deeper (uncertainty about 0.003) at about 2$''$, or HN$^{13}$C/HNC and $^{13}$CO/C$^{18}$O in matched beams, would separate opacity from abundance. (5) The temperature of the N core is unconstrained; HCN 4–3/3–2 and thermometer lines at $\le0.5''$ would give it. (6) One galaxy with an energetically minor AGN: calibration needs resolved ratios and $^{13}$C species in a sample of AGN of known luminosity. Distance, flux scale and image-plane (not uv-plane) beam matching are further assumptions. Dense gas near S (Appendix B) with ratios of 1–1.4 is a single-array, 3–5σ result.

## 5 Conclusions

In NGC 3256 the HCN/HCO$^+$ J=3–2 brightness-temperature ratio is 0.99 at the starburst nucleus and 0.44 at the low-luminosity AGN at 50 pc, a contrast of 2.1 (+0.1/−0.3) that falls with beam size to 1.44–1.50 in a flux-complete 2.2$''$ beam, with an untested missing-flux term at $\le0.6''$. The difference is mostly in HCN at the peak, with HNC and CN also weak relative to HCO$^+$ at S. No CON-level vibrational HCN surface brightness is seen at either nucleus. Excitation alone is mildly disfavoured at 2.2$''$, but abundance, opacity and $^{12}$C/$^{13}$C are not separated. In this one galaxy the HCN/HCO$^+$ ratio does not flag the AGN; a calibration needs resolved multi-line and isotopologue data in a sample of AGN.

---

## Appendix A. Tests supporting the main results

**A.1 Flux recovery and array consistency (Section 3.1).** The Band 6 ratio agrees within 9% across three observing blocks of the long-baseline and compact 12-m configurations; flux-recovery and ladder-ratio systematics are $\lesssim10\%$. At 0.5–1.5$''$ the long-baseline configuration loses 15–20% of the flux at S (HCO$^+$ $-16\%$, HCN $-21\%$) but only about 4% of the HCN flux at N; inside 0.5$''$ the arrays agree to 2–8%. The result does not depend on the aperture, window or baseline within the budget of Appendix C.

**A.2 Line wings (Section 3.1).** An apparent rise of HCN/HCO$^+$ in the line wings is beam mixing of the N core with surrounding gas, and the deficit at S is present at all velocities; earlier wing-ratio candidates are not retained.

**A.3 Continuum, velocity and observing blocks (Section 3.1).** At equal HCO$^+$ the ratio is 0.63 at N and 0.43 at S; this shows that the nuclei differ at equal HCO$^+$ but not that the ratio is independent of continuum brightness within a nucleus, which is untested (requested). The ratio is flat in velocity within each nucleus (0.95–1.09 at N, 0.36–0.43 at S; window, beam and errors pending) and the three observing blocks agree within 9% (individual values pending). The H40$\alpha$ comparison is marginal (1.5–3.3σ) and not used.

**A.4 Saturation (Section 3.4).** The limit HCN/H$^{13}$CN $>8.3$ at N is far below any plausible $^{12}$C/$^{13}$C and says nothing about saturation; opacity is discussed in Section 3.4.

**A.5 Excitation, temperature and density (Section 3.4).** Ladder fits for N and S, including degeneracies between $T_k$ and $n$ and beam filling, are consistent with the one-zone result. Thermometer lines (CO 2–1, $^{13}$CO cores) give $T_k(\mathrm{S})\ge26$–40 K. A nested-sampling fit gives $T_c(\mathrm{S})=17\pm1$ K and an unconstrained $T_c(\mathrm{N})$. Two-zone models and a cold optically thick branch are disfavoured but depend on priors; H$^{13}$CN 4–3 is uninformative except at $n\ge10^6$ cm$^{-3}$.

**A.6 J-dependence: 1–0 versus 3–2 and 4–3 (Section 3.2).** At matched flux-complete beams the contrast is $1.46\pm0.02$ and $1.48\pm0.04$ (1–0, 3–2) at 2.2$''$, $1.62\pm0.06$ and $1.60\pm0.07$ at 1.5$''$; $D=1.00\pm0.05$ (12-m+7-m) and $1.07\pm0.05$ (a different estimator), adopted $1.04\pm0.07$ with the uncertainty spanning the two estimators. The 3–2/1–0 HCN ratio at S/N is $0.93\pm0.05$. This does not test the core. **4–3.** Matching HCN 4–3 and HCO$^+$ 4–3 (two programmes) to 3–2 at one beam, the 4–3/3–2 $T_b$ ratios at 1.0$''$ (raw) are HCN $0.646\pm0.011$ (N), $0.471\pm0.034$ (S); HCO$^+$ $0.663\pm0.013$ (N), $0.838\pm0.020$ (S). Two bounding corrections for missing 4–3 flux (a uniform pedestal; and a fractional correction as if all lost flux sat in the core) give contrasts of 3.34 and 2.26 at 1.0$''$ and 3.26 and 2.37 at 1.5$''$ (raw: 3.67, 3.85), against 3–2 contrasts of 2.12 and 2.03. The ratio of 4–3 and 3–2 contrasts is 1.73, 1.58 and 1.07 (raw, pedestal, fractional) at 1.0$''$. The recoveries rely on a single-tuning 4.5$''$ observation with an unexplained 1.5$\times$ absolute flux offset; the 1.6 blended contrast and the uv-plane two-component fits (1.4–1.9) are not core measurements. The compact-to-extended behaviour of S HCN (4–3 flux in the high-resolution data is 0.61 of the single tuning) is the crux. The continuum edge ramp biases the HCN 4–3 measurement by 0–10%.

**A.7 Isotopologues and the joint fit (Section 3.4).** Bayes factors (abundance-only over excitation-only, natural log), always positive but strongly error-model dependent:

| prior | errors halved | default | errors doubled |
|---|---|---|---|
| $^{12}$C/$^{13}$C 40–100 | 6.70 (BF 809) | 3.71 (41) | 1.92 (6.8) |
| 40–100, physical | 6.68 (794) | 3.56 (35) | 1.78 (6.0) |
| 20–200 | 6.40 (603) | 3.27 (26) | 1.43 (4.2) |
| 20–200, physical | 6.46 (642) | 3.08 (22) | 1.37 (3.9) |

The "physical" priors are $\log n\ge3.8$, size $\le100$ pc and $X(\mathrm{HCO^+})\le10^{-8}$. Best-fit 1–0 contrasts are excitation-only 1.09, column-only 1.13, abundance-only 1.33 and full 1.29 (observed 1.45); the 3–2 contrasts are 1.23, 1.35, 1.57, 1.53 (observed 1.50). $N(\mathrm{HCN})/N(\mathrm{HCO^+})=3.2\,[2.5,4.0]$ at N and $2.0\,[1.6,2.5]$ at S; moving the ratio at S from 0.4 to 1 needs a factor of 3–4 in abundance ratio.

**Table A3.** $^{13}$C/$^{12}$C ratios at N and S (1–0, 2.2$''$ aperture) and double ratio (13/12)$_N$/(13/12)$_S$ by estimator (T template, A aperture, P peak pixel; 1σ errors; MC with Gaussian errors).

| line | est. | N | S | double ratio, nominal | MC median [16–84%] | P(≤1) |
|---|---|---|---|---|---|---|
| H$^{13}$CN/HCN | T | 0.0217±0.0079 | 0.0070±0.0068 | 3.1 | 2.6 [1.3–7.0] | 0.08 |
| H$^{13}$CN/HCN | A | 0.0204±0.0058 | 0.0060±0.0091 | 3.4 | 2.3 [1.1–6.6] | 0.09 |
| H$^{13}$CN/HCN | P | 0.0319±0.0044 | 0.0216±0.0061 | 1.5 | 1.5 [1.1–2.1] | 0.08 |
| H$^{13}$CO$^+$/HCO$^+$ | T | 0.0192±0.0039 | 0.0114±0.0092 | 1.7 | 1.5 [0.9–3.7] | 0.22 |
| H$^{13}$CO$^+$/HCO$^+$ | A | 0.0276±0.0049 | 0.0144±0.0054 | 1.9 | 1.9 [1.3–3.1] | 0.04 |
| H$^{13}$CO$^+$/HCO$^+$ | P | 0.0306±0.0031 | 0.0327±0.0043 | 0.9 | 0.9 [0.8–1.1] | 0.65 |

The MC ratio is ill defined where S is consistent with 0 (S≤0 in 11–26% of draws). Opacity at N for thin $^{13}$C lines (single zone, same excitation): $^{12}$C/$^{13}$C = 40, 60, 90, 135 gives $\tau_{12}(\mathrm{HCN})=0$, 0.6, 1.6, 2.8 and $\tau_{12}(\mathrm{HCO^+})=0.23$, 1.2, 2.3, 3.7. H$^{13}$CN 4–3 gives HCN/H$^{13}$CN $>9$ at N and $>6.7$ at S. An earlier draft quoted double ratios of 5.8 and 3.1 from a 384-variant budget; team A could not reproduce them and they are withdrawn.

**A.8 Pumping and background (Sections 3.3, 3.4).** Continuum and dust background limits ($T_\mathrm{bg}\lesssim3.2$ K) disfavour radiative pumping of the envelope. Limits at coarser beams are in Table C3.

**A.10 12-m-only J=3–2 species contrasts (Section 3.2).** HCN/HCO$^+$, HNC/HCO$^+$, CN/HCO$^+$ at 3–2: 1.4$''$ $2.14\pm0.05$, $1.84\pm0.14$, $2.56\pm0.09$ (N 0.833, 0.477, 0.893; S 0.389, 0.259, 0.348); 2.5$''$ $1.93\pm0.05$, $2.26\pm0.21$, $2.99\pm0.14$. HNC/HCN is $0.86\pm0.06$ (1.4$''$) and $1.17\pm0.10$ (2.5$''$), CN/HCN 1.2–1.55. These lack 7-m data and are upper values of the contrast; the HNC value at 1.4$''$ is below HCN. Normalised to the 1.3 mm continuum, HCN, HNC and CN are 2.4–3.1 times and HCO$^+$ 1.3–1.4 times weaker at S (Table C5), a factor that scales with the unmeasured $T_{d,N}/T_{d,S}$ (0.5–2) and is not used.

## Appendix B. Side findings

**B.1 Dense gas in the S outflow.** A blue lobe and a red ridge near S show HCN/HCO$^+\approx1$–1.4 in a single array at 3–5σ, with a dense-gas mass of about $10^7\,M_\odot$ and a rate of 15–60 $M_\odot$ yr$^{-1}$. This is not established in a second array and is absent at J=1–0; we note it as possible shocked, HCN-rich outflow gas, consistent with Nishimura et al. (2024).

**B.2 Methanol.** CH$_3$OH is detected, extended and cold; an enhancement in the bridge was not reproduced.

**B.3 Disk and clump gas.** The HCN/HCO$^+$ ratio in disk and clump gas declines with J (about 11% at 2–1 and 6–10% at 1–0; 3–2 not measurable).

## Appendix C. Measurement conventions and tables

**Conventions.** $T_b$ is the Rayleigh–Jeans brightness temperature with the $\nu^2$ factor applied. The velocity window is $\pm350$ km s$^{-1}$ about 2780 km s$^{-1}$ (N 2785, S 2750 in some tables), except in the fits of Section 3.3 ($\pm150$ km s$^{-1}$). The baseline is a constant (per-pixel median of line-free channels at $<-450$ and $>+600$ km s$^{-1}$); a linear baseline was used for the 4–3, CN and species tables. Apertures have radius equal to the beam FWHM centred on the HCO$^+$ 3–2 peak (N: 156.96346°, −43.90388°; S: 156.96339°, −43.90533°). Noise is from line-free channels. The 0.26$''$ beam is the common circular convolution of a native 0.246$''\times0.232''$ beam. 1$''$ = 194 pc ($D=40$ Mpc).

**Table C1. N/S contrast of HCN/HCO$^+$ J=3–2 by beam.**

| beam | aperture (long-baseline) | stat | single pixel | peak channel | compact 12-m | 12-m+7-m ($\pm350/\pm150$) | missing-flux term |
|---|---|---|---|---|---|---|---|
| 0.26$''$ (50 pc) | 2.26 | 0.06 | 2.09 | 2.13 | – | – | – (untested) |
| 0.28$''$ (55 pc) | 2.32 | 0.06 | 2.14 | 2.24 | – | – | – (untested) |
| 0.62$''$ (120 pc) | 2.12 | 0.04 | 2.24 | 2.04 | – | – | – (untested) |
| 1.0$''$ (194 pc) | 1.99 | 0.04 | 2.16 | 1.84 | 1.94 | – | −0.06 |
| 1.29$''$ (250 pc) | 1.90 | 0.04 | 2.10 | 1.83 | 1.84 | – | −0.06 |
| 1.5$''$ (291 pc) | 1.83 | 0.04 | 2.04 | 1.83 | 1.81 | 1.53 / 1.63 | −0.3 / −0.2 |
| 2.2$''$ (427 pc) | 1.64 | 0.03 | 1.88 | 2.00 | 1.65 | 1.44 / 1.50 | −0.20 / −0.14 |

Adopted: $2.1\,(+0.1/-0.3)$ (0.26–0.62$''$; statistical 0.05 and estimators $\pm0.1$ in quadrature give $+0.1$; the one-sided missing-flux term, untested below 1$''$, is taken as up to $-0.3$ from the 1.5–2.2$''$ result), $2.0\,(+0.1/-0.2)$ (1.0$''$; $-0.06$ measured between the two 12-m configurations, the 12-m+7-m effect at 1.0$''$ is unmeasured), $1.9\pm0.15$ (1.29$''$), $1.5\pm0.1$ (2.2$''$, 12-m+7-m, $\pm150$ km s$^{-1}$; 1.44 at $\pm350$). The flux scale (14%) cancels in $C$. Changing the window from $\pm350$ to $\pm150$ km s$^{-1}$ changes the flux-complete 2.2$''$ contrast by 0.06 and the 1.5$''$ one by 0.10; a blue-only versus two-sided baseline changes $R$ by $\le0.01$. Individual ratios at 0.26$''$ are $R_N=0.99\pm0.14$, $R_S=0.44\pm0.06$ including the flux scale. The 1–0 contrast also depends on aperture (1.3–1.5).

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

Figures 1–3, 5–7 are the teams' figures; Figure 4 was drawn by the paper team from the flux-complete tabulated values (`paper/scripts/fig4_family.py`).

## References

Aalto, S., et al. 2015, A&A 584, A42. Butterworth, J., et al. 2022, A&A 667, A131. Falstad, N., et al. 2021, A&A (CON vibrational HCN; volume and bibcode pending). Ledger, B., et al. 2021, MNRAS (details pending). Pereira-Santaella, M., et al. 2024 (JWST, NGC 3256; details pending). Harada, N., et al. 2018, ApJ 855, 49. Imanishi, M., et al. 2016, AJ 152, 218. Izumi, T., et al. 2016, ApJ 818, 42. Michiyama, T., et al. 2018, ApJ 868, 95. Nishimura, Y., et al. 2024, A&A 686, A48.
