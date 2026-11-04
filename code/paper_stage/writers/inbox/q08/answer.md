# q08: direct N/S contrast of (H13CN/H13CO+) at 2.2" (J=1-0)

Script: `q08_calc.py` (this folder). It reads `/workspace/out/q47/ratios10.json`, which `scripts/q47_ratios10.py` writes (Q47b).
Data: 2015.1.00993.S X95d (H13CN, H13CO+) and X961 (HCN, HCO+). This is the flux-complete 2.1x1.7" project. All channels are convolved to a circular 2.2" beam, in r=2.2" apertures.
Estimators: T = 12C-template amplitude fit, A = ±150 km/s integral, P = peak pixel.
Errors: max(formal, empirical from 40 off-source apertures).

## 1. Measured 13C/12C flux ratios (1-0, 2.2") and 12C ratios
Ratios and errors are taken from the json. H13CN/H13CO+ = (H13CN/HCN)·I_HCN / ((H13CO+/HCO+)·I_HCO+).

| est | H13CN/HCN N | S | H13CO+/HCO+ N | S | (H13CN/H13CO+)_N | _S (Jy ratio) | N/S |
|---|---|---|---|---|---|---|---|
| T | 0.0217±0.0079 | 0.0070±0.0068 | 0.0192±0.0039 | 0.0114±0.0092 | 0.95 | 0.35 | 2.7 |
| A | 0.0204±0.0058 | 0.0060±0.0091 | 0.0276±0.0049 | 0.0144±0.0054 | 0.62 | 0.24 | 2.6 |
| P | 0.0319±0.0044 | 0.0216±0.0061 | 0.0306±0.0031 | 0.0327±0.0043 | 0.94 | 0.38 | 2.5 |

The Jy ratio converts to Tb by ×(ν_H13CO+/ν_H13CN)²=1.0095, which cancels in N/S.

## 2. Contrast (H13CN/H13CO+)_N / (H13CN/H13CO+)_S
Errors: Monte Carlo with Gaussian draws of the four ratios, 4e5 draws.

| est | point value | MC median [16,84%] | note |
|---|---|---|---|
| T | 2.7 | 2.4 [0.75, 7.7] | 24% of draws have a non-positive S ratio and are dropped; 95% lower bound 0.13; P(>1)=0.78 |
| A | 2.6 | 1.7 [0.72, 5.4] | 26% of draws dropped; P(>1)=0.73 |
| P | 2.5 | 2.5 [1.8, 3.7] | linear propagation 2.5±0.9; P(>1)=1.00; 2.5% quantile 1.3 |

The 12C HCN/HCO+ contrast in the same beam and apertures is 1.48 (T, A) and 1.57 (P), with N≈0.84-0.90 and S≈0.57-0.58 (Jy ratio).

**Honest reading.** The T and A contrasts are not constrained. S H13CN is consistent with zero (0.007±0.007 and 0.006±0.009), so the 1σ range spans 1 to 5-8. The P estimator gives 2.5 (+1.2/-0.7), but it uses the brightest pixel in each 13C cube. That choice biases the S ratio up, which is why P sits at about 0.022 against about 0.007 for T and A. The P estimate is not independent of the noise.
The value used in the paper should be about 2.5 with a 1σ range of about 0.7-7 (T/A), and it is at most about 1-2σ above 1.

## 3. Implied N/S contrast of HCN/HCO+
If both 13C lines are thin and 12C/13C is the same at N and S, then N(HCN)/N(HCO+) differs between N and S by the contrast in section 2, i.e. about 2.5 (point values 2.5-2.7). This is larger than the observed 12C contrast at 2.2", 1.5, and consistent with the 2.1±0.2 at 50-120 pc. It would imply that the 12C HCN/HCO+ ratio at N is suppressed by opacity (HCN thicker) in the 2.2" beam. The numbers are not significant at the 2σ level for T and A (see above), so this is a hint, not a result.
Two caveats:
- The thin assumption is not demonstrated. The RADEX posteriors for the column+abundance solution predict f=(13C/12C)×(12C/13C)=2.7 at N, i.e. H13CN is not thin there (Q55).
- The ratio is a 2.2" ratio, whereas the core contrast is measured at 0.2-1".

## 4. The items I cannot support from our work
- **12C/13C=80 for N: not ours.** We never adopted 80 for N. The Q47 and Q55 models marginalise over 12C/13C = 32/56/100 (q47_dyn13, q55e_joint). Their posteriors sit at the 100 edge and only constrain X_S/X_N = 0.98 (0.65-1.58, `out/q47/xratio.json`). All our predictions use 100, i.e. the table above is quoted as 13C/12C ratios. The ratios of ratios do not depend on 12C/13C if it is common to N and S. If the paper uses 80 for N, it needs its own citation. I have not found a source in our work, and have not verified one.
- **Δχ² = 5-8 → 1.7-2.1σ: not ours.** I found no such conversion in our scripts or outputs. The conversions I can compute (`q08_calc.py`, two-sided Gaussian-equivalent σ from the χ² tail):
  - 1 dof: Δχ² = 5 → 2.24σ, 8 → 2.83σ.
  - 2 dof: Δχ² = 5 → 1.74σ, 8 → 2.36σ.
  - 3 dof: 5 → 1.37σ, 8 → 2.00σ.
  - 4 dof: 5 → 1.06σ, 8 → 1.69σ.

  The quoted 1.7-2.1σ is not matched exactly by any single dof. 2 dof reproduces the lower end (1.74σ) but gives 2.4σ at the upper end. 3 dof gives 1.4-2.0σ. Our own profile Δχ² statements (Q2, Q42, Q49) are per-parameter, 1-dof, e.g. "T_N profile Δχ² < 5 over 10-150 K". At 1 dof, Δχ² = 5-8 would be 2.2-2.8σ, not 1.7-2.1σ. The authors should state which dof they used, or quote the Δχ² as it is.
