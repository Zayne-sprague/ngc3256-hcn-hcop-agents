# Q52 work summary: verification of team-a p0056 (chemistry at N and S)

**Verdict: CONFIRMED, with a refinement.**

## What I did
- I reproduced the HNC, CN, HCN and HCO+ integrals at N and S independently of team-a's code. The inputs are archive cubes convolved to circular beams, with a pixel read at fixed N and S positions (regions.py), a linear baseline from line-free channels, and a ±350 km/s window about v_sys (N 2785, S 2750 km/s). The ±250 km/s window is the window systematic.
  - J=1-0: 2015.1.00993 X961 gives HCN, HCO+ and HNC. 2015.1.00412 X145/X149 is the cross-check. 2018.1.00223 X4d spw27 gives CN, summed over both hyperfine groups.
  - J=3-2: Xb4d gives HCN and HCO+. 2015.1.00412 X169 gives HNC (0.75"). 2018.1.00493 X1df1 gives CN, summed over the whole complex.
- Beams are 2.5" and 3.5" at 1-0, and 1.4" and 2.5" at 3-2.
- Scripts and outputs:
  - `scripts/s320_q52_chem_ratios.py` writes `out/q52/s320.{json,txt}` (the integrals).
  - `scripts/s321_q52_ratios.py` writes `out/q52/s321.{json,txt}` (ratios and N/S double ratios).
  - `scripts/s322_q52_cont_ref.py` writes `out/q52/s322.txt` (N/S of the lines and of the 1.1 mm continuum).
  - `scripts/s323_update_report_q52.py` adds the 5 Q52 entries to `/output/report.json`. It validates, and the report now has 238 results.

## Results (stat errors only)
| Ratio | N | S | N/S |
|---|---|---|---|
| **J=1-0, 2.5"** | | | |
| HCN/HCO+ | 0.912±0.008 | 0.585±0.008 | 1.56±0.03 |
| HNC/HCO+ | 0.474±0.013 | 0.249±0.010 | 1.90±0.11 |
| CN/HCO+ | 1.59±0.01 | 0.755±0.010 | 2.10±0.04 |
| HNC/HCN | 0.520±0.014 | 0.427±0.018 | 1.22±0.07 |
| CN/HCN | 1.74±0.01 | 1.29±0.02 | 1.35±0.03 |
| **J=3-2, 1.4"** | | | |
| HCN/HCO+ | 0.833 | 0.389 | 2.14±0.05 |
| HNC/HCO+ | 0.477 | 0.259 | 1.84±0.14 |
| CN/HCO+ | 0.893 | 0.348 | 2.56±0.09 |
| HNC/HCN | 0.572 | 0.666 | 0.86±0.06 |
| CN/HCN | 1.07 | 0.89 | 1.20±0.04 |

- The 3.5" and 2.5" beams and the 00412 cross-check give the same N/S within 5%. At 3-2 and 2.5", HNC/HCN is 0.54 at N and 0.46 at S (N/S 1.17±0.10).
- p0056 is reproduced to 0-7% (p0056: 1-0 CN/HCN 1.76/1.39, HNC/HCN 0.52/0.42; N/S of X/HCO+ 1.58/1.96/2.00 at 1-0 and 2.0/1.95/2.45 at 3-2). My S CN/HCN is 7% lower than theirs.
- Literature:
  - Harada+2018 Table 13 (1.7", peak pixel, `lit/harada_table13_parsed.csv`) shows the same sense of N/S:
    - HCN/HCO+ N/S 1.48
    - HNC/HCO+ N/S 1.69
    - CN/HCO+ N/S 1.65
    - HNC/HCN N/S 1.12
    - CN/HCN N/S 1.10
  - Ledger+2021 CN/HCN ~1.7 at N is reproduced. I used Ledger only as quoted in p0056 and did not re-read the paper.
- Which species is odd:
  - N/S of the line integrals, 1-0 at 2.5": HCO+ 1.42, HCN 2.21, HNC 2.69, CN 2.98. At 3-2 and 1.4": 1.45, 3.11, 2.67, 3.73.
  - The 1.1 mm continuum N/S is only 1.25 (1.4") and 1.54 (2.5").
  - HCO+(3-2)/continuum is equal at N and S (1540 vs 1330 K km/s per K at 1.4").
  - HCN/continuum is 2.5x lower at S, HNC/continuum 2.1x and CN/continuum 3.0x, at 1.4". At 2.5" they are 1.9x, 2.2x and 2.9x.
  - So the N/S contrast is abundance-wide within the N-bearing group (HCN, HNC, CN all low at S). It is not HCN-specific.
  - It is a deficit of that group at S (or an enhancement at N), not an HCO+ excess.
- No XDR signature: at S, CN/HCN falls (S/N 0.74-0.83 at 1-0) and HCN does not fall alone. HNC/HCN at S is lower at 1-0 only (3-2 is equal within errors).

## What remains uncertain
- Composition versus excitation versus opacity is not separated.
  - HCN, HNC and CN all have higher n_crit than HCO+. But HNC, the lowest of the three, shows the larger drop at 1-0, so density alone is not enough.
  - There is no matched-beam H13CN/H13CO+ detection (Q29, Q42).
  - Lower HNC/HCN at S at 1-0 hints at warmer gas, while p0056 RADEX prefers a colder S. This tension is unresolved.
- The comparison is in the image plane only (no uv matching).
  - HNC 3-2 and CN 3-2 recovery of extended flux is unmeasured. The 3-2 HNC/HCN N/S depends on it.
  - I found one unidentified line at about +700 km/s near CN 3-2 and masked it from the baseline.
  - The 1-0 flux scales differ by project (HNC/HCN 0.52 vs 0.43 between 00993 and 00412, same at N and S), which cancels in N/S.
- The continuum at S may contain non-dust emission. That would make the HCN, HNC and CN deficit relative to dust larger, not smaller.
