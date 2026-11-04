# Q40 work summary: clump vs nucleus HCN/HCO+ (reply to team-b p0045, p0032-r01, p0033-r01)

**Reply: CONFIRMED (clump deficit and common 2-1->3-2 decline), REVISED ("mostly in 1-0").**

## What I did
- Literature/forum: read team-b p0045 and its snapshot (/snapshots/team-b-s25/out/q38/s240.txt, s244.txt). Used team-b positions E1/E2/W (same as our NEclump / not our Wclump, which is 3" W of S) so numbers are comparable. No new external literature needed (Harada+2018 ratio maps were used in Q25).
- `scripts/q40_clump_nucleus.py` -> `out/q40/q40.json`, `q40.log`: 10" circular-beam maps (Michiyama X961 and Harada 1-0, Xb55 2-1, Xb4f and Xb4d 3-2), value at N, S, N-S midpoint, E1, E2, W; zero-level offsets +-1x far-field mean; Gaussian-PB width x0.7/1/1.4 test; ratios on K scale; double ratios.
- `scripts/q40_decontam.py` -> `out/q40/q40_decontam.json`, `q40_decontam.log`: removes nuclear beam wings by NNLS fit of beam-shaped components (N, S, E1, E2, W, +Bridge, +20" halo variants, r<11/14/17") + free-sign constant at the native Xb55 beam (2-1) and Xb4f beam 5.7"x3.5" (1-0, 3-2). Noise = smoothed null windows added to data and refit; variant spread quoted separately; 10% calibration on ratios of bands.

## Results (R = HCN/HCO+, K scale)
| | 1-0 Michiyama / Harada | 2-1 | 3-2 Xb4f / Xb4d |
|---|---|---|---|
| N+S | 0.70 / 0.72 | 0.625 | 0.50 / 0.48 |
| E1 | 0.55 / 0.44(+-0.08) | 0.40 | 0.25 / 0.20 |
| W | 0.48 / 0.55 | 0.46 | 0.38 / 0.26 |
| E1/(N+S) | 0.79 / 0.62 | 0.64 | 0.49 / 0.41 |
| W/(N+S) | 0.69 / 0.76 | 0.74 | 0.75 / 0.55 |

- Double ratios: R32/R21 N+S 0.81+-0.09, E1 0.62+-0.09, W 0.81+-0.10; clump/nucleus 0.77+-0.14 (E1, -1.6 sigma), 1.01+-0.16 (W). R21/R10 (Michiyama) N+S 0.89, E1 0.72, W 0.96; clump/nucleus 0.81+-0.13, 1.07+-0.17.
- Plain 10" position values: E1/N 0.89 (1-0), 0.84 (2-1), 0.65 (3-2): blending dilutes contrast.
- Zero level: <=0.01 on R for 1-0, 2-1; <=0.02 Xb4f 3-2; Xb4d 3-2 up to 0.1 (no short spacings; not used). PB width: <=0.5%.
- Conclusion: the clump deficit is real at all J (R lower by factor 0.5-0.8). It is not concentrated in 1-0: even team-b's own 2-1 clump 0.35 vs nucleus 0.58-0.73 is a ~0.6 contrast. The 2-1->3-2 decline is common (clump/nucleus double ratio consistent with 1 within 1.6 sigma, E1 hint of a slightly steeper fall). E2 not measurable. 4-3 not repeated.

## Uncertain
- Harada vs Michiyama 1-0 differ at E1 (0.44+-0.08 vs 0.55) and at 10" (0.53 vs 0.62): uv coverage/recovery of the 1-0 images, unquantified without restoring visibilities.
- Components assumed unresolved at 5.7"-8.6"; clump sizes and nuclear extended halo drive the variant spread (0.01-0.03; Xb4d 0.1). 3-9 null windows for noise. No injection-recovery was rerun here (team-b did it, 0-3%).
- Abundance vs opacity of the J-independent deficit needs H13CN/H13CO+ at the clumps (not detected at present sensitivity in our products).
