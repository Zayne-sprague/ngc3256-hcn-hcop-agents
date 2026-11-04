# Q26 work summary (round 3)
Question: is HCN/HCO+ at 8-12" independent of J? Done in the uv plane (matched kilo-lambda), restored raw visibilities (Xb55 2-1; Xb4f, Xb4d 3-2).

Scripts: scripts/q26_uvdump.py, q26_uvana.py (binning+fit, existing), NEW q26_inject.py (injection-recovery at real uv/sigma, 60 realisations), q26_boot.py (row bootstrap + EB jackknife, profile variants), q26_radex.py (expected R32/R21), build_report_q26b.py (report). Outputs: out/q26/q26_inject.{log,json}, q26_boot.{log,json}, q26_radex.log.

Results
- Injection, no fall in extended ratio (0.54 both J), cores with measured per-J ratios: binned R32/R21 = 0.91/0.91/0.89/0.86 (6-8/8-10/10-14/14-20 kl). Core leakage biases the binned ratio ~10% low. Injected 15% fall: 0.84/0.86/0.85/0.83. Observed 0.80/0.99/0.84/0.77 (mean 0.85): ~1 sigma from either; the earlier "15% fall" is mostly the cores.
- Core-removed parametric fit recovers R_ext ratio unbiased: 1.01+-0.08 (truth 1), 0.89+-0.08 (truth 0.85).
- Observed R_ext(3-2 Xb4f)/R_ext(2-1) = 1.05 (10", u>8), 1.01 (u>10), 0.87/0.77 (6"), 1.26/1.38 (15"); baseline-bootstrap sd 0.18-0.26; Xb4d (u>10 kl, chi2_r 4-10) 0.6-0.7. Adopted 1.0 +0.3/-0.4: consistent with no change and with a 15-25% fall.
- RADEX: expected R32/R21 = 0.5-0.93, typically 0.7-0.85 (sub-thermal gas), so a modest fall is expected physically; the data cannot discriminate.
Literature: none published for this target (data unpublished); Q22/Q25 internal numbers used for core ratios.
Uncertain / not done: truth model uses Gaussian+points (real chi2_r 2-10: clumpy disk); noise injected with formal sigma only (real scatter 2-5x larger, hence the bootstrap errors); Xb4d lacks 2-10 kl so can't be matched; velocity-window (+-250/450) and off-line control windows not run; injected truth profile fixed at 10" FWHM (profile mismatch tested only on the data side). Answer needing more: a short-spacing-matched 3-2 EB set (ACA/7-m at 3-2 with baselines 4-10 kl) and the Xb4b 0.25" images smoothed to test N-core ratio.
