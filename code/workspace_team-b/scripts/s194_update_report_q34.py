"""Q34 s194: write Q34 results into /output/report.json (idempotent).  usage: python3 -I s194_update_report_q34.py"""
import json
R=json.load(open('/output/report.json')); R['results']=[x for x in R['results'] if not str(x.get('quantity','')).startswith('Q34')]
EV=('/workspace/scripts/s191_q34_casa_image.py (+.log); /workspace/scripts/s192_q34_measure.py -> /workspace/out/q34/s192.txt, s192.json; /workspace/scripts/s193_q34_radex.py -> /workspace/out/q34/s193.txt, s193.json; '
    'tapered cubes /workspace/out/q34/{HCN32,HCOp32}_{1p5,2p2}_{comb,c12}.pbcor.fits; /workspace/scripts/s189_q29_contrast_budget.py (previous image-plane budget)')
COM=['Band 6 side restored from raw visibilities: 2017.1.00379.S Xb4d 12-m EB uid___A002_Xda2c82_Xe07 (0.81 GB split MS) + 3 Xb4f 7-m EBs (0.28 GB); only 3 of the 14 Xb4f 7-m EBs and 1 of 2 Xb4d 12-m EBs, so the combined uv coverage is not the full archive one',
     'Band 3 side is NOT re-imaged from visibilities: restores of X961/X145 failed (no pipeline manifest) and X94a is 57.7 GB (>30 GB limit). 1-0 values are the short-spacing-complete archive cubes X145 (1.5") and X961 (2.2"), which agree to <=5% (Q29), convolved to the same circular beam; the match is therefore in beam and in short-spacing completeness, not an identical uv-range cut',
     'tclean: mosaic gridder (12m+7m primary beams), natural weights, uvtaper 1.1"/1.95" + restoringbeam 1.5"/2.2" circular, auto-multithresh, 20 km/s channels, pbcor; CASA reported "cannot convolve to new beam" warnings at rescaling in 6 of 8 runs but the delivered cubes carry the exact target beam',
     'p0034 conventions: aperture r = beam FWHM on the HCO+ 3-2 peak within 0.4" of nominal, +-350 km/s about 2780 km/s (and +-150), RJ T_b, null-free noise from line-free channels, D=2.8 sigma window dependence quoted separately; 10% absolute flux per band cancels in N/S contrasts, not in absolute r31 (14%)',
     'No Band 6 -> Band 3-coverage whole-disk control (taper to the 1-0 uv range) was possible without Band 3 visibilities']
ADD=[dict(quantity='Q34: uv-complete re-imaging of Band 6 HCN/HCO+ 3-2 at common circular beams (1.5", 2.2"): N/S contrast of HCN/HCO+ at 1-0 vs 3-2, D = contrast(3-2)/contrast(1-0)',
 target='N, S, Bridge; 1.5" (X145 1-0) and 2.2" (X961 1-0)',
 value=('2.2": 1-0 N/S 1.467+-0.022; 3-2 (Xb4d 12m + Xb4f 7m, restored from visibilities) 1.441+-0.023 vs 1.650 from the Xb4d image alone (12m only: 1.451); D = 0.983+-0.022 (+-350 km/s), 1.038+-0.016 (+-150). '
        '1.5": 1-0 1.616+-0.059; 3-2 1.532+-0.026 (Xb4d image 1.809); D = 0.948+-0.038 (+-350), 1.011+-0.029 (+-150). N/Bridge D=0.965+-0.018 (2.2"), 0.86+-0.03 (1.5"); S/Bridge D=0.98+-0.02 (2.2"). '
        'ADOPTED: D = 1.00+-0.04 (stat 0.02, +-350 vs +-150 window 0.03, residual systematic 0.02), i.e. the 1-0 and 3-2 N/S HCN/HCO+ contrasts agree. The previous D=1.10+-0.06 (p0034) is an artefact of the Xb4d 12-m 3-2 flux loss, which is larger at S (flux gain from adding the 7 m array at 2.2": HCN N +6%, S +34%, Bridge +15%; HCO+ +7%, +17%, +10%) and had inflated the 3-2 contrast. '
        'The 3-2 N/S contrast at 2.2" falls from 1.65 to 1.44 once short spacings are included; the 3-2 N/S contrast quoted at 50-120 pc (Q26, 2.1+-0.2) is a higher-resolution quantity and is not changed by this test.'),
 unit='dimensionless',uncertainty='D: +-0.04 adopted; 3-2 contrast stat 0.02-0.03',method='s191 (CASA tclean from calibrated split MS) -> s192 aperture measurement with p0034 conventions; compared with s183 1-0 values and the Xb4d image-plane values.',
 evidence=EV,confidence='medium-high (Band 6 matched, Band 3 image products)',assumptions=COM),
 dict(quantity='Q34: uv-complete r31 = T_b(3-2)/T_b(1-0) for HCN and HCO+ at N, S, Bridge, and N/S of r31',
 target='N, S, Bridge, 2.2" and 1.5" circular beams',
 value=('2.2" (+-350): HCN N 0.263+-0.002, S 0.261+-0.004, Bridge 0.258+-0.003; HCO+ N 0.351+-0.002, S 0.343+-0.003, Bridge 0.333+-0.003 (Xb4d-image values were HCN 0.247/0.194/0.225, HCO+ 0.330/0.292/0.303). '
        '1.5": HCN N 0.275, S 0.279, Bridge 0.252; HCO+ N 0.390, S 0.376, Bridge 0.307. N/S of r31 (flux scales cancel): HCN 1.007+-0.020 (2.2"), 0.98+-0.03 (1.5"); HCO+ 1.025+-0.010, 1.04+-0.02 (+-150 km/s: HCN 1.06+-0.014 / 1.05, HCO+ 1.02 / 1.04). '
        'Previous image-plane values 1.27/1.09 (HCN/HCO+) are therefore NOT supported: with short spacings r31 is the same at N and S within 0-6%, no evidence for an excitation difference between the nuclei; the Bridge has the lowest r31 only in HCO+ at 1.5" (0.31 vs 0.38-0.39) and is equal at 2.2" (0.33-0.34). '
        'r31(HCN) < r31(HCO+) everywhere by 25% (HCN 3-2/1-0 = 0.26, HCO+ 0.34-0.35). Absolute r31 +-14% (10% flux per band).'),
 unit='dimensionless (T_b ratio)',uncertainty='absolute +-14% (flux scale); N/S +-0.02-0.03 (stat) +-0.03 (window)',method='s192',evidence=EV,confidence='medium-high for N/S; medium for absolute',assumptions=COM),
 dict(quantity='Q34: flux recovery of the Band 6 12-m (Xb4d) 3-2 image at N, S, Bridge relative to 12m+7m at matched beams (two-sided measured recovery)',
 target='N, S, Bridge; r=beam and r=6" apertures; 1.5" and 2.2"',
 value=('combined(12m+7m)/Xb4d image flux at 2.2" r=2.2": HCN N 1.06, S 1.34, Bridge 1.15; HCO+ N 1.07, S 1.17, Bridge 1.10 (1.5": HCN 1.00/1.25/1.24; HCO+ 1.04/1.11/1.10). 12m-only retapered/combined: r=beam 0.89-0.94 (2.2"), 0.91-0.96 (1.5"); r=6" 0.81-0.86. '
        'This reproduces the earlier Xb4f/Xb4d bound at S (1.30-1.36) for HCN, and is a measurement rather than a one-sided limit.'),
 unit='ratio',uncertainty='+-0.01-0.02 stat; 12m EB subset: 1 of 2',method='s192',evidence=EV,confidence='medium-high',assumptions=COM),
 dict(quantity='Q34: excitation vs abundance for N vs S with the uv-complete 3-2 (RADEX s186 framework)',
 target='N, S, Bridge; 2.2", 1.5"',
 value=('Posteriors (1-0+3-2) practically unchanged: N Tk 22-40 K (broad 11-158), log n 3.9-4.0, log A=N(HCN)/N(HCO+) 0.70 (5.0); S Tk 112-141 K (35-251), log n 4.0, log A 0.40 (2.5); Bridge Tk 80-112 K, log A 0.5-0.6. The HCN/HCO+ column ratio N/S (x2, 0.3 dex) stays, as does the lower N Tk, but both are driven by the 1-0 HCN/HCO+ ratio contrast and by the (1-0 vs 3-2) HCN and HCO+ line opacity, weakly. '
        'ln BF N vs S (vs free), 2.2", 1-0+3-2: same Tk,n with free columns (abundance-only) +1.02; same A, free Tk,n (excitation-only) +0.21; excitation+column only (same A and N) +1.31 (old Xb4d-based: +1.30, -0.04, +0.97); 1.5": +0.26, -0.42, -0.10 (old +0.61, -1.20, -1.04). '
        'The difference abundance-only minus excitation-only is +0.8 (2.2") and +0.7 (1.5"): at best weak (<1 in lnBF) preference for a common density/temperature with different HCN/HCO+ abundance at N vs S over a common abundance with different excitation, not decisive. '
        'A full "same everything" (filling only) is disfavoured ln BF -4.6 (2.2") and -9.8 (1.5"), far less than before (-14, -21), because the Xb4d recovery bias at S is gone. If the data are read as within 2 sigma (D=1.00+-0.04, r31 N/S=1.0), nothing in the data requires different excitation; N/S differences are in the HCN/HCO+ ratio at both J, which is an abundance (or opacity/column) signal.'),
 unit='ln Bayes factor, K, cm-3 dex',uncertainty='grid quantization 0.1 dex; ln BF depends on band-scale and per-line error priors (Q29 sensitivity: +-1)',method='s193 (s186 machinery, 3-2 replaced by s192 uv-complete values)',evidence=EV,confidence='medium-low (priors; not decisive)',assumptions=COM)]
for a in ADD:
    a['assumptions']=[dict(name='Assumption %d'%(i+1),value=t) for i,t in enumerate(COM)]
    a['confidence']={'medium-high':'medium'}.get(a['confidence'].split(' ')[0],a['confidence'].split(' ')[0]) if not a['confidence'].startswith('medium-low') else 'low'
    a['confidence']=a['confidence'] if a['confidence'] in ('high','medium','low','tentative') else 'medium'
R['results'].extend(ADD); json.dump(R,open('/output/report.json','w'),indent=1)
