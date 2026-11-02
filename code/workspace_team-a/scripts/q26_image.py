# CASA 6 script (run via rawdata.py casa / --then). Re-images restored 2017.1.00379.S visibilities of
# HCN/HCO+ 2-1 (Xb55, 7-m), 3-2 (Xb4f 7-m; Xb4d 12-m compact) with a COMMON uv range and taper and a common
# restoring beam, cleaned (multiscale) so the u>uvmin filtering does not leave a dirty-beam bowl in the annuli.
import glob, json, os, itertools
from casatasks import tclean, exportfits, rmtables
OUT = '/workspace/out/q26/casa'
os.makedirs(OUT, exist_ok=True)
groups = {}
for s in glob.glob('/archive/derived/*/*/summary.json'):
    j = json.load(open(s)); d = os.path.dirname(s)
    ous = j['ous'].split('/')[-1]
    if ous in ('Xb55', 'Xb4f', 'Xb4d'):
        groups.setdefault(ous, []).extend(sorted(glob.glob(d + '/*.split.ms')))
print(groups)
LINES = {'HCN': {'2-1': 177.2611, '3-2': 265.8864}, 'HCO+': {'2-1': 178.3751, '3-2': 267.5576}}
JOBS = [('Xb55', '2-1'), ('Xb4f', '3-2'), ('Xb4d', '3-2')]
UVMIN = [6, 8]
TAPER = '8arcsec'; BEAM = '10arcsec'
for ous, J in JOBS:
    vis = sorted(set(groups.get(ous, [])))
    if not vis: print('no data', ous); continue
    for sp, rf in LINES.items():
        for uvmin in UVMIN:
            name = f'{OUT}/{ous}_{sp.replace("+","p")}_{J}_u{uvmin}'
            rmtables(name + '*')
            tclean(vis=vis, imagename=name, specmode='cube', restfreq=f'{rf[J]}GHz', outframe='LSRK',
                   veltype='radio', start='1775km/s', width='20km/s', nchan=100,
                   imsize=[96, 96], cell='2arcsec', weighting='natural', uvrange=f'>{uvmin}klambda',
                   uvtaper=[TAPER], restoringbeam=BEAM, gridder='standard', pblimit=0.2,
                   deconvolver='multiscale', scales=[0, 3, 8], niter=3000, threshold='0mJy', nsigma=3.0,
                   usemask='auto-multithresh', sidelobethreshold=2.0, noisethreshold=4.0,
                   interactive=False, pbcor=True, restart=True)
            for ext in ('image', 'image.pbcor', 'pb', 'psf', 'residual'):
                try: exportfits(name + '.' + ext, name + '.' + ext + '.fits', overwrite=True, dropdeg=True)
                except Exception as e: print('export fail', ext, e)
