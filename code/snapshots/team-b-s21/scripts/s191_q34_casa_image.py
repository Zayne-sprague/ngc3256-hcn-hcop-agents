"""Q34 CASA (host): image Band 6 HCN/HCO+ 3-2 visibilities (2017.1.00379.S Xb4d 12-m EB Xe07 + Xb4f 7-m EBs) on a common kλ range and common circular beams.
Variants per line and target beam (1.5", 2.2"): 'comb' = 12m+7m (shortest baseline ~7.7 kλ), 'c12' = 12m only (>=13 kλ, short-spacing-poor, ~ X94a-like),
'a7' = 7m only (sanity). uvtaper + restoringbeam = target circular beam, mosaic gridder (12m/7m primary beams), natural weights, auto-multithresh, pbcor, FITS export.
Cubes: 20 km/s channels, radio LSRK, 2180-3340 km/s. Run: rawdata.py casa scripts/s191_q34_casa_image.py"""
import os, time, glob, shutil
from casatasks import tclean, exportfits, impbcor
OUT='/workspace/out/q34/'; WK='/workspace/out/q34/work/'; os.makedirs(WK,exist_ok=True)
# /archive/derived is read-only: tclean needs a writable MS, so copy the (small) split MSs
MSD=OUT+'ms/'; os.makedirs(MSD,exist_ok=True)
def cp(pat):
  o=[]
  for m in sorted(glob.glob(pat)):
    d=MSD+os.path.basename(m)
    if not os.path.exists(d): shutil.copytree(m,d)
    o.append(d)
  return o
D12=cp('/archive/derived/team-b/20261026T165750-8b30/*.split.ms')
D07=cp('/archive/derived/team-b/20261026T173354-761f/*.split.ms')
LINES={'HCN32':'265.88618GHz','HCOp32':'267.55762GHz'}
BEAMS={'1p5':(1.5,1.1),'2p2':(2.2,1.95)}   # target beam, taper FWHM
VARS=[('comb',D12+D07),('c12',D12)]   # 7-m alone cannot reach a 1.5/2.2" beam (~5")
for bn,(bm,tp) in BEAMS.items():
  for vn,vis in VARS:
    for ln,rf in LINES.items():
      tag='%s_%s_%s'%(ln,bn,vn); fn=OUT+tag+'.pbcor.fits'
      if os.path.exists(fn): continue
      t0=time.time(); im=WK+tag
      for f in glob.glob(im+'.*'): shutil.rmtree(f,ignore_errors=True)
      try:
        tclean(vis=vis,imagename=im,field='ngc_3256',specmode='cube',gridder='mosaic',imsize=240,cell='0.3arcsec',phasecenter='ICRS 10:27:51.23 -43.54.16.6',
               restfreq=rf,start='2180km/s',width='20km/s',nchan=58,outframe='LSRK',veltype='radio',weighting='natural',uvtaper=['%.2farcsec'%tp],
               restoringbeam=['%.2farcsec'%bm],pblimit=0.2,niter=30000,usemask='auto-multithresh',sidelobethreshold=2.0,noisethreshold=4.0,lownoisethreshold=1.5,
               negativethreshold=0.0,minbeamfrac=0.3,nsigma=2.5,threshold='0mJy',interactive=False,perchanweightdensity=True,cycleniter=2000,parallel=False,pbcor=True)
        exportfits(im+'.image.pbcor',fn,overwrite=True,dropdeg=True)
        exportfits(im+'.pb',OUT+tag+'.pb.fits',overwrite=True,dropdeg=True)
        exportfits(im+'.image',OUT+tag+'.image.fits',overwrite=True,dropdeg=True)
        print('DONE',tag,'%.1f min'%((time.time()-t0)/60),flush=True)
      except Exception as e: print('FAIL',tag,repr(e),flush=True)
      for f in glob.glob(im+'.*'): shutil.rmtree(f,ignore_errors=True)
