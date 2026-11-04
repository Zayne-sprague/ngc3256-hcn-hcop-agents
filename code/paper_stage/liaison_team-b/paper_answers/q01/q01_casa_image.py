"""q01 (paper stage) CASA (host): uv-matched Band 6 HCN/HCO+ 3-2 cubes at circular 0.6" and 1.0" beams, 2017.1.00379.S.
Visibilities: Xb4b 12-m EB Xd271e2_Xb367 (0.25") + Xb4d 12-m EB Xe07 (0.8") + Xb4f 7-m EBs X1cf5, X2a80, X42f0.
Variants: 'comb' = all arrays (uv >~7.7 kl, includes the 7-m), 'c12' = 12-m only (uv >~13 kl, short-spacing-poor).
Same recipe as s191_q34_casa_image.py (mosaic gridder, natural weights, uvtaper + restoringbeam = circular target beam, auto-multithresh, 20 km/s channels, pbcor).
Run: rawdata.py casa paper_answers/q01/q01_casa_image.py   (also run automatically after the Xb4b restore, --then)"""
import os, time, glob, shutil
from casatasks import tclean, exportfits
OUT='/workspace/out/q01/'; WK=OUT+'work/'; MSD=OUT+'ms/'
for d in (OUT,WK,MSD): os.makedirs(d,exist_ok=True)
def cp(pat):
  o=[]
  for m in sorted(glob.glob(pat)):
    d=MSD+os.path.basename(m)
    if not os.path.exists(d): shutil.copytree(m,d)
    o.append(d)
  return o
B4B=cp('/archive/derived/team-b/*/uid___A002_Xd271e2_Xb367.split.ms')
D12=cp('/archive/derived/team-b/20261026T165750-8b30/*.split.ms')
D07=cp('/archive/derived/team-b/20261026T173354-761f/*.split.ms')
print('MS:',B4B,D12,D07,flush=True)
LINES={'HCN32':'265.88618GHz','HCOp32':'267.55762GHz'}
# (target circular beam, taper FWHM): taper ~ sqrt(beam^2 - nat^2), nat ~0.3" for Xb4b+Xb4d natural weights
BEAMS={'0p6':(0.6,0.52),'1p0':(1.0,0.95)}
VARS=[('comb',B4B+D12+D07),('c12',B4B+D12)]
for bn,(bm,tp) in BEAMS.items():
  for vn,vis in VARS:
    for ln,rf in LINES.items():
      tag='%s_%s_%s'%(ln,bn,vn); fn=OUT+tag+'.pbcor.fits'
      if os.path.exists(fn): continue
      t0=time.time(); im=WK+tag
      for f in glob.glob(im+'.*'): shutil.rmtree(f,ignore_errors=True)
      try:
        tclean(vis=vis,imagename=im,field='ngc_3256',specmode='cube',gridder='mosaic',imsize=200,cell='0.1arcsec',phasecenter='ICRS 10:27:51.23 -43.54.16.6',
               restfreq=rf,start='2180km/s',width='20km/s',nchan=58,outframe='LSRK',veltype='radio',weighting='natural',uvtaper=['%.2farcsec'%tp],
               restoringbeam=['%.2farcsec'%bm],pblimit=0.2,niter=30000,usemask='auto-multithresh',sidelobethreshold=2.0,noisethreshold=4.0,lownoisethreshold=1.5,
               negativethreshold=0.0,minbeamfrac=0.3,nsigma=2.5,threshold='0mJy',interactive=False,perchanweightdensity=True,cycleniter=2000,parallel=False,pbcor=True)
        exportfits(im+'.image.pbcor',fn,overwrite=True,dropdeg=True)
        exportfits(im+'.pb',OUT+tag+'.pb.fits',overwrite=True,dropdeg=True)
        exportfits(im+'.image',OUT+tag+'.image.fits',overwrite=True,dropdeg=True)
        print('DONE',tag,'%.1f min'%((time.time()-t0)/60),flush=True)
      except Exception as e: print('FAIL',tag,repr(e),flush=True)
      for f in glob.glob(im+'.*'): shutil.rmtree(f,ignore_errors=True)
