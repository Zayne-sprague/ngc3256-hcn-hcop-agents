# Q23 round 2: T_b vs beam with (a) 0.18" cube, (b) 0.82" cube convolved to b>=1" (full line coverage, incl. N red wing), (c) both + feather-like dT (7m, SD)
import numpy as np, json, sys, warnings, pandas as pd; warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/scripts')
from astropy import units as u
from radio_beam import Beam
import importlib.util
spec=importlib.util.spec_from_file_location('mf','/workspace/scripts/q23b_missing_flux.py')
# reuse helpers without re-running the main loop: re-define minimal pieces here
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve_fft
src=open('/workspace/scripts/q23b_missing_flux.py').read().split('res={}')[0]
exec(src)
mf=json.load(open('/workspace/out/q23b/missing_flux.json')); vp=json.load(open('/workspace/out/q23b/vmatch_pos.json'))
wd=pd.read_csv('/workspace/out/q5/q5e_width.csv')
D={k:load(k,0) for k in ['mid']}
VGf=np.arange(2250,3000,5.08)
out={}
for nm,(x0,y0) in CORES.items():
    r=wd[(wd.line=='HCN')&(wd.src==nm)&(wd.aperture_r==0)].iloc[0]; v0,fw=float(r.v0),float(r.fwhm)
    dT=np.interp(VGf,np.array(mf[nm]['VG']),np.array(mf[nm]['dT1_7m_minus_12m'])+np.array(mf[nm]['dT2_sd_minus_7m']))
    sub,cx,cy=cut(D['mid'],x0,y0,8.0)
    o={}
    for b in [1.0,1.7]:
        tgt=Beam(b*u.arcsec,b*u.arcsec,0*u.deg)
        if b<0.82: continue
        S=conv(sub,D['mid']['beam'],D['mid']['pix'],tgt)
        T=np.interp(VGf,D['mid']['v'],TRJ(interp_at(S,cx,cy),D['mid']['nu'],b,b))
        m=np.abs(VGf-v0)<=fw/2; g=np.exp(-0.5*((VGf-v0)/(fw/2.3548))**2)
        rec={}
        for lab,TT in [('12m_0.82cube',T),('+7m+SD',T+dT)]:
            rec[lab]={'peak':float(TT.max()),'v_peak':float(VGf[TT.argmax()]),'HCN_mean_FWHM':float(TT[m].mean()),'HCN_gauss_wt':float((TT*g).sum()/g.sum())}
        o[str(b)]=rec
    # 0.18 cube values (from vmatch_pos) with dT added (same velocity-grid mean)
    cb=vp['cores'][nm]['beams']
    h18={}
    for b,rec in cb.items():
        ci=rec['at_core_interp']
        h18[b]={'peak':ci['peak_RJ'],'HCN_mean_FWHM':ci['HCN']['T_mean_FWHM'],'cov':ci['HCN']['cov_FWHM'],'dT_mean_FWHM(7m+SD)':float(dT[np.abs(VGf-v0)<=fw/2].mean()),'dT_at_peak':float(dT[np.argmin(np.abs(VGf-ci['v_peak']))])}
    o['hi18']=h18; o['dense']=[v0,fw]
    out[nm]=o
    print(nm,json.dumps({k:v for k,v in o.items() if k!='hi18'},default=float)[:900]); print(nm,json.dumps(h18,default=float)[:1500])
json.dump(out,open('/workspace/out/q23b/corrected.json','w'),indent=1)
