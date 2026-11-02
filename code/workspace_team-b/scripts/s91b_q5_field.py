# Q5: empirical (spatial-field) noise of the matched-filter amplitude. For each region's line shape/centre, the single-line (+linear baseline) amplitude is computed
# in every pixel of the cube, averaged over the same beam-diameter aperture, and the robust scatter (1.4826 MAD) is taken over off-source pixels
# (5-13 arcsec from the centre, > 2.5 arcsec from every region position, within the field where the pb-corrected noise is usable). Includes correlated noise, baseline ripple, pb noise gradient.
import sys; sys.path.insert(0,'/workspace/scripts')
from s91_q5_lib import *
from common import files, cube_info, WCS
from regions import REGIONS
from scipy.ndimage import convolve
import json
CASES=[('Xb47','spw22',342.88285,'CS 7-6',700),('Xb47','spw20',344.31061,'SO 8_8-7_7',700),('Xb4f','spw22',251.82577,'SO 6_5-5_4',700),('Xb4d','spw31',251.82577,'SO 6_5-5_4',700)]
out={}
for ous,spw,f0,lab,win in CASES:
    fn=files(f'*{ous}*{spw}.cube.I.pbcor.fits')[0]; h,hd,w,nu=cube_info(fn); cube=np.array(h[0].data[0],dtype=np.float32)*1e3; ps=abs(hd['CDELT1'])*3600
    wc=WCS(hd).celestial; ny,nx=cube.shape[1:]; yy,xx=np.mgrid[:ny,:nx]
    x0,y0=[float(t) for t in wc.world_to_pixel_values(156.96346,-43.90461)]; rr=np.hypot(xx-x0,yy-y0)*ps
    z=np.load(f'/workspace/out/q5/spec_{ous}_{spw}.npz'); apr=float(z['apr_arcsec'])/ps; r=int(np.ceil(apr)); kyx=np.mgrid[-r:r+1,-r:r+1]; ker=(np.hypot(*kyx)<=apr).astype(float); ker/=ker.sum()
    mask=(rr>5)&(rr<13)
    for reg in REGIONS:
        ra,dec,_=REGIONS[reg]
        shp,vc,fw,snr=line_shape(ous,reg,mode='profile',tname='ref')
        v=vel(nu,f0)-vc; m=abs(v)<win; t=shp(v[m]); X=np.array([t,np.ones(m.sum()),v[m]/win]).T; P=np.linalg.pinv(X)
        amp=np.tensordot(P[0],np.nan_to_num(cube[m]),axes=(0,0))*1e-3         # Jy/beam km/s per pixel
        amp=convolve(amp,ker,mode='nearest')
        mk=mask.copy()
        for k2,(ra2,dec2,_) in REGIONS.items():
            xp,yp=[float(t) for t in wc.world_to_pixel_values(ra2,dec2)]; mk&=np.hypot(xx-xp,yy-yp)*ps>2.5
        a=amp[mk&np.isfinite(amp)&(amp!=0)]
        out[f'{ous}|{lab}|{reg}']=dict(mad=float(1.4826*np.median(abs(a-np.median(a)))),std=float(a.std()),median=float(np.median(a)),n=int(len(a)),p99=float(np.percentile(a,99)),max=float(a.max()))
    print(ous,lab,{k.split('|')[2]:round(v['mad'],3) for k,v in out.items() if k.startswith(f'{ous}|{lab}')},flush=True)
json.dump(out,open('/workspace/out/q5/field_noise.json','w'),indent=1)
