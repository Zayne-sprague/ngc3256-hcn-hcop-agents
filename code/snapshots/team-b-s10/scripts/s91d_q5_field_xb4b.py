# Q5: off-source scatter of the SO 6_5-5_4 aperture-matched amplitude for Xb4b (0.25") by random aperture positions (5-13" from centre, >2.5" from region positions), S and N shapes.
import sys; sys.path.insert(0,'/workspace/scripts')
from s91_q5_lib import *
from common import files, cube_info, WCS
from regions import REGIONS
import json
rng=np.random.default_rng(3); out={}
fn=files('*Xb4b*spw31.cube.I.pbcor.fits')[0]; h,hd,w,nu=cube_info(fn); ps=abs(hd['CDELT1'])*3600; wc=WCS(hd).celestial
z=np.load('/workspace/out/q5/spec_Xb4b_spw31.npz'); apr=float(z['apr_arcsec'])/ps
x0,y0=[float(t) for t in wc.world_to_pixel_values(156.96346,-43.90461)]
pos=[]
reg_pix=[[float(t) for t in wc.world_to_pixel_values(r[0],r[1])] for r in REGIONS.values()]
while len(pos)<600:
    r=rng.uniform(5,13)/ps; a=rng.uniform(0,2*np.pi); x=x0+r*np.cos(a); y=y0+r*np.sin(a)
    if all(np.hypot(x-p[0],y-p[1])*ps>2.5 for p in reg_pix): pos.append((x,y))
m=int(np.ceil(apr))+1; yy,xx=np.mgrid[-m:m+1,-m:m+1]; ker=np.hypot(xx,yy)<=apr
spectra=[]
for x,y in pos:
    xi,yi=int(round(x)),int(round(y)); patch=np.array(h[0].data[0,:,yi-m:yi+m+1,xi-m:xi+m+1],dtype=np.float32); spectra.append(np.nanmean(patch[:,ker],axis=1)*1e3)
spectra=np.array(spectra)
for reg in ('S','N'):
    shp,vc,fw,snr=line_shape('Xb4b',reg,mode='profile',tname='ref'); v=vel(nu,251.82577)-vc; mm=abs(v)<700
    t=shp(v[mm]); X=np.array([t,np.ones(mm.sum()),v[mm]/700]).T; P=np.linalg.pinv(X)
    a=np.nan_to_num(spectra[:,mm])@P[0]*1e-3
    out[reg]=dict(mad=float(1.4826*np.median(abs(a-np.median(a)))),std=float(a.std()),median=float(np.median(a)),n=len(a),frac_above_0p043=float((a>=0.043).mean()))
print(out); json.dump(out,open('/workspace/out/q5/field_noise_xb4b_SO.json','w'))
