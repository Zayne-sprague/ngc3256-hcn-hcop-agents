# Q22: spw-edge response across all OUS/spws: mean pbcor peak-pixel level in first/last 10 channels (mJy/b), minus the local reference (ch 60-120 median from same edge), at N,S,Bridge,W; continuum at region from line-free cont image.
import numpy as np, glob, json, os, warnings, logging
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
from astropy.io import fits
from astropy.wcs import WCS
res={}
def cont(ous,ra,dec):
    f=glob.glob(f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.spw*_*.cont.I.pbcor.fits')[0]
    h=fits.open(f); d=np.squeeze(h[0].data); w=WCS(h[0].header).celestial
    x,y=w.world_to_pixel_values(ra,dec); return float(d[int(round(float(y))),int(round(float(x)))])*1e3
print('%-5s %-3s %-6s %7s | %s'%('ous','spw','reg','cont','lowedge(0-10) lowref(60-120) | highedge(0-10) hiref'))
for fn in sorted(glob.glob('/workspace/out/regspec/*.npz')):
    ous,spw=os.path.basename(fn).split('.')[:2]; z=np.load(fn)
    nu=z['nu']; nu=nu/1e9 if nu.max()>1e6 else nu
    for r in ['N','S','Bridge','W']:
        if r+'_pix' not in z.files: continue
        y=z[r+'_pix']*1e3; ra,dec=z[r+'_xy']; c=cont(ous,ra,dec)
        lo=y[:10].mean(); lr=np.median(y[60:120]); hi=y[-10:].mean(); hr=np.median(y[-120:-60]); n=np.nanstd(y[len(y)//2-60:len(y)//2+60])
        res[f'{ous}.{spw}.{r}']=dict(cont=c,lo=float(lo),lo_ref=float(lr),hi=float(hi),hi_ref=float(hr),rms=float(n),f0=float(nu[0]),f1=float(nu[-1]))
        print(f'{ous:5s} {spw:5s} {r:6s} {c:7.2f} | {lo-lr:6.2f} {lr:6.2f} | {hi-hr:6.2f} {hr:6.2f}  rms~{n:.1f}  [{nu[0]:.1f}-{nu[-1]:.1f}]')
json.dump(res,open('/workspace/out/q22/s131.json','w'),indent=0)
