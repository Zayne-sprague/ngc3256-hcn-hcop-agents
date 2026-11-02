# Compare the pipeline clean mask (Xb4d cubes) with the emission regions: fraction of the |v-VSYS|<350 channels at each region (beam ap centre) inside the mask,
# and in the Sblue wing channels (-400..-150); also the line-integrated signal outside the mask in the 20"x20" field.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
for spw,f0 in (('spw25',267.55763),('spw27',265.88618),('spw29',None),('spw31',None)):
    fm=files(f'*Xb4d.ngc_3256_sci.{spw}.cube.I.mask.fits')[0]; fp=files(f'*Xb4d.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]
    hm,hd,w,nu=cube_info(fm); hp=fits.open(fp,memmap=True)[0]
    wc=WCS(hd).celestial; msk=np.squeeze(hm[0].data)
    print(f'## Xb4d {spw}: mask fraction (all) = {np.nanmean(msk>0.5):.4f}; channels with any mask: {int((msk.reshape(msk.shape[0],-1).max(1)>0.5).sum())}/{msk.shape[0]}')
    if f0 is None: continue
    v=vrad(nu,f0*1e9)-VSYS
    for k,(ra,dec,_) in REGIONS.items():
        x,y=wc.world_to_pixel_values(ra,dec); x,y=int(round(float(x))),int(round(float(y)))
        col=msk[:,y-2:y+3,x-2:x+3].max(axis=(1,2))>0.5
        a=(abs(v)<350); b=(v>=-400)&(v<=-150); c=(v>=150)&(v<=400)
        print(f'  {k:7s} masked frac core {col[a].mean():.2f} blue-wing {col[b].mean():.2f} red-wing {col[c].mean():.2f}')
    # signal outside mask: integrate pbcor over |v|<350 where mask==0 vs where mask==1 (full image)
    sel=np.where(abs(v)<350)[0]; cube=np.array(hp.data[0,sel[0]:sel[-1]+1],dtype=np.float32); mm=msk[sel[0]:sel[-1]+1]>0.5
    dv=abs(np.median(np.diff(v)))
    inside=np.nansum(np.where(mm,cube,0))*dv; outside=np.nansum(np.where(~mm,cube,0))*dv
    print(f'  sum over |v|<350: inside mask {inside:.1f}, outside mask {outside:.1f} Jy/beam km/s*pix (noise-dominated outside)')
