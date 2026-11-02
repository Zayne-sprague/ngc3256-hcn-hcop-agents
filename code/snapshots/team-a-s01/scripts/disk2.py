# Disk test: is the HCO+ 2-1 (Xb55, 9.4x6.3") emission at r>8" just the beam wing of the compact emission, or genuinely extended gas?
# Convolve the Xb4f HCO+ 3-2 (and Xb4d) moment-0 maps to the Xb55 beam, regrid to the Xb55 grid, scale to the central (r<6") 2-1 flux, compare residuals.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from matched import kernel
from scipy.signal import fftconvolve
from scipy.ndimage import map_coordinates
import pandas as pd
RA0,DE0=156.96346,-43.90461
def pos(hd):
    w=WCS(hd).celestial; ny,nx=hd['NAXIS2'],hd['NAXIS1']; Y,X=np.mgrid[:ny,:nx]; ra,de=w.pixel_to_world_values(X,Y)
    return (ra-RA0)*np.cos(np.deg2rad(DE0))*3600,(de-DE0)*3600
hB=fits.getheader(f('Xb55','spw18')); TB=(hB['BMAJ']*3600,hB['BMIN']*3600,hB['BPA']); dxB,dyB=pos(hB); rB=np.hypot(dxB,dyB)
zB=np.load('/workspace/out/m0_Xb55_spw18.npz'); mB,eB=zB['m0'],zB['e0']
res={}
for ous,spw,nm in [('Xb4f','spw16','HCO+ 3-2 Xb4f'),('Xb4d','spw25','HCO+ 3-2 Xb4d')]:
    hd=fits.getheader(f(ous,spw)); nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
    if ous=='Xb4d':
        import glob
        z=np.load('/workspace/out/m0_Xb4d_spw25.npz')
    else: z=np.load(f'/workspace/out/m0_{ous}_{spw}.npz')
    m0=np.nan_to_num(z['m0']); k=kernel(hd,TB,nb)
    om=lambda b: np.pi*b[0]*b[1]/(4*np.log(2))
    sm=fftconvolve(m0,k,mode='same')*om(TB)/om(nb)
    w=WCS(hd).celestial; x,y=w.world_to_pixel_values(RA0+dxB/3600/np.cos(np.deg2rad(DE0)),DE0+dyB/3600)
    reg=map_coordinates(sm,[y,x],order=1,mode='nearest')
    bowlB=np.nanmedian(mB[(rB>14)&(rB<22)])
    c=(rB<6)&np.isfinite(mB)
    s=np.nansum((mB-bowlB)[c])/np.nansum(reg[c]); resid=(mB-bowlB)-s*reg
    out=[]
    for a,b in [(6,8),(8,10),(10,12),(12,14),(14,18)]:
        m=(rB>=a)&(rB<b)&np.isfinite(mB)
        out.append((a,b,np.nanmean((mB-bowlB)[m]),np.nanmean(s*reg[m]),np.nanmean(resid[m]),np.nanstd(resid[(rB>14)&(rB<22)])/np.sqrt(max(1,m.sum()/5))))
    print(nm,'scale 2-1/3-2 =',round(s,3)); print(pd.DataFrame(out,columns=['r_in','r_out','obs_2-1_bowlcorr','scaled_conv_3-2','resid','e_resid']).round(2).to_string())
    # the SW island
    isl=(rB>8)&(np.hypot(dxB+7.2,dyB+3.6)<5)
    print('  SW island (r>8, within 5" of (-7.2,-3.6)): obs %.2f, model %.2f, resid %.2f Jy/b km/s'%(np.nanmean((mB-bowlB)[isl]),np.nanmean(s*reg[isl]),np.nanmean(resid[isl])))
