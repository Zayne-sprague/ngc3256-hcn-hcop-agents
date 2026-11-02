# Objective "extended disk" definition. Mask = S/N>5 islands (after subtracting the local zero-level 'bowl', median of the r=14-22" ring) in the Xb4f (5.7x3.5") HCO+ 3-2 +-350 km/s map,
# lying outside R_in = 8" (all compact N/S/Bridge/NE/E peaks of Xb4d at S/N>8 are inside 8"). Islands, their integrated flux, and stacked HCO+/HCN/CH3OH spectra are reported;
# if empty: noise-level limit in the 8-14" annulus.  Output out/disk_islands.csv, out/disk_stack.npz, out/disk_stack.csv
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy import ndimage as ndi
import pandas as pd
RA0,DE0=156.96346,-43.90461; RIN=8.; C=299792.458
def grid(ous,spw):
    hd=fits.getheader(f(ous,spw)); w=WCS(hd).celestial; ny,nx=hd['NAXIS2'],hd['NAXIS1']; Y,X=np.mgrid[:ny,:nx]
    ra,de=w.pixel_to_world_values(X,Y); dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600
    return hd,w,dx,dy,np.hypot(dx,dy)
rows=[]; out={}
for ous,spw,name in [('Xb4f','spw16','HCO+ 3-2'),('Xb4f','spw18','HCN 3-2'),('Xb55','spw18','HCO+ 2-1'),('Xb55','spw16','HCN 2-1'),('Xb47','spw18','HCO+ 4-3')]:
    z=np.load(f'/workspace/out/m0_{ous}_{spw}.npz'); m0,e0=z['m0'],z['e0']; hd,w,dx,dy,r=grid(ous,spw)
    ring=(r>14)&(r<22)&np.isfinite(m0); bowl=np.nanmedian(m0[ring]); e_ring=1.4826*np.nanmedian(np.abs(m0[ring]-bowl))
    sn=(m0-bowl)/np.sqrt(e0**2+0*e_ring**2); sn_emp=(m0-bowl)/max(e_ring,np.nanmedian(e0[r<14]))
    mask=(sn_emp>5)&(r>RIN)&(r<20)&np.isfinite(m0)
    lab,n=ndi.label(mask)
    pix=abs(hd['CDELT2'])*3600
    for i in range(1,n+1):
        m=lab==i; ys,xs=np.where(m); k=np.argmax(np.where(m,m0,-1e9)); yk,xk=np.unravel_index(k,m0.shape)
        rows.append(dict(ous=ous,line=name,island=i,npix=int(m.sum()),area_arcsec2=float(m.sum()*pix**2),peak=float(m0[yk,xk]-bowl),SN_emp=float(sn_emp[yk,xk]),dRA=float(dx[yk,xk]),dDec=float(dy[yk,xk]),r=float(r[yk,xk]),sum_Jykms_pix=float(np.nansum(m0[m]-bowl)),bowl=float(bowl),e_ring=float(e_ring)))
    ann=(r>=RIN)&(r<14)&np.isfinite(m0)
    rows.append(dict(ous=ous,line=name,island=0,npix=int(ann.sum()),peak=float(np.nanmean(m0[ann]-bowl)),SN_emp=float(np.nanmean(m0[ann]-bowl)/(e_ring/np.sqrt(max(1,ann.sum()/ (np.pi*hd['BMAJ']*hd['BMIN']*3600**2/4/np.log(2)/pix**2))))),dRA=np.nan,dDec=np.nan,r=np.nan,bowl=float(bowl),e_ring=float(e_ring)))
    out[f'{ous}_{spw}_sn']=sn_emp
df=pd.DataFrame(rows); pd.set_option('display.width',250,'display.max_rows',200); print(df.round(2).to_string()); df.round(3).to_csv('/workspace/out/disk_islands.csv',index=False)
