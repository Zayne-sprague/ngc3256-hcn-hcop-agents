# Q5: spectra at FIXED sky positions (same for every spw of an OUS; no per-cube peak picking, which biases weak lines positive).
# Position of each region = refined peak in the HCO+ (J=4-3 Xb47, J=3-2 Band 6) cube of that OUS (out/regspec/*_xy). Aperture = circle of diameter = geometric-mean
# beam FWHM of the reference cube (mean in pbcor Jy/beam), identical for all spws of the OUS. Also a peak-pixel (nearest pixel) spectrum.
# usage: python3 -I s90_q5_extract.py
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
REF={'Xb47':'spw18','Xb4f':'spw16','Xb4d':'spw25','Xb4b':'spw25'}
SPW={'Xb47':['16','18','20','22'],'Xb4f':['16','18','20','22'],'Xb4d':['25','27','29','31'],'Xb4b':['25','27','29','31']}
for ous in ('Xb47','Xb4f','Xb4d','Xb4b'):
    ref=np.load(f'/workspace/out/regspec/{ous}.{REF[ous]}.cube.pbcor.npz')
    apr_as=0.5*np.sqrt(float(ref['bmaj'])*float(ref['bmin']))
    pos={k:ref[k+'_xy'] for k in REGIONS if k+'_xy' in ref}
    for sp in SPW[ous]:
        fn=files(f'*{ous}*spw{sp}.cube.I.pbcor.fits')[0]
        h,hd,w,nu=cube_info(fn); wc=WCS(hd).celestial; ps=abs(hd['CDELT1'])*3600; nz=len(nu); ny,nx=hd['NAXIS2'],hd['NAXIS1']
        out={'nu':nu,'bmaj':hd['BMAJ']*3600,'bmin':hd['BMIN']*3600,'apr_arcsec':apr_as}
        for k,(ra,dec) in pos.items():
            x,y=[float(t) for t in wc.world_to_pixel_values(ra,dec)]; r=apr_as/ps; m=int(np.ceil(r+2))
            if x<m or y<m or x>nx-m or y>ny-m: continue
            x0,x1,y0,y1=int(x)-m,int(x)+m+1,int(y)-m,int(y)+m+1
            patch=np.array(h[0].data[0,:,y0:y1,x0:x1],dtype=np.float32); yy,xx=np.mgrid[y0:y1,x0:x1]
            ap=np.hypot(xx-x,yy-y)<=r
            out[k+'_ap']=np.nanmean(patch[:,ap],axis=1); out[k+'_pix']=patch[:,int(round(y))-y0,int(round(x))-x0]
        np.savez(f'/workspace/out/q5/spec_{ous}_spw{sp}.npz',**out); print(ous,sp,[k for k in pos if k+'_ap' in out],f'apr={apr_as:.2f}"',flush=True)
