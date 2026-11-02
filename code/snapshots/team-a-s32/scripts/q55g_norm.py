"""Q55g: absolute normalisation of the N/S contrast: line integrals (Jy km/s, same regions/windows as the ratios) divided by the 1.3 mm
dust-dominated continuum (Xb4d spw25_27_29_31 cont, pbcor) in the same circular beam (1.1" for 3-2, 1.0" for 4-3, 2.5" for 1-0).
N/S ratio of line/continuum for HCN, HCO+, HNC, CN. Continuum Tb and Td dependence stated in the report. Output out/q55/norm_cont.csv"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve
from scipy.ndimage import map_coordinates
import q4_lib as q
from matched import kernel
CF='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'
h=fits.open(CF); hd=h[0].header; img=np.nan_to_num(np.squeeze(h[0].data)).astype(np.float32); w=WCS(hd).celestial
nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); print('cont native beam',nb,'freq',hd.get('CRVAL3'))
def cont(reg,bm):
    tb=(bm,bm,0.); k=kernel(hd,tb,nb); sc=q.OM(tb)/q.OM(nb) if hasattr(q,'OM') else None
    a=fftconvolve(img,k.astype(np.float32),mode='same')*sc
    px,py=w.world_to_pixel_values(q.REG.loc[reg,'RA_deg'],q.REG.loc[reg,'Dec_deg'])
    return float(map_coordinates(a,[[py],[px]],order=1)[0])*1e3   # mJy/beam
import pandas as pd
M32=pd.read_csv('/workspace/out/q53/matched32.csv'); R=pd.read_csv('/workspace/out/q54/ratios.csv'); T43=pd.read_csv('/workspace/out/q53/matched43.csv')
rows=[]
def add(line,bm,Iv):
    cn,cs=cont('N',bm),cont('S',bm)
    rows.append(dict(line=line,beam=bm,I_N=Iv[0],I_S=Iv[1],cont_N_mJy=cn,cont_S_mJy=cs,NS_line=Iv[0]/Iv[1],NS_cont=cn/cs,NS_line_over_cont=(Iv[0]/cn)/(Iv[1]/cs)))
for bm in (1.1,1.4,1.9,2.5):
    m=M32[(M32.ous=='Xb4d')&(M32.beam_as==bm)]
    if len(m)==0: continue
    n=m[m.reg=='N'].iloc[0]; s=m[m.reg=='S'].iloc[0]
    add('HCN 3-2',bm,(n.IHCN,s.IHCN)); add('HCO+ 3-2',bm,(n.IHCO,s.IHCO))
    r=R[(R.num=='HNC3-2_00412')&(R.den=='HCN3-2_Xb4d')&(R.beam==bm)]
    if len(r)==2: add('HNC 3-2',bm,(r[r.reg=='N'].Inum.iloc[0],r[r.reg=='S'].Inum.iloc[0]))
    r=R[(R.num=='CN3-2_X1df1')&(R.den=='HCN3-2_Xb4d')&(R.beam==bm)]
    if len(r)==2: add('CN 3-2 (N=3-2 7/2-5/2 group)',bm,(r[r.reg=='N'].Inum.iloc[0],r[r.reg=='S'].Inum.iloc[0]))
for bm in (1.0,1.4,1.9,2.5):
    m=T43[T43.beam_as==bm]
    n=m[m.reg=='N'].iloc[0]; s=m[m.reg=='S'].iloc[0]
    add('HCN 4-3',bm,(n.IHCN,s.IHCN)); add('HCO+ 4-3',bm,(n.IHCO,s.IHCO))
V=pd.read_csv('/workspace/out/q54/vs_hco.csv')
r10=pd.read_csv('/workspace/out/q55/r10_core.csv')
for bm in (1.0,1.4,1.9,2.5):
    a=r10[r10.beam==bm]; add('HCN 1-0 (00965)',bm,(a[a.reg=='N'].IHCN.iloc[0],a[a.reg=='S'].IHCN.iloc[0])); add('HCO+ 1-0 (00965)',bm,(a[a.reg=='N'].IHCO.iloc[0],a[a.reg=='S'].IHCO.iloc[0]))
D=pd.DataFrame(rows); D.to_csv('/workspace/out/q55/norm_cont.csv',index=False)
pd.set_option('display.width',200); print(D.round(3).to_string())
