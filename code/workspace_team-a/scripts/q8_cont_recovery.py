# Q8-3: recovery proxy from the line-free continuum images (pbcor 2D): all smoothed to the same beam, sampled at N/S/Bridge.
# Xb4d vs Xb4f at the SAME frequency (257.7 GHz) = pure array-recovery ratio of the continuum (no spectral-index term).
# Xb47 (346.5 GHz) / Xb4f gives rec47/rec4f * (346.5/257.7)^alpha (alpha unknown, compared N vs S only under alpha_N=alpha_S: ASSUMPTION).
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *; from matched import kernel
from scipy.signal import fftconvolve
import pandas as pd
RA0,DE0=156.96346,-43.90461
POS={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15)}
F={'Xb47':f('Xb47','spw16_18_20_22',kind='cont.I.pbcor'),'Xb4f':f('Xb4f','spw16_18_20_22',kind='cont.I.pbcor'),'Xb4d':f('Xb4d','spw25_27_29_31',kind='cont.I.pbcor'),'Xb4b':f('Xb4b','spw25_27_29_31',kind='cont.I.pbcor'),'Xb55':f('Xb55','spw16_18_20_22',kind='cont.I.pbcor')}
from scipy.ndimage import map_coordinates
OM=lambda b: np.pi*b[0]*b[1]/(4*np.log(2))
def img(o,tb,crop=14):
    h=fits.open(F[o],memmap=True); hd=h[0].header; a=np.squeeze(h[0].data).astype(float); w=WCS(hd).celestial
    nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); pix=abs(hd['CDELT2'])*3600
    cx,cy=[float(x) for x in w.world_to_pixel_values(RA0,DE0)]; hp=int(crop/pix); x0,y0=int(cx)-hp,int(cy)-hp
    sub=np.nan_to_num(a[max(y0,0):int(cy)+hp+1,max(x0,0):int(cx)+hp+1]); x0=max(x0,0); y0=max(y0,0)
    k=kernel(hd,tb,nb)
    if k is None: return None
    s=fftconvolve(sub,k,mode='same')*OM(tb)/OM(nb)
    out={}
    for r,(x,y) in POS.items():
        px,py=[float(q) for q in w.world_to_pixel_values(RA0+x/3600/np.cos(np.deg2rad(DE0)),DE0+y/3600)]
        out[r]=float(map_coordinates(s,[[py-y0],[px-x0]],order=1)[0])*1e3
    return out,nb,hd.get('RESTFRQ',None)
rows=[]
for bn,tb in {'Xb4f(6.2x3.7,82)':(6.2,3.7,82.),'Xb55(9.5x6.4,-88)':(9.5,6.4,-88.)}.items():
    res={o:img(o,tb) for o in F}
    for o,r in res.items():
        if r is None: continue
        for reg,v in r[0].items(): rows.append(dict(beam=bn,array=o,region=reg,mJy_beam=v))
R=pd.DataFrame(rows); R.to_csv('/workspace/out/q8/continuum_matched.csv',index=False)
P=R.pivot_table(index=['beam','array'],columns='region',values='mJy_beam'); print(P.round(3))
for bn in R.beam.unique():
    p=P.loc[bn]
    for a in p.index:
        if a=='Xb4f' or 'Xb4f' not in p.index: continue
        print(bn,a,'/Xb4f',(p.loc[a]/p.loc['Xb4f']).round(3).to_dict(),' N/S double ratio',round((p.loc[a,'N']/p.loc['Xb4f','N'])/(p.loc[a,'S']/p.loc['Xb4f','S']),3))
