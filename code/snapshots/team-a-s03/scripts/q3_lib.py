# Q3 library: Xb4d/Xb4b cubes cropped around the nuclei, optionally convolved channel-by-channel to a common circular beam; empirical per-pixel noise of window sums from line-free channels.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
from scipy.ndimage import uniform_filter1d
VS=2775.   # fiducial systemic velocity (km/s, radio LSRK): Sakamoto+2014 Vsys=2775; S-local line centre is ~2722-2740 (p0001)
KEY={'HCN':'HCN 3-2@%s','HCO+':'HCO+ 3-2@%s'}
# line-free velocity ranges (radio km/s, in each line's own frame): exclude Vsys+-650, HOC+ 3-2 (-1000 km/s of HCO+), HCN v2=1f (+401 of HCO+), band edges
FREE={'HCN':[(1380,2100)],'HCO+':[(1280,1680),(1880,2100)]}
def cube(sp,ous='Xb4d',crop=14.,fw=None):
    ous_,spw,f0=LINES[KEY[sp]%ous]; h=fits.open(f(ous_,spw),memmap=True); hd=h[0].header
    fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C; dv=abs(np.median(np.diff(v))); pix=abs(hd['CDELT2'])*3600
    cx,cy=[float(a) for a in WCS(hd).celestial.world_to_pixel_values(RA0,DE0)]; hp=int(crop/pix)
    x0=int(cx)-hp; y0=int(cy)-hp; x1=int(cx)+hp+1; y1=int(cy)+hp+1
    a=np.nan_to_num(np.array(h[0].data[0,:,y0:y1,x0:x1],np.float32))
    nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
    if fw:
        k=kernel(hd,(fw,fw,0.),nb).astype(np.float32); sc=OM((fw,fw,0))/OM(nb)
        a=np.array([fftconvolve(c,k,mode='same')*sc for c in a],np.float32)
    dx,dy=offsets(hd,y1-y0,x1-x0,x0,y0)
    return dict(v=v,dv=dv,a=a,dx=dx,dy=dy,pix=pix,hd=hd,x0=x0,y0=y0,beam=(fw,fw,0) if fw else nb,f0=f0,sp=sp)
def wsum(c,lo,hi):
    m=(c['v']>=lo)&(c['v']<hi); return c['a'][m].sum(0)*c['dv'],int(m.sum())
def wnoise(c,lo,hi):
    """per-pixel sigma (Jy/beam km/s) of a window sum of the same channel count, from sliding sums over line-free channels (robust MAD)"""
    n=int(((c['v']>=lo)&(c['v']<hi)).sum()); v=c['v']; S=[]
    cs=np.cumsum(np.concatenate([np.zeros((1,)+c['a'].shape[1:],np.float32),c['a']]),0)
    for i in range(len(v)-n):
        if any(v[i]>=a and v[i+n-1]<=b for a,b in FREE[c['sp']]): S.append(cs[i+n]-cs[i])
    S=np.array(S)*c['dv']; med=np.median(S,0)
    return 1.4826*np.median(np.abs(S-med),0),len(S)
def samp(c,img,x,y):
    from scipy.ndimage import map_coordinates
    w=WCS(c['hd']).celestial
    px,py=[float(a) for a in w.world_to_pixel_values(RA0+x/3600/np.cos(np.deg2rad(DE0)),DE0+y/3600)]
    return float(map_coordinates(img,[[py-c['y0']],[px-c['x0']]],order=1,mode='nearest')[0])
def spec(c,x,y,rad=0.):
    """spectrum (Jy/beam) at offset x,y; rad>0: mean in circle (arcsec)"""
    if rad<=0:
        from scipy.ndimage import map_coordinates
        w=WCS(c['hd']).celestial
        px,py=[float(a) for a in w.world_to_pixel_values(RA0+x/3600/np.cos(np.deg2rad(DE0)),DE0+y/3600)]
        return np.array([map_coordinates(ch,[[py-c['y0']],[px-c['x0']]],order=1,mode='nearest')[0] for ch in c['a']])
    m=np.hypot(c['dx']-x,c['dy']-y)<=rad; return c['a'][:,m].mean(1)
def chan_noise(c):
    m=np.zeros(len(c['v']),bool)
    for a,b in FREE[c['sp']]: m|=(c['v']>=a)&(c['v']<=b)
    return m
