# Q3 round-2 library: Band-6 cube access (frequency axis), circular-beam smoothing, aperture spectra, baselines, bootstrap-free noise from line-free channels.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q3_pos import *
OUT='/workspace/out/q3b/'
F0={'HCN':265.886434,'HCO+':267.557626}
SPW={'HCN':'spw27','HCO+':'spw25'}
def fcube(ous,spw,fw=1.05,crop=9.,kind='cube.I.pbcor',sep='.',center=(0.,0.)):
    """crop (arcsec half-size about RA0,DE0 + center) cube, convolved channel-wise to a circular beam fw (None = native). Returns dict with freq (GHz), a[ch,y,x] (Jy/beam of the smoothed beam), dx,dy, pix, hd, x0,y0"""
    h=fits.open(f(ous,spw,kind,sep),memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9; pix=abs(hd['CDELT2'])*3600
    cx,cy=[float(a) for a in WCS(hd).celestial.world_to_pixel_values(RA0,DE0)]; hp=int(crop/pix)
    x0=max(0,int(cx)-hp); x1=int(cx)+hp+1; y0=max(0,int(cy)-hp); y1=int(cy)+hp+1
    a=np.nan_to_num(np.array(h[0].data[0,:,y0:y1,x0:x1] if h[0].data.ndim==4 else h[0].data[:,y0:y1,x0:x1],np.float32))
    nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); beam=nb
    if fw:
        k=kernel(hd,(fw,fw,0.),nb)
        if k is not None:
            k=k.astype(np.float32); sc=OM((fw,fw,0))/OM(nb); a=np.array([fftconvolve(c,k,mode='same')*sc for c in a],np.float32); beam=(fw,fw,0)
    dx,dy=offsets(hd,y1-y0,x1-x0,x0,y0)
    return dict(fr=fr,a=a,dx=dx,dy=dy,pix=pix,hd=hd,x0=x0,y0=y0,beam=beam,ous=ous,spw=spw)
def aper(c,x,y,rad):
    """mean spectrum within circle rad (arcsec; <=0 -> bilinear point sample) at offset (x,y)"""
    if rad>0:
        m=np.hypot(c['dx']-x,c['dy']-y)<=rad; return c['a'][:,m].mean(1)
    from scipy.ndimage import map_coordinates
    w=WCS(c['hd']).celestial; px,py=[float(a) for a in w.world_to_pixel_values(RA0+x/3600/np.cos(np.deg2rad(DE0)),DE0+y/3600)]
    return np.array([map_coordinates(ch,[[py-c['y0']],[px-c['x0']]],order=1,mode='nearest')[0] for ch in c['a']])
def vel(fr,f0): return (1-fr/f0)*C
def robust_sigma(x): x=x[np.isfinite(x)]; return 1.4826*np.median(np.abs(x-np.median(x)))
def baseline(v,s,free,order=1):
    """subtract poly baseline fitted to channels in the velocity ranges `free`; returns (s_corr, per-channel sigma from residual)"""
    m=np.zeros(len(v),bool)
    for a,b in free: m|=(v>=a)&(v<=b)
    p=np.polyfit(v[m]-v[m].mean(),s[m],order); b=np.polyval(p,v-v[m].mean()); r=s-b
    return r,robust_sigma(r[m])
def rebin(v,s,edges):
    """mean of s in velocity bins (edges), returns centres, means, counts"""
    cen=[];mm=[];nn=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        k=(v>=lo)&(v<hi); cen.append(0.5*(lo+hi)); mm.append(s[k].mean() if k.any() else np.nan); nn.append(k.sum())
    return np.array(cen),np.array(mm),np.array(nn)
def vcube(sp,ous='Xb4d',fw=1.05,crop=9.):
    """HCN or HCO+ 3-2 cube with velocity axis (radio, own rest frequency)"""
    kind,sep=('cube.I.pbcor','.')
    c=fcube(ous,SPW[sp],fw,crop,kind,sep); c['v']=vel(c['fr'],F0[sp]); c['dv']=abs(np.median(np.diff(c['v']))); c['sp']=sp
    c['free']=FREE[sp]; return c
