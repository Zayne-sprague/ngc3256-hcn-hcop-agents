# Q25: curves of growth (cumulative flux inside r, 10" beam) of HCN/HCO+ 1-0 (Harada, Michiyama) and 2-1 (Xb55), and Fourier transfer ratio Harada/Michiyama.
import sys; sys.path.insert(0,'/workspace/scripts')
from q25_ratio10 import *
from scipy.ndimage import map_coordinates
RS=[4,8,12,16,20,25,30]
G={}
for nm,p in P.items():
    for ln,f0 in REST.items():
        d=load10(p,f0); s=sm(d); r=np.hypot(d['dx']+0.075,d['dy']); omp=OM(tb)/d['pix']**2
        G[(nm,ln)]=[np.nansum(s[r<R])/omp for R in RS]; G[(nm,ln,'map')]=(d,s)
for k in ('HCN 2-1','HCO+ 2-1'):
    d=load(k,nullspec=False); s=smooth(d,10.); r=np.hypot(d['dx']+0.075,d['dy']); omp=OM(tb)/d['pix']**2
    G[('Xb55',k)]=[np.nansum(s[r<R])/omp for R in RS]
print('cumulative flux Jy km/s inside r =',RS)
for k,v in G.items():
    if len(k)==2: print(f'{k[0]:9s} {k[1]:9s}',np.round(v,2))
print('Harada/Michiyama cumulative: HCN',np.round(np.array(G[("Harada","HCN 1-0")])/np.array(G[("Michiyama","HCN 1-0")]),3),' HCO+',np.round(np.array(G[("Harada","HCO+ 1-0")])/np.array(G[("Michiyama","HCO+ 1-0")]),3))
for nm in ('Harada','Michiyama','Xb55'):
    h=np.array(G[(nm,'HCN 1-0' if nm!='Xb55' else 'HCN 2-1')]);o=np.array(G[(nm,'HCO+ 1-0' if nm!='Xb55' else 'HCO+ 2-1')]); print(nm,'cumulative HCN/HCO+',np.round(h/o,3))
# Fourier amplitude ratio on common grid (0.25" pixels, 128" box), in uv annuli
def regrid(d,s,n=512,pix=0.25):
    x=(np.arange(n)-n/2)*pix; X,Y=np.meshgrid(x,x)   # arcsec offsets from RA0,DE0 (x east)
    w=WCS(d['hd']).celestial; ra=RA0+(X)/3600./np.cos(np.deg2rad(DE0)); de=DE0+Y/3600.
    px,py=w.world_to_pixel_values(ra,de); return map_coordinates(np.nan_to_num(s),[py,px],order=1,cval=0.)
print('\nFourier amplitude ratio Harada/Michiyama of the 10" maps (relative to u<2 kl); u in kl (angular)')
UB=[(0,1),(1,2),(2,4),(4,6),(6,9),(9,13)]
for ln in REST:
    A=[regrid(*G[(nm,ln,'map')]) for nm in ('Harada','Michiyama')]
    Fa=[np.abs(np.fft.fftshift(np.fft.fft2(a))) for a in A]; u=np.fft.fftshift(np.fft.fftfreq(512,d=0.25/206265.))/1e3; U=np.hypot(*np.meshgrid(u,u))
    print(ln,[ (b,round(float(np.mean(Fa[0][(U>=b[0])&(U<b[1])])/np.mean(Fa[1][(U>=b[0])&(U<b[1])])),3)) for b in UB])
