# Q15c: velocity structure of the red wing: aperture spectra, per-velocity-bin centroids (does the ridge track N, or sit ~1.8" S of N?), p-v cut along the S-N axis. Xb4d 1.05" beam.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from scipy.optimize import curve_fit
from q3_pos import POS
O='/workspace/out/q15/'; log=open(O+'q15c.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
z=np.load(O+'cube_Xb4d.npz'); dx,dy,dv=z['dx'],z['dy'],float(z['dv'])
Sx,Sy=POS['S']; Nx,Ny=POS['N']
def binspec(cube,v,x,y,rad,edges):
    m=np.hypot(dx-x,dy-y)<=rad; s=cube[:,m].mean(1); out=[]
    for a,b in zip(edges[:-1],edges[1:]):
        k=(v>=a)&(v<b); out.append(s[k].mean()*1e3)
    return np.array(out)
edges=np.arange(2380,3281,40.)
cen=0.5*(edges[1:]+edges[:-1])
POSI={'S':(Sx,Sy),'N':(Nx,Ny),'ridge(S+3.3)':(Sx,Sy+3.3),'mid N-S (Bridge)':(0.15,0.0),'ridge(N-1.0)':(0.0,Ny-1.0),'lobe(S-1.93)':POS['Soutflow']}
for sp,cube,v in (('HCN',z['HCN'],z['vH']),('HCO+',z['HCO'],z['vO'])):
    pr(f'## {sp} Xb4d: mean spectrum within r=0.5" (mJy/beam), 40 km/s bins; per-channel rms (3.9 km/s ch) ~{np.median(z["rmsH" if sp=="HCN" else "rmsO"])*1e3:.2f} mJy/b; 40-km/s bin noise ~{np.median(z["rmsH" if sp=="HCN" else "rmsO"])*1e3/np.sqrt(40/dv/1.0):.2f} (formal, correlated pixels ignored)')
    pr('v_bin '+' '.join(f'{c:6.0f}' for c in cen))
    for nm,(x,y) in POSI.items(): pr(f'{nm:18s}'+' '.join(f'{a:6.1f}' for a in binspec(cube,v,x,y,0.5,edges)))
# centroid fits per velocity bin: beam-fixed 2D Gaussian (1.05" FWHM) on bin map (sum of channels -> Jy/b km/s), within 4" of N
def g2(xy,A,x0,y0,b):
    x,y=xy; s=1.05/2.3548; return A*np.exp(-((x-x0)**2+(y-y0)**2)/(2*s**2))+b
def vbin_map(cube,v,a,b): m=(v>=a)&(v<b); return cube[m].sum(0)*dv
pr('## per-velocity-bin point-source (beam) fits within 4.0" of the N nucleus: centroid dN relative to N (positive = north) and amplitude (Jy/b km/s); two components allowed? single beam-sized component, free background')
for sp,cube,v in (('HCN',z['HCN'],z['vH']),('HCO+',z['HCO'],z['vO'])):
    for a in range(2840,3290,60):
        M=vbin_map(cube,v,a,a+60); m=np.hypot(dx-Nx,dy-(Ny-1.5))<=3.5
        best=None
        for y0s in (Ny,Ny-1.8):
            try:
                p,cv=curve_fit(g2,(dx[m],dy[m]),M[m],p0=[0.1,Nx,y0s,0],bounds=([0,-3,Ny-4.5,-0.2],[3,3,Ny+1.5,0.2]))
                r=((M[m]-g2((dx[m],dy[m]),*p))**2).sum()
                if best is None or r<best[0]: best=(r,p,np.sqrt(np.diag(cv)))
            except Exception as e: pass
        if best:
            p=best[1]; e=best[2]*np.sqrt(m.sum()/(m.sum()/(1.1331*1.05**2/0.15**2)))
            pr(f'{sp:5s} v {a:4d}-{a+60:4d} ({a-2730:+4d}..{a+60-2730:+4d}): dE={p[1]:+.2f} dN(from N)={p[2]-Ny:+.2f} (+-{e[2]:.2f}) dN(from S)={p[2]-Sy:+.2f} amp {p[0]:.3f}')
# p-v along the S-N axis: strip |dE-0|<0.6, rows dN each 0.5"
pr('## p-v cut along x in [-0.6,+0.6] (all dE rel phase centre), HCN, mJy/beam per 40-km/s bin; row = offset from S along N (positive north)')
for sp,cube,v in (('HCN',z['HCN'],z['vH']),('HCO+',z['HCO'],z['vO'])):
    pr(f'## {sp}  v_bin '+' '.join(f'{c:5.0f}' for c in cen))
    for yo in np.arange(-5.0,8.1,0.75):
        y=Sy+yo; m=(np.abs(dx-0.0)<=0.6)&(np.abs(dy-y)<=0.375)
        s=cube[:,m].mean(1); row=[s[(v>=a)&(v<b)].mean()*1e3 for a,b in zip(edges[:-1],edges[1:])]
        pr(f'dN_S={yo:+5.2f} (y={y:+.2f}) '+' '.join(f'{a:5.0f}' for a in row))
