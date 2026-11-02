# Q7: disk 2-1 HCN/HCO+ in annuli about the N-S midpoint (Xb55, smoothed to 10"), zero level variants. Independent of team-b s36 (own code, q1_lib moment maps).
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
H=load('HCN 2-1',nullspec=False); O=load('HCO+ 2-1',nullspec=False)
tb=10.
sH,sO=smooth(H,tb),smooth(O,tb)
xm,ym=0.0,0.0; r=np.hypot(H['dx']-xm,H['dy']-ym); pix=H['pix']
omp=OM((tb,tb,0))/pix**2
print('map half-size',r.max()/1.41,'arcsec; pixel',pix)
def zl(far):
    m=(r>=far[0])&(r<far[1])&(r<40); return np.nanmean(sH[m]),np.nanmean(sO[m]),m.sum()
for far in [(30,45),(26,40),(None,None)]:
    z=zl(far) if far[0] else (0,0,0)
    print(f'## zero level: far-field {far} -> HCN {z[0]:.3f} HCO+ {z[1]:.3f} Jy/b km/s (n={z[2]})')
    for a,b in [(0,4),(4,8),(8,12),(12,17)]:
        m=(r>=a)&(r<b); nH=m.sum()
        fH=(sH[m]-z[0]).sum()/omp; fO=(sO[m]-z[1]).sum()/omp
        # zero-level error propagated: far-field mean uncertainty ~ pixel rms/sqrt(n_eff beams)
        print(f'  {a:2d}-{b:2d}": HCN {fH:7.2f} HCO+ {fO:7.2f} Jy km/s  R={fH/fO:.3f}')
# sensitivity to zero level: shift by +-1 sigma of the far-field mean (use pixel rms/sqrt(n beams in far field))
m=(r>=30)&(r<45); sd=[np.nanstd(sH[m]),np.nanstd(sO[m])]; nb=m.sum()/omp; ez=[s/np.sqrt(nb) for s in sd]
print('far-field (30-45") pixel rms',sd,' n_beams',nb,' zero-level 1-sigma',ez)
z=zl((30,45)); a,b=8,12; mm=(r>=a)&(r<b)
for dz in (-1,0,1):
    fH=(sH[mm]-z[0]-dz*ez[0]).sum()/omp; fO=(sO[mm]-z[1]-dz*ez[1]).sum()/omp; print(f'  8-12" zero level shifted {dz:+d} sigma (HCN & HCO+ together): R={fH/fO:.3f}')
