# Q22: (a) spw16/spw18 overlap (sky 176.366-176.469 GHz) and spw20/spw22 overlap (164.536-165.110): same sky frequencies in two spws -> is the ramp a property of the sky/continuum subtraction or of the spw edge?
# (b) exponential edge-ramp fit in spw22 (low-freq edge) vs region continuum.
import numpy as np, json, warnings, logging, glob
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import curve_fit
L=lambda s: np.load(f'/workspace/out/regspec/Xb55.spw{s}.cube.pbcor.npz')
z={s:L(s) for s in (16,18,20,22)}; nu={s:z[s]['nu']/1e9 if z[s]['nu'].max()>1e6 else z[s]['nu'] for s in z}
regs=['N','S','Bridge','E1','E2','W','Sblue']
h=fits.open(glob.glob('/data/member.uid___A001_X1273_Xb55.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits')[0]); d=np.squeeze(h[0].data); w=WCS(h[0].header).celestial
cont={}
for r in regs:
    x,y=w.world_to_pixel_values(*z[22][r+'_xy']); cont[r]=float(d[int(round(float(y))),int(round(float(x)))])*1e3
out={'cont':cont,'overlap':{},'fit':{}}
def ov(a,b,name):
    lo=max(nu[a][0],nu[b][0]); hi=min(nu[a][-1],nu[b][-1])
    print(f'\noverlap {name}: {lo:.4f}-{hi:.4f} GHz')
    # bin by 10 channels in overlap in common frequency grid
    fgrid=np.arange(lo,hi,0.01953)
    for r in regs[:3]+['W']:
        ya=np.interp(fgrid,nu[a],z[a][r+'_pix']*1e3); yb=np.interp(fgrid,nu[b],z[b][r+'_pix']*1e3)
        bs=10 if len(fgrid)>=30 else 5; nb=len(fgrid)//bs*bs; A=ya[:nb].reshape(-1,bs).mean(1); B=yb[:nb].reshape(-1,bs).mean(1)
        print(f'  {r:6s} spw{a}-binned',np.round(A,1)[:8],'\n         spw{b}-binned'.replace('{b}',str(b)),np.round(B,1)[:8],' mean diff(b-a) %.2f +- %.2f'%((B-A).mean(),(B-A).std()/np.sqrt(len(A))))
        out['overlap'][f'{name}.{r}']=dict(a=A.tolist(),b=B.tolist())
ov(16,18,'16/18'); ov(20,22,'20/22')
# (b)
f=lambda x,A,t,c: A*np.exp(-x/t)+c
print('\nspw22 low-edge exponential fit y=A exp(-ch/tau)+c, ch=0..110, noise 2.0 mJy/ch, 1.3x corr: ')
ch=np.arange(110.)
for r in regs:
    y=z[22][r+'_pix'][:110]*1e3
    try:
        p,c=curve_fit(f,ch,y,p0=[-3,30,0],bounds=([-50,3,-5],[0.5,300,5]),sigma=np.full(110,2.1*1.3),absolute_sigma=True,maxfev=5000); e=np.sqrt(np.diag(c))
    except Exception as ex: print(r,ex); continue
    print(f'  {r:6s} cont {cont[r]:6.2f}  A={p[0]:6.2f}+-{e[0]:.2f}  tau={p[1]:6.1f}+-{e[1]:.1f} ch ({p[1]*1.9534:.0f} MHz)  c={p[2]:5.2f}+-{e[2]:.2f}  A/cont={p[0]/cont[r]:.3f}')
    out['fit'][r]=dict(A=float(p[0]),eA=float(e[0]),tau=float(p[1]),etau=float(e[1]),c=float(p[2]),ec=float(e[2]),cont=cont[r])
print('\nspw22 10-ch binned profile from low edge (mJy/b):')
for r in regs: print(f'  {r:6s}',np.round(z[22][r+'_pix'][:200].reshape(20,10).mean(1)*1e3,1))
json.dump(out,open('/workspace/out/q22/s132.json','w'))
