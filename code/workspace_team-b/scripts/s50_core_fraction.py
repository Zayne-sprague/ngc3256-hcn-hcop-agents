# Q4 round 2, fix 4: propagate the N-core flux fraction (not a single value) into the H13CN-derived opacity bound, with a two-component (core + envelope) ladder model.
# f3 = fraction of the HCN 3-2 flux of the Xb47 beam aperture (r<1.7") inside r<0.5" (Xb4b 0.25" annuli, s43/s45). Xb4b resolves out extended flux (Xb4d/Xb4b recovers 65%/47% of Xb4f, s26), so the same core flux is also referred to the
# Xb4d (0.95") total in r<1.7" to give a lower f3.  f4 = core fraction of the HCN 4-3 flux: f4 = f3 * L_c / L_tot with L = Tb(4-3)/Tb(3-2) (RJ, K km/s): L_tot measured at the common 6.5" beam (Xb47 4-3 vs Xb4d/Xb4f 3-2; the only matched-beam ladder), L_c from RADEX (thick 0.85-1.07; mid 0.9; thin 0.53-0.79; s44).
# Then bound on the core: f4/r_c + (1-f4)/X_env <= 1/L_H13  (L_H13 = 8.3-9.5 from s41), r_c = HCN/H13CN of the core = (1-e^-tau)/(1-e^-tau/X_c); report tau_max(core) as a band over f3, L_c, X.
import sys,warnings,json; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from common import *
from regions import REGIONS
from scipy.optimize import brentq
# (1) f3 from independent totals
b=pd.read_csv('/workspace/out/s43_profiles.csv'); n=b[b.region=='N'].reset_index(drop=True); ps=0.047; bm=1.1331*0.2456*0.2373/ps**2
fl=n.I_HCN*n.npix/bm; core=fl[:2].sum(); tot_b=fl[:3].sum()+fl[3]*(1.7**2-1.0**2)/(2.0**2-1.0**2)
print('Xb4b HCN 3-2 flux (Jy km/s): r<0.5 core %.2f; r<1.7 total (Xb4b) %.2f -> f3(Xb4b)=%.2f'%(core,tot_b,core/tot_b))
def aperture_flux(ous,spw,f0,r,half=None):
    h=fits.open(D+f'member.uid___A001_X1273_{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits',memmap=True); hd=h[0].header
    nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])+1-hd['CRPIX3'])*hd['CDELT3']; v=C*(1-nu/f0)-VSYS; ok=np.where(abs(v)<350)[0]
    psx=abs(hd['CDELT1'])*3600; x,y=WCS(hd).celestial.world_to_pixel_values(REGIONS['N'][0],REGIONS['N'][1]); hw=int(r/psx)+6; x0,y0=int(round(float(x))),int(round(float(y)))
    cu=np.array(h[0].data[0,ok[0]:ok[-1]+1,y0-hw:y0+hw+1,x0-hw:x0+hw+1],dtype=np.float32); ref=np.array(h[0].data[0,:40,y0-hw:y0+hw+1,x0-hw:x0+hw+1],dtype=np.float32)
    dv=abs(np.diff(v)).mean(); yy,xx=np.mgrid[-hw:hw+1,-hw:hw+1]; rr=np.hypot(xx-(float(x)-x0),yy-(float(y)-y0))*psx
    I=np.nansum(cu,0)*dv; base=np.nanmedian(ref.sum(0)*0)  # cubes are continuum subtracted
    bmaj,bmin=hd['BMAJ']*3600,hd['BMIN']*3600; A=1.1331*bmaj*bmin/psx**2
    return np.nansum(I[rr<r])/A
for ous,spw,lab in (('Xb4d','spw27','Xb4d 0.94x0.74"'),('Xb4f','spw18','Xb4f 5.7x3.5"')):
    for r in (1.7,2.5,3.5):
        print(f'  HCN 3-2 flux r<{r}" {lab}: {aperture_flux(ous,spw,265.886434e9,r):.2f} Jy km/s')
Fd=aperture_flux('Xb4d','spw27',265.886434e9,1.7); f3d=core/Fd; print(f'-> f3 referred to Xb4d total r<1.7" = {f3d:.2f} (Xb4b-only: {core/tot_b:.2f}); a core flux of 5.8 Jy km/s is {core/aperture_flux("Xb4f","spw18",265.886434e9,2.5):.2f} of the Xb4f r<2.5" flux')
F3={'Xb4b-only':core/tot_b,'Xb4d-total':f3d}
# (2) L_tot from the matched 6.5" beam (s30): Jy/beam km/s, N 'pos'
S43={'HCN':17.941,'HCO+':42.261}; S32={'HCN':20.362,'HCO+':37.438}   # Xb47 J4 6.5as, Xb4d J3 6.5as
Ltot={m:S43[m]/S32[m]*(FJ32/FJ43)**2 for m,FJ32,FJ43 in (('HCN',265.886434,354.505478),('HCO+',267.557626,356.734240))}
print('L_tot = Tb(4-3)/Tb(3-2) at 6.5" (N, pos):',{k:round(v,3) for k,v in Ltot.items()},'; flux-scale (Band 6 vs 7) systematic ~15% on this ratio')
Lt=Ltot['HCN']
# (3) bound
def tau_from_r(r,X):  # tau such that (1-e^-t)/(1-e^-t/X)=r, r<X
    g=lambda t:(1-np.exp(-t))/(1-np.exp(-t/X))-r
    return brentq(g,1e-3,1e5) if (g(1e-3)>0 and r>1.0) else (0. if r>=X else np.nan)
rows=[]
for f3 in (0.15,0.20,0.25,0.30,0.36,0.40,0.50):
    for Lc,lab in ((1.0,'thick L_c=1.0'),(0.9,'mid 0.9'),(0.57,'thin 0.57')):
        for Lt_ in (Lt*0.85,Lt,Lt*1.15):
            f4=min(f3*Lc/Lt_,1.0)
            for X in (30,50,100):
                for LH in (8.3,9.5):
                    den=1/LH-(1-f4)/X
                    rmin=f4/den if den>0 else np.inf
                    tmax=tau_from_r(rmin,X) if np.isfinite(rmin) else np.nan
                    rows.append(dict(f3=f3,Lc=lab,Ltot=round(Lt_,3),f4=f4,X=X,LH13=LH,r_core_min=rmin,tau_core_max=tmax,thick_ok=(tmax>=10) if np.isfinite(tmax) else False))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q4r2/s50_core_bound.csv',index=False)
pd.set_option('display.width',250)
print('\nf4 = core fraction of the HCN 4-3 flux (central L_tot):'); print(df[(df.Ltot==round(Lt,3))].pivot_table(index='f3',columns='Lc',values='f4').round(2).to_string())
print('\ntau_core(HCN 4-3) upper limit [median over LH13 8.3/9.5; central L_tot]; rows f3, cols (core type, X):')
c=df[df.Ltot==round(Lt,3)].pivot_table(index='f3',columns=['Lc','X'],values='tau_core_max',aggfunc='median'); print(c.round(1).to_string())
print('\nwhole band (f3 0.15-0.50, Ltot +-15%, L_c types, X 30-100, LH13 8.3-9.5): tau_core,max range & fraction of scenarios in which a tau>=10 core is allowed')
for lab,g in df.groupby('Lc'): print(f'  {lab}: tau_max min {g.tau_core_max.min():.1f}  median {g.tau_core_max.median():.1f}  max {g.tau_core_max.max():.1f}; tau>=10 allowed in {100*g.thick_ok.mean():.0f}% of scenarios')
for X,g in df.groupby('X'): print(f'  X={X}: tau_max median {g.tau_core_max.median():.1f}; thick(>=10) allowed {100*g.thick_ok.mean():.0f}%')
sel=df[(df.f3.isin([0.20,0.36,0.50]))&(df.Ltot==round(Lt,3))&(df.LH13==9.5)]
print('\nthick-core (L_c=1.0) tau_max by f3 and X at LH13=9.5 (central):'); print(sel[sel.Lc.str.startswith('thick')].pivot_table(index='f3',columns='X',values='tau_core_max').round(1).to_string())
print('and f3 referred to Xb4d total (%.2f), Xb4b-only (%.2f)'%(f3d,core/tot_b))
json.dump(dict(F3=F3,Ltot=Ltot),open('/workspace/out/q4r2/s50_core_fraction.json','w'),indent=1)
